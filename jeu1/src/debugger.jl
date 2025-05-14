include("io.jl")
include("solver2.jl")

using JSON
using JuMP
using CPLEX

function build_model(game::PalisadeGame)
    h, w = game.height, game.width
    k = game.region_size
    V = [(i, j) for i in 1:h, j in 1:w]  # 所有格子
    E = []

    # 构造边（只考虑右和下，避免重复）
    for i in 1:h
        for j in 1:w
            if j < w
                push!(E, ((i,j), (i,j+1)))
            end
            if i < h
                push!(E, ((i,j), (i+1,j)))
            end
        end
    end

    R = collect(1:(h*w)÷k)  # 区域编号

    # 开始建模
    m = Model(CPLEX.Optimizer)

    # ========= 变量 =========

    @variable(m, wvars[e in E], Bin, base_name = "w")  # 墙体变量（边 e = ((i,j),(i',j'))）
    @variable(m, x[i=1:h, j=1:w, r in R], Bin, base_name = "x")  # 区域归属变量
    @variable(m, y[i=1:h, j=1:w, r in R], Bin, base_name = "y")  # 区域根节点变量
    @variable(m, f[i=1:h, j=1:w, i2=1:h, j2=1:w, r in R] >= 0, base_name = "f")  # 连通性流变量（可取负）

    # ========= 约束 =========

    # 每格只能属于一个区域
    @constraint(m, [i=1:h, j=1:w], sum(x[i,j,r] for r in R) == 1, base_name="one_region")

    # 每个区域恰好包含 k 个格子
    @constraint(m, [r in R], sum(x[i,j,r] for (i,j) in V) == k, base_name="region_size")

    # 每个区域有一个根
    @constraint(m, [r in R], sum(y[i,j,r] for (i,j) in V) == 1, base_name="region_root_unique")
    @constraint(m, [i=1:h, j=1:w, r in R], y[i,j,r] <= x[i,j,r], base_name="root_in_region")

    # 相邻格子若不在同一区域，则存在墙
    for ((i1,j1), (i2,j2)) in E, r in R
        key = ((i1,j1),(i2,j2))
        @constraint(m, x[i1,j1,r] - x[i2,j2,r] <= wvars[key], base_name="wall_diff_pos")
        @constraint(m, x[i2,j2,r] - x[i1,j1,r] <= wvars[key], base_name="wall_diff_neg")
    end

    for ((i1,j1),(i2,j2)) in E
        @constraint(m, sum(x[i1,j1,r] + x[i2,j2,r] for r in R) <= 1 + (1 - wvars[minmax((i1,j1),(i2,j2))]),
                    base_name="internal_wall_link_$(i1)$(j1)$(i2)$(j2)")
    end 
    
    # clue墙数约束（修复边界墙恒为1的情况）
    for ((i,j), clue) in game.clues
        terms = Any[]
        if i > 1
            push!(terms, wvars[minmax((i,j), (i-1,j))])
        else
            push!(terms, 1)
        end
        if i < h
            push!(terms, wvars[minmax((i,j), (i+1,j))])
        else
            push!(terms, 1)
        end
        if j > 1
            push!(terms, wvars[minmax((i,j), (i,j-1))])
        else
            push!(terms, 1)
        end
        if j < w
            push!(terms, wvars[minmax((i,j), (i,j+1))])
        else
            push!(terms, 1)
        end
        @constraint(m, sum(terms) == clue, base_name="clue")
    end

    # 连通性流约束：流平衡（根点出流为k-1，其他点净流=-1）
    for (i,j) in V, r in R
        neighbors = [(i+di, j+dj) for (di,dj) in ((0,1), (1,0), (0,-1), (-1,0)) if 1 <= i+di <= h && 1 <= j+dj <= w]
        inflow = sum(f[i2,j2,i,j,r] for (i2,j2) in neighbors)
        outflow = sum(f[i,j,i2,j2,r] for (i2,j2) in neighbors)
        @constraint(m, outflow - inflow == k*y[i,j,r] - x[i,j,r], base_name="flow_balance_$i$j$r")
    end

    # 流只能通过没有墙、且流过的格子必须属于对应区域
    for ((i1,j1),(i2,j2)) in E, r in R
        key = ((i1,j1),(i2,j2))
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*(1 - wvars[key]), base_name="flow_wall")
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*x[i1,j1,r], base_name="flow_xfrom")
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*x[i2,j2,r], base_name="flow_xto")
    end

    # ========= 目标函数 =========
    @objective(m, Min, 0)

    return m
end

sol_path = "palisade_feasible_solution.json"

"""
check_solution_conflicts(model::Model, json_path::String)
读取 JSON 格式的解，并检查模型中不满足的约束。

返回：
- 一个 Vector，包含违反约束的 (name, residual) 元组
"""
function check_solution_conflicts(model::Model, json_path::String)
    println("载入解并检查约束...")

    # === 读取 JSON 解 ===
    solution = JSON.parsefile(json_path)

    # === 设置变量值 ===
    function set_value_wall(dict, var_prefix)
        for (k_str, v) in dict
            k = eval(Meta.parse(k_str))  # 转换为元组 ((1,1),(2,1))
            set_start_value(getindex(model[var_prefix], k), v)
        end
    end    

    function set_value(dict, var_prefix)
        for (k_str, v) in dict
            # 解析字符串为整数索引（不构造 Tuple）
            indices = parse.(Int, split(strip(k_str, ['(', ')', ' ']), ','))  # 去掉括号和空格
            try
                set_start_value(getindex(model[var_prefix], indices...), v)
            catch e
                @warn "变量 $var_prefix[$(join(indices,','))] 不存在或赋值失败：$e"
            end
        end
    end

    function fill_missing_zero!(model::Model, var_prefix::Symbol, assigned_dict::Dict)
        var_container = model[var_prefix]
        total = 0
        filled = 0
        for idx in keys(var_container)
            total += 1
            strkey = string(idx)  # 直接转成 "(i,j,...)"
            if !haskey(assigned_dict, strkey)
                set_start_value(var_container[idx], 0.0)
                filled += 1
            end
        end
        println("变量 $var_prefix 缺省值补 0 完成：$filled / $total")
    end    

    set_value_wall(solution["wvars"], :wvars)
    set_value(solution["xvars"], :x)
    set_value(solution["yvars"], :y)
    fill_missing_zero!(model, :f, solution["fvars"])
    set_value(solution["fvars"], :f)
    

    # === 检查所有约束 ===
    function eval_affine_function(func::MOI.ScalarAffineFunction{Float64}, model::Model)
        val = func.constant
        for term in func.terms
            varref = JuMP.VariableRef(model, term.variable)  # 从 MOI index 得到 JuMP 变量
            variable_val = JuMP.start_value(varref)
            if isnan(variable_val)
                error("Variable $(varref) has no start_value set!")
            end
            val += term.coefficient * variable_val
        end
        return val
    end    
    
    bad_constraints = []
    
    for con in all_constraints(model; include_variable_in_set_constraints = true)
        con_name = name(con)
    
        func = MOI.get(model, MOI.ConstraintFunction(), con)
        set = MOI.get(model, MOI.ConstraintSet(), con)
    
        if !(func isa MOI.ScalarAffineFunction{Float64})
            println("Skipping non-affine constraint: $con_name")
            continue
        end
    
        lhs = eval_affine_function(func, model)
        violation = 0.0
    
        if isa(set, MOI.EqualTo)
            rhs = set.value
            violation = abs(lhs - rhs)
        elseif isa(set, MOI.LessThan)
            rhs = set.upper
            violation = max(lhs - rhs, 0)
        elseif isa(set, MOI.GreaterThan)
            rhs = set.lower
            violation = max(rhs - lhs, 0)
        else
            println("Unknown set type: $set")
            continue
        end
    
        if violation > 1e-4
            push!(bad_constraints, (con_name, violation))
        end
    end
    

    if isempty(bad_constraints)
        println("所有约束均满足！")
    else
        println("发现违反的约束 $(length(bad_constraints)) 个：")
        for (n, v) in bad_constraints
            println(" - $n violated by $v")
        end
    end

    return bad_constraints
end
