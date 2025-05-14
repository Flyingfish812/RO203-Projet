using CPLEX

function build_model(game::PalisadeGame)
    h, w = game.height, game.width
    k = game.region_size
    V = [(i, j) for i in 1:h, j in 1:w]  # 所有格子
    E = []

    # Edges (only consider right and down to avoid duplicates)
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

    R = collect(1:(h*w)÷k)  # Region numbers

    m = Model(CPLEX.Optimizer)

    # ========= Variables =========

    @variable(m, wvars[e in E], Bin, base_name = "w")  # Edges (side e = ((i,j),(i',j')))
    @variable(m, x[i=1:h, j=1:w, r in R], Bin, base_name = "x")  # Region membership variables
    @variable(m, y[i=1:h, j=1:w, r in R], Bin, base_name = "y")  # Region root variables
    @variable(m, f[i=1:h, j=1:w, i2=1:h, j2=1:w, r in R] >= 0, base_name = "f")  # Connectivity flow variables

    # ========= Constraints =========

    # Each cell can only belong to one region
    @constraint(m, [i=1:h, j=1:w], sum(x[i,j,r] for r in R) == 1, base_name="one_region")

    # Each region contains exactly k cells
    @constraint(m, [r in R], sum(x[i,j,r] for (i,j) in V) == k, base_name="region_size")

    # Each region has one root
    @constraint(m, [r in R], sum(y[i,j,r] for (i,j) in V) == 1, base_name="region_root_unique")
    @constraint(m, [i=1:h, j=1:w, r in R], y[i,j,r] <= x[i,j,r], base_name="root_in_region")

    # If neighboring cells are not in the same region, there is a wall
    for ((i1,j1), (i2,j2)) in E, r in R
        key = ((i1,j1),(i2,j2))
        @constraint(m, x[i1,j1,r] - x[i2,j2,r] <= wvars[key], base_name="wall_diff_pos")
        @constraint(m, x[i2,j2,r] - x[i1,j1,r] <= wvars[key], base_name="wall_diff_neg")
        @constraint(m, x[i1,j1,r] + x[i2,j2,r] + wvars[key] <= 2, base_name="wall_bound")
    end
    
    # Clue constraints
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

    # Connectivity constraints: flow balance (root with outflow k-1, the others with a flow -1)
    for (i,j) in V, r in R
        neighbors = [(i+di, j+dj) for (di,dj) in ((0,1), (1,0), (0,-1), (-1,0)) if 1 <= i+di <= h && 1 <= j+dj <= w]
        inflow = sum(f[i2,j2,i,j,r] for (i2,j2) in neighbors)
        outflow = sum(f[i,j,i2,j2,r] for (i2,j2) in neighbors)
        @constraint(m, outflow - inflow == k*y[i,j,r] - x[i,j,r], base_name="flow_balance_$i$j$r")
    end

    # Flow can only go through walls
    for ((i1,j1),(i2,j2)) in E, r in R
        key = ((i1,j1),(i2,j2))
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*(1 - wvars[key]), base_name="flow_wall")
        @constraint(m, f[i2,j2,i1,j1,r] <= (k-1)*(1 - wvars[key]), base_name="flow_wall")
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*x[i1,j1,r], base_name="flow_xfrom")
        @constraint(m, f[i1,j1,i2,j2,r] <= (k-1)*x[i2,j2,r], base_name="flow_xto")
    end

    # ========= Objective =========
    @objective(m, Min, 0)

    return m
end

function extractSolution(game::PalisadeGame, ds::DisjointSet, model::Model)
    h, w = game.height, game.width
    k = game.region_size
    R = collect(1:(h*w)÷k)

    # 重置墙体
    game.h_walls .= Int8(-1)
    game.v_walls .= Int8(-1)

    # 获取墙体变量名称集合
    for i in 1:h
        for j in 1:w
            if j < w
                val = round(Int, value(model[:wvars][((i,j),(i,j+1))]))
                game.v_walls[i, j+1] = val
            end
            if i < h
                val = round(Int, value(model[:wvars][((i,j),(i+1,j))]))
                game.h_walls[i+1, j] = val
            end
        end
    end

    # 边界墙体保持为 1（必要）
    game.h_walls[1, :] .= 1
    game.h_walls[end, :] .= 1
    game.v_walls[:, 1] .= 1
    game.v_walls[:, end] .= 1

    # 重建并查集
    ds.parent = collect(1:h*w)
    ds.size = ones(Int, h*w)

    # 还原每个格子属于哪个区域
    grid_region = fill(0, h, w)
    for r in R
        for i in 1:h, j in 1:w
            if round(Int, value(model[:x][i,j,r])) == 1
                grid_region[i,j] = r
            end
        end
    end

    # 以区域编号分组合并格子（更新并查集）
    for r in R
        cells = [(i,j) for i in 1:h, j in 1:w if grid_region[i,j] == r]
        if length(cells) <= 1
            continue
        end
        base = cells[1]
        for c in cells[2:end]
            connect(game, ds, base, c)
        end
    end

    return game, ds
end

function export_solution(model::Model, filepath::String)
    println("正在导出模型解至文件：$filepath")

    solution = Dict{String, Dict{String, Float64}}()

    for var in all_variables(model)
        val = try
            value(var)
        catch
            continue
        end

        if isnothing(val)
            continue
        end

        val = round(val, digits=6)
        var_name = name(var)

        # 拆解变量名和索引（形如 x[1,2,3] → name = x, key = (1,2,3)）
        if occursin("[", var_name)
            varname = first(split(var_name, "[")) * "vars"
            index_str = "[" * join(split(var_name, "[")[2:end], "[")  # 保留括号部分
            index_str = strip(index_str, ['[', ']'])  # 变成 "1,2,3"
            index_str = "(" * index_str * ")"         # 变成元组字符串 "(1,2,3)"
        else
            varname = var_name
            index_str = "()"
        end

        if !haskey(solution, varname)
            solution[varname] = Dict{String, Float64}()
        end
        solution[varname][index_str] = val
    end

    open(filepath, "w") do io
        JSON3.write(io, solution; indent=2)
    end

    println("导出完成，变量组：", join(keys(solution), ", "))
end

export_path = "solution.json"