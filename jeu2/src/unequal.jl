#!/usr/bin/env julia
using JuMP
using CPLEX
using Random

# Update domains based on a cell's value
function updateDomains!(grid::Matrix{Int}, domains::Matrix{Set{Int}}, i::Int, j::Int, 
                       horizontal::Matrix{Int}, vertical::Matrix{Int})
    n = size(grid, 1)
    value = grid[i, j]
    
    if value == 0
        return
    end
    
    # Remove the value from domains in the same row
    for col in 1:n
        if col != j && value in domains[i, col]
            delete!(domains[i, col], value)
        end
    end
    
    # Remove the value from domains in the same column
    for row in 1:n
        if row != i && value in domains[row, j]
            delete!(domains[row, j], value)
        end
    end
    
    # Apply inequality constraints
    # Check right
    if j < n && horizontal[i, j] != 0
        if horizontal[i, j] == 1  # >
            for v in collect(domains[i, j+1])
                if v >= value
                    delete!(domains[i, j+1], v)
                end
            end
        else  # <
            for v in collect(domains[i, j+1])
                if v <= value
                    delete!(domains[i, j+1], v)
                end
            end
        end
    end
    
    # Check left
    if j > 1 && horizontal[i, j-1] != 0
        if horizontal[i, j-1] == 1  # Cell to the left > current cell
            for v in collect(domains[i, j-1])
                if v <= value
                    delete!(domains[i, j-1], v)
                end
            end
        else  # Cell to the left < current cell
            for v in collect(domains[i, j-1])
                if v >= value
                    delete!(domains[i, j-1], v)
                end
            end
        end
    end
    
    # Check below
    if i < n && vertical[i, j] != 0
        if vertical[i, j] == 1  # v
            for v in collect(domains[i+1, j])
                if v >= value
                    delete!(domains[i+1, j], v)
                end
            end
        else  # ^
            for v in collect(domains[i+1, j])
                if v <= value
                    delete!(domains[i+1, j], v)
                end
            end
        end
    end
    
    # Check above
    if i > 1 && vertical[i-1, j] != 0
        if vertical[i-1, j] == 1  # Cell above > current cell
            for v in collect(domains[i-1, j])
                if v <= value
                    delete!(domains[i-1, j], v)
                end
            end
        else  # Cell above < current cell
            for v in collect(domains[i-1, j])
                if v >= value
                    delete!(domains[i-1, j], v)
                end
            end
        end
    end
end

# Constraint propagation to reduce domains
function constraintPropagation!(grid::Matrix{Int}, domains::Matrix{Set{Int}}, 
                               horizontal::Matrix{Int}, vertical::Matrix{Int})
    n = size(grid, 1)
    changed = true
    
    while changed
        changed = false
        
        # Look for cells with only one possible value
        for i in 1:n
            for j in 1:n
                if grid[i, j] == 0 && length(domains[i, j]) == 1
                    grid[i, j] = first(domains[i, j])
                    updateDomains!(grid, domains, i, j, horizontal, vertical)
                    changed = true
                end
            end
        end
        
        # Look for values that can only appear in one cell of a row
        for i in 1:n
            for val in 1:n
                possible_cols = [j for j in 1:n if val in domains[i, j]]
                if length(possible_cols) == 1 && grid[i, possible_cols[1]] == 0
                    j = possible_cols[1]
                    grid[i, j] = val
                    domains[i, j] = Set([val])
                    updateDomains!(grid, domains, i, j, horizontal, vertical)
                    changed = true
                end
            end
        end
        
        # Look for values that can only appear in one cell of a column
        for j in 1:n
            for val in 1:n
                possible_rows = [i for i in 1:n if val in domains[i, j]]
                if length(possible_rows) == 1 && grid[possible_rows[1], j] == 0
                    i = possible_rows[1]
                    grid[i, j] = val
                    domains[i, j] = Set([val])
                    updateDomains!(grid, domains, i, j, horizontal, vertical)
                    changed = true
                end
            end
        end
    end
end

# Find cell with minimum domain size
function findMinDomainCell(grid::Matrix{Int}, domains::Matrix{Set{Int}})
    n = size(grid, 1)
    min_size = n + 1
    min_i, min_j = 0, 0
    
    for i in 1:n
        for j in 1:n
            if grid[i, j] == 0 && length(domains[i, j]) > 0 && length(domains[i, j]) < min_size
                min_size = length(domains[i, j])
                min_i, min_j = i, j
                if min_size == 2  # Early exit if we find a cell with only 2 options
                    return min_i, min_j
                end
            end
        end
    end
    
    return min_i, min_j
end

# Check if placing a value is consistent with the rules
function isConsistent(grid::Matrix{Int}, i::Int, j::Int, val::Int, 
                     horizontal::Matrix{Int}, vertical::Matrix{Int})
    n = size(grid, 1)
    
    # Check row
    for col in 1:n
        if grid[i, col] == val
            return false
        end
    end
    
    # Check column
    for row in 1:n
        if grid[row, j] == val
            return false
        end
    end
    
    # Check inequality constraints
    # Right neighbor
    if j < n && horizontal[i, j] != 0
        if horizontal[i, j] == 1 && grid[i, j+1] != 0 && val <= grid[i, j+1]
            return false
        elseif horizontal[i, j] == -1 && grid[i, j+1] != 0 && val >= grid[i, j+1]
            return false
        end
    end
    
    # Left neighbor
    if j > 1 && horizontal[i, j-1] != 0
        if horizontal[i, j-1] == 1 && grid[i, j-1] != 0 && grid[i, j-1] <= val
            return false
        elseif horizontal[i, j-1] == -1 && grid[i, j-1] != 0 && grid[i, j-1] >= val
            return false
        end
    end
    
    # Bottom neighbor
    if i < n && vertical[i, j] != 0
        if vertical[i, j] == 1 && grid[i+1, j] != 0 && val <= grid[i+1, j]
            return false
        elseif vertical[i, j] == -1 && grid[i+1, j] != 0 && val >= grid[i+1, j]
            return false
        end
    end
    
    # Top neighbor
    if i > 1 && vertical[i-1, j] != 0
        if vertical[i-1, j] == 1 && grid[i-1, j] != 0 && grid[i-1, j] <= val
            return false
        elseif vertical[i-1, j] == -1 && grid[i-1, j] != 0 && grid[i-1, j] >= val
            return false
        end
    end
    
    return true
end

# Backtracking algorithm
function backtrack!(grid::Matrix{Int}, domains::Matrix{Set{Int}}, 
                   horizontal::Matrix{Int}, vertical::Matrix{Int})
    n = size(grid, 1)
    
    # Check if the grid is complete
    if !any(grid .== 0)
        return true
    end
    
    # Find cell with minimum domain size (minimum remaining values heuristic)
    i, j = findMinDomainCell(grid, domains)
    
    if i == 0 || j == 0 || isempty(domains[i, j])
        return false  # No valid cell found or cell has empty domain
    end
    
    # Try each value in the domain
    for val in sort(collect(domains[i, j]))  # Sort for deterministic behavior
        # Check if the value is consistent with the constraints
        if isConsistent(grid, i, j, val, horizontal, vertical)
            # Make a copy of the current state
            grid_copy = copy(grid)
            domains_copy = deepcopy(domains)
            
            # Place the value
            grid[i, j] = val
            domains[i, j] = Set([val])
            
            # Update domains based on this assignment
            updateDomains!(grid, domains, i, j, horizontal, vertical)
            
            # Apply constraint propagation
            constraintPropagation!(grid, domains, horizontal, vertical)
            
            # Check if any domain became empty
            if !any(isempty.(domains[grid .== 0]))
                # Recurse
                if backtrack!(grid, domains, horizontal, vertical)
                    return true
                end
            end
            
            # If we reach here, this assignment failed, restore state
            grid .= grid_copy
            domains .= domains_copy
        end
    end
    
    return false  # No solution found
end

# Heuristic solver for unequal Sudoku using backtracking and constraint propagation
function heuristicSolve(instance_grid::Matrix{Int}, instance_horizontal::Matrix{Int}, instance_vertical::Matrix{Int})
    start_time = time()
    n = size(instance_grid, 1)
    
    # Create a deep copy of the grid to work with
    grid = copy(instance_grid)
    
    # Initialize domains for each cell (possible values)
    domains = [Set(1:n) for i in 1:n, j in 1:n]
    
    # Set fixed values and update domains
    for i in 1:n
        for j in 1:n
            if grid[i, j] > 0
                domains[i, j] = Set([grid[i, j]])
                # Update domains of other cells in the same row and column
                updateDomains!(grid, domains, i, j, instance_horizontal, instance_vertical)
            end
        end
    end
    
    # Apply constraint propagation to reduce domains
    constraintPropagation!(grid, domains, instance_horizontal, instance_vertical)
    
    # Use backtracking to solve the puzzle
    success = backtrack!(grid, domains, instance_horizontal, instance_vertical)
    
    solve_time = time() - start_time
    
    return success, solve_time, grid
end



# Solves the unequal Sudoku (Futoshiki) puzzle using CPLEX
function cplexSolve(instance_grid::Matrix{Int}, instance_horizontal::Matrix{Int}, instance_vertical::Matrix{Int})
    n = size(instance_grid, 1)
    
    # Start timing
    start_time = time()
    
    # Create optimization model
    model = Model(CPLEX.Optimizer)
    
    # Decision variables: x[i,j,k] = 1 if cell (i,j) contains number k
    @variable(model, x[1:n, 1:n, 1:n], Bin)
    
    # Each cell must contain exactly one value
    for i in 1:n
        for j in 1:n
            @constraint(model, sum(x[i,j,k] for k in 1:n) == 1)
        end
    end
    
    # Each row must contain each number exactly once
    for i in 1:n
        for k in 1:n
            @constraint(model, sum(x[i,j,k] for j in 1:n) == 1)
        end
    end
    
    # Each column must contain each number exactly once
    for j in 1:n
        for k in 1:n
            @constraint(model, sum(x[i,j,k] for i in 1:n) == 1)
        end
    end
    
    # Fixed cells from instance_grid
    for i in 1:n
        for j in 1:n
            if instance_grid[i,j] > 0
                @constraint(model, x[i,j,instance_grid[i,j]] == 1)
            end
        end
    end
    
    # Horizontal inequality constraints
    for i in 1:n
        for j in 1:n-1
            if instance_horizontal[i,j] == 1  # >
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) >= sum(k * x[i,j+1,k] for k in 1:n) + 1)
            elseif instance_horizontal[i,j] == -1  # <
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) + 1 <= sum(k * x[i,j+1,k] for k in 1:n))
            end
        end
    end
    
    # Vertical inequality constraints
    for i in 1:n-1
        for j in 1:n
            if instance_vertical[i,j] == 1  # v
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) >= sum(k * x[i+1,j,k] for k in 1:n) + 1)
            elseif instance_vertical[i,j] == -1  # ^
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) + 1 <= sum(k * x[i+1,j,k] for k in 1:n))
            end
        end
    end
    
    # Add a dummy objective (since we just need a feasible solution)
    @objective(model, Min, 0)
    
    # Set a time limit (e.g., 60 seconds)
    set_time_limit_sec(model, 60.0)
    
    # Suppress solver output
    set_silent(model)
    
    # Solve the model
    optimize!(model)
    
    # End timing
    solve_time = time() - start_time
    
    # Check if a solution was found
    if termination_status(model) == MOI.OPTIMAL
        # Extract solution
        solution_grid = zeros(Int, n, n)
        for i in 1:n
            for j in 1:n
                for k in 1:n
                    if value(x[i,j,k]) > 0.5
                        solution_grid[i,j] = k
                    end
                end
            end
        end
        return true, solve_time, solution_grid
    else
        return false, solve_time, zeros(Int, n, n)
    end
end



# Generates a random Latin square of size n×n
function generateLatinSquare(n::Int)
    # Start with a standard Latin square (Cayley table of the cyclic group)
    grid = zeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            grid[i, j] = mod((i + j - 2), n) + 1
        end
    end
    
    # Apply random row and column permutations to make it more random
    rows = randperm(n)
    cols = randperm(n)
    
    permuted_grid = zeros(Int, n, n)
    for i in 1:n
        for j in 1:n
            permuted_grid[i, j] = grid[rows[i], cols[j]]
        end
    end
    
    return permuted_grid
end

# Determines inequalities based on the Latin square
function determineInequalities(grid::Matrix{Int})
    n = size(grid, 1)
    
    # Initialize inequality matrices
    horizontal_ineq = zeros(Int, n, n-1)  # > is 1, < is -1, none is 0
    vertical_ineq = zeros(Int, n-1, n)    # v is 1, ^ is -1, none is 0
    
    # Fill horizontal inequalities
    for i in 1:n
        for j in 1:n-1
            if grid[i, j] > grid[i, j+1]
                horizontal_ineq[i, j] = 1  # >
            else
                horizontal_ineq[i, j] = -1  # <
            end
        end
    end
    
    # Fill vertical inequalities
    for i in 1:n-1
        for j in 1:n
            if grid[i, j] > grid[i+1, j]
                vertical_ineq[i, j] = 1  # v
            else
                vertical_ineq[i, j] = -1  # ^
            end
        end
    end
    
    return horizontal_ineq, vertical_ineq
end

# Creates instance data with a given density of clues
function createInstanceData(grid::Matrix{Int}, horizontal_ineq::Matrix{Int}, vertical_ineq::Matrix{Int}, density::Float64)
    n = size(grid, 1)
    
    # Total number of possible clues (grid cells + horizontal inequalities + vertical inequalities)
    total_clues = n*n + n*(n-1) + (n-1)*n
    
    # Number of clues to include based on density
    num_clues = Int(round(density * total_clues))
    
    # Create a list of all possible clue positions
    clue_positions = []
    
    # Grid cells
    for i in 1:n
        for j in 1:n
            push!(clue_positions, (:grid, i, j))
        end
    end
    
    # Horizontal inequalities
    for i in 1:n
        for j in 1:n-1
            push!(clue_positions, (:horizontal, i, j))
        end
    end
    
    # Vertical inequalities
    for i in 1:n-1
        for j in 1:n
            push!(clue_positions, (:vertical, i, j))
        end
    end
    
    # Shuffle and select a subset of positions
    shuffle!(clue_positions)
    selected_positions = clue_positions[1:num_clues]
    
    # Initialize empty instance data
    instance_grid = zeros(Int, n, n)
    instance_horizontal = zeros(Int, n, n-1)
    instance_vertical = zeros(Int, n-1, n)
    
    # Fill in the selected clues
    for (type, i, j) in selected_positions
        if type == :grid
            instance_grid[i, j] = grid[i, j]
        elseif type == :horizontal
            instance_horizontal[i, j] = horizontal_ineq[i, j]
        else # type == :vertical
            instance_vertical[i, j] = vertical_ineq[i, j]
        end
    end
    
    return instance_grid, instance_horizontal, instance_vertical
end

# Checks if a puzzle has a unique solution
function hasUniqueSolution(instance_grid::Matrix{Int}, instance_horizontal::Matrix{Int}, instance_vertical::Matrix{Int})
    # First, solve the puzzle normally
    success, time, solution = cplexSolve(instance_grid, instance_horizontal, instance_vertical)
    
    if !success
        return false, time, solution  # No solution found
    end
    
    # Create a new model to check for a different solution
    n = size(instance_grid, 1)
    model = Model(CPLEX.Optimizer)
    
    # Decision variables
    @variable(model, x[1:n, 1:n, 1:n], Bin)
    
    # Same constraints as in cplexSolve
    # Each cell must contain exactly one value
    for i in 1:n
        for j in 1:n
            @constraint(model, sum(x[i,j,k] for k in 1:n) == 1)
        end
    end
    
    # Each row must contain each number exactly once
    for i in 1:n
        for k in 1:n
            @constraint(model, sum(x[i,j,k] for j in 1:n) == 1)
        end
    end
    
    # Each column must contain each number exactly once
    for j in 1:n
        for k in 1:n
            @constraint(model, sum(x[i,j,k] for i in 1:n) == 1)
        end
    end
    
    # Fixed cells
    for i in 1:n
        for j in 1:n
            if instance_grid[i,j] > 0
                @constraint(model, x[i,j,instance_grid[i,j]] == 1)
            end
        end
    end
    
    # Horizontal inequality constraints
    for i in 1:n
        for j in 1:n-1
            if instance_horizontal[i,j] == 1  # >
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) >= sum(k * x[i,j+1,k] for k in 1:n) + 1)
            elseif instance_horizontal[i,j] == -1  # <
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) <= sum(k * x[i,j+1,k] for k in 1:n) - 1)
            end
        end
    end
    
    # Vertical inequality constraints
    for i in 1:n-1
        for j in 1:n
            if instance_vertical[i,j] == 1  # v
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) >= sum(k * x[i+1,j,k] for k in 1:n) + 1)
            elseif instance_vertical[i,j] == -1  # ^
                @constraint(model, sum(k * x[i,j,k] for k in 1:n) <= sum(k * x[i+1,j,k] for k in 1:n) - 1)
            end
        end
    end
    
    # Add constraint to ensure we find a different solution
    different_solution_constraint = @constraint(model, sum(x[i,j,solution[i,j]] for i in 1:n, j in 1:n if solution[i,j] > 0) <= sum(solution .> 0) - 1)
    
    # Set objective (dummy)
    @objective(model, Min, 0)
    
    # Suppress solver output
    set_silent(model)
    
    # Solve the model with a time limit
    set_time_limit_sec(model, 30.0)
    optimize!(model)
    
    # If no second solution found, the puzzle has a unique solution
    if termination_status(model) == MOI.INFEASIBLE
        return true, time, solution  # Unique solution
    else
        return false, time, solution  # Multiple solutions
    end
end

# Generates a valid unequal Sudoku (Futoshiki) instance
function generateInstance(n::Int, density::Float64; max_attempts::Int=120)
    for attempt in 1:max_attempts
        # Generate a random Latin square
        grid = generateLatinSquare(n)
        
        # Determine inequalities
        horizontal_ineq, vertical_ineq = determineInequalities(grid)
        
        # Create instance with given density
        instance_grid, instance_horizontal, instance_vertical = createInstanceData(grid, horizontal_ineq, vertical_ineq, density)
        
        # Check if the instance has a unique solution
        is_unique, solve_time, solution = hasUniqueSolution(instance_grid, instance_horizontal, instance_vertical)
        
        if is_unique
            return instance_grid, instance_horizontal, instance_vertical, grid, solve_time
        end
        
        println("Attempt $attempt: Instance doesn't have a unique solution. Trying again...")
    end
    
    error("Failed to generate a valid instance after $max_attempts attempts")
end

# Function to generate and save multiple instances
function generateDataSet(; output_dir::String="../data", n_values::Vector{Int}=[4, 5, 6], 
                          density_values::Vector{Float64}=[0.2, 0.3, 0.4], num_instances::Int=10)
    # Create output directory if it doesn't exist
    mkpath(output_dir)
    
    # Generate instances for each combination of parameters
    for n in n_values
        for density in density_values            
            for instance_num in 1:num_instances
                # Generate instance
                instance_grid, instance_horizontal, instance_vertical, solution_grid, solve_time = 
                    generateInstance(n, density)
                
                # Create a unique filename
                filename = "unequal_n$(n)_d$(Int(density*100))_$(instance_num).txt"
                filepath = joinpath(output_dir, filename)
                
                # Save the instance to a file
                open(filepath, "w") do file
                    # Write instance metadata
                    write(file, "# Unequal Sudoku instance\n")
                    write(file, "# n: $n\n")
                    write(file, "# density: $density\n")
                    write(file, "# solve_time: $solve_time\n\n")
                    
                    # Write instance grid
                    write(file, "GRID\n")
                    for i in 1:n
                        line = join(instance_grid[i,:], ",")
                        write(file, line * "\n")
                    end
                    
                    # Write horizontal inequalities
                    write(file, "\nHORIZONTAL\n")
                    for i in 1:n
                        line = join(instance_horizontal[i,:], ",")
                        write(file, line * "\n")
                    end
                    
                    # Write vertical inequalities
                    write(file, "\nVERTICAL\n")
                    for i in 1:n-1
                        line = join(instance_vertical[i,:], ",")
                        write(file, line * "\n")
                    end
                    
                    # Write solution grid
                    write(file, "\nSOLUTION\n")
                    for i in 1:n
                        line = join(solution_grid[i,:], ",")
                        write(file, line * "\n")
                    end
                end
            end
        end
    end
    println("Dataset generation complete.")
end

# Read an instance from a file
function readInputFile(path::String)
    lines = readlines(path)
    
    # Extract metadata
    n = nothing
    mode = "NONE"
    grid_data = []
    horizontal_data = []
    vertical_data = []
    solution_data = []
    
    for line in lines
        # Skip empty lines and comments
        if isempty(strip(line)) || startswith(strip(line), "#")
            continue
        end
        
        # Check for section headers
        if line == "GRID"
            mode = "GRID"
            continue
        elseif line == "HORIZONTAL"
            mode = "HORIZONTAL"
            continue
        elseif line == "VERTICAL"
            mode = "VERTICAL"
            continue
        elseif line == "SOLUTION"
            mode = "SOLUTION"
            continue
        end
        
        # Extract n from the first data line if not already set
        if mode == "GRID" && isempty(grid_data)
            n = length(split(line, ","))
        end
        
        # Parse data based on current mode
        if mode == "GRID"
            push!(grid_data, parse.(Int, split(line, ",")))
        elseif mode == "HORIZONTAL"
            push!(horizontal_data, parse.(Int, split(line, ",")))
        elseif mode == "VERTICAL"
            push!(vertical_data, parse.(Int, split(line, ",")))
        elseif mode == "SOLUTION"
            push!(solution_data, parse.(Int, split(line, ",")))
        end
    end
    
    # Convert to matrices
    grid = Matrix{Int}(hcat(grid_data...)')
    horizontal = Matrix{Int}(hcat(horizontal_data...)')
    vertical = Matrix{Int}(hcat(vertical_data...)')
    solution = Matrix{Int}(hcat(solution_data...)')
    
    return grid, horizontal, vertical, solution
end

# Function to display a grid with inequalities
function displayInstance(instance_grid::Matrix{Int}, instance_horizontal::Matrix{Int}, instance_vertical::Matrix{Int})
    n = size(instance_grid, 1)
    
    for i in 1:n
        # Print row of numbers
        for j in 1:n
            print(instance_grid[i,j] == 0 ? "." : instance_grid[i,j], " ")
            
            # Print horizontal inequality if not the last column
            if j < n
                if instance_horizontal[i,j] == 1
                    print("> ")
                elseif instance_horizontal[i,j] == -1
                    print("< ")
                else
                    print("  ")
                end
            end
        end
        println()
        
        # Print vertical inequalities if not the last row
        if i < n
            for j in 1:n
                if instance_vertical[i,j] == 1
                    print("v ")
                elseif instance_vertical[i,j] == -1
                    print("^ ")
                else
                    print("  ")
                end
                
                # Add spacing for alignment
                if j < n
                    print("  ")
                end
            end
            println()
        end
    end
end

# Function to display a solution
function displaySolution(solution::Matrix{Int})
    n = size(solution, 1)
    
    for i in 1:n
        for j in 1:n
            print(solution[i,j], " ")
        end
        println()
    end
end

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()

    dataFolder = "../data/"
    resFolder = "../res/"

    # Array which contains the name of the resolution methods
    resolutionMethod = ["cplex", "heuristic"]
    #resolutionMethod = ["cplex", "heuristique"]

    # Array which contains the result folder of each resolution method
    resolutionFolder = resFolder .* resolutionMethod

    # Create each result folder if it does not exist
    for folder in resolutionFolder
        if !isdir(folder)
            mkdir(folder)
        end
    end
            
    global isOptimal = false
    global solveTime = -1

    # For each instance
    for file in filter(x->occursin(".txt", x), readdir(dataFolder))  
        
        println("-- Resolution of ", file)
        grid, horizontal, vertical, expected_solution = readInputFile(dataFolder * file)
        
        # For each resolution method
        for methodId in 1:size(resolutionMethod, 1)
            
            outputFile = resolutionFolder[methodId] * "/" * file

            # If the instance has not already been solved by this method
            if !isfile(outputFile)
                
                fout = open(outputFile, "w")  

                resolutionTime = -1
                isOptimal = false
                
                # If the method is cplex
                if resolutionMethod[methodId] == "cplex"
                    
                    # Solve it and get the results
                    isOptimal, resolutionTime, solution = cplexSolve(grid, horizontal, vertical)
                    
                    # If a solution is found, write it
                    if isOptimal
                        println(fout, "solution = [")
                        for i in 1:size(solution, 1)
                            println(fout, "    ", join(solution[i,:], ", "), ",")
                        end
                        println(fout, "]")
                    end

                # If the method is one of the heuristics
                else
                    
                    isSolved = false

                    # Start a chronometer 
                    startingTime = time()
                    
                    # While the grid is not solved and less than 100 seconds are elapsed
                    while !isOptimal && resolutionTime < 100
                        
                        # Solve it and get the results
                        isOptimal, resolutionTime, solution = heuristicSolve(grid, horizontal, vertical)

                        # Stop the chronometer
                        resolutionTime = time() - startingTime
                        
                    end

                    # Write the solution (if any)
                    if isOptimal
                        println(fout, "solution = [")
                        for i in 1:size(solution, 1)
                            println(fout, "    ", join(solution[i,:], ", "), ",")
                        end
                        println(fout, "]")
                    end 
                end

                println(fout, "solveTime = ", resolutionTime) 
                println(fout, "isOptimal = ", isOptimal)
                
                # Write problem information if needed
                println(fout, "original_grid = ", grid)
                println(fout, "horizontal = ", horizontal)
                println(fout, "vertical = ", vertical)
                
                close(fout)
            end

            # Display the results obtained with the method on the current instance
            include(outputFile)
            println(resolutionMethod[methodId], " optimal: ", isOptimal)
            println(resolutionMethod[methodId], " time: " * string(round(solveTime, sigdigits=2)) * "s\n")
        end         
    end 
end

"""
Create a pdf file which contains a performance diagram associated to the results of the ../res folder
Display one curve for each subfolder of the ../res folder.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "solveTime" and a variable "isOptimal"
"""
function performanceDiagram(outputFile::String)

    resultFolder = "../res/"
    
    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    folderName = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file
        
        # If it is a subfolder
        if isdir(path)
            
            folderName = vcat(folderName, file)
             
            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # Array that will contain the resolution times (one line for each subfolder)
    results = Array{Float64}(undef, subfolderCount, maxSize)

    for i in 1:subfolderCount
        for j in 1:maxSize
            results[i, j] = Inf
        end
    end

    folderCount = 0
    maxSolveTime = 0

    # For each subfolder
    for file in readdir(resultFolder)
            
        path = resultFolder * file
        
        if isdir(path)

            folderCount += 1
            fileCount = 0

            # For each text file in the subfolder
            for resultFile in filter(x->occursin(".txt", x), readdir(path))

                fileCount += 1
                include(path * "/" * resultFile)

                if isOptimal
                    results[folderCount, fileCount] = solveTime

                    if solveTime > maxSolveTime
                        maxSolveTime = solveTime
                    end 
                end 
            end 
        end
    end 

    # Sort each row increasingly
    results = sort(results, dims=2)

    println("Max solve time: ", maxSolveTime)

    # For each line to plot
    for dim in 1: size(results, 1)

        x = Array{Float64, 1}()
        y = Array{Float64, 1}()

        # x coordinate of the previous inflexion point
        previousX = 0
        previousY = 0

        append!(x, previousX)
        append!(y, previousY)
            
        # Current position in the line
        currentId = 1

        # While the end of the line is not reached 
        while currentId != size(results, 2) && results[dim, currentId] != Inf

            # Number of elements which have the value previousX
            identicalValues = 1

             # While the value is the same
            while results[dim, currentId] == previousX && currentId <= size(results, 2)
                currentId += 1
                identicalValues += 1
            end

            # Add the proper points
            append!(x, previousX)
            append!(y, currentId - 1)

            if results[dim, currentId] != Inf
                append!(x, results[dim, currentId])
                append!(y, currentId - 1)
            end
            
            previousX = results[dim, currentId]
            previousY = currentId - 1
            
        end

        append!(x, maxSolveTime)
        append!(y, currentId - 1)

        # If it is the first subfolder
        if dim == 1

            # Draw a new plot
            plot(x, y, label = folderName[dim], legend = :bottomright, xaxis = "Time (s)", yaxis = "Solved instances",linewidth=3)

        # Otherwise 
        else
            # Add the new curve to the created plot
            savefig(plot!(x, y, label = folderName[dim], linewidth=3), outputFile)
        end 
    end
end 

"""
Create a latex file which contains an array with the results of the ../res folder.
Each subfolder of the ../res folder contains the results of a resolution method.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "solveTime" and a variable "isOptimal"
"""
function resultsArray(outputFile::String)
    
    resultFolder = "../res/"
    dataFolder = "../data/"
    
    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    # Open the latex output file
    fout = open(outputFile, "w")

    # Print the latex file output
    println(fout, raw"""\documentclass{article}

\usepackage[french]{babel}
\usepackage [utf8] {inputenc} % utf-8 / latin1 
\usepackage{multicol}

\setlength{\hoffset}{-18pt}
\setlength{\oddsidemargin}{0pt} % Marge gauche sur pages impaires
\setlength{\evensidemargin}{9pt} % Marge gauche sur pages paires
\setlength{\marginparwidth}{54pt} % Largeur de note dans la marge
\setlength{\textwidth}{481pt} % Largeur de la zone de texte (17cm)
\setlength{\voffset}{-18pt} % Bon pour DOS
\setlength{\marginparsep}{7pt} % Séparation de la marge
\setlength{\topmargin}{0pt} % Pas de marge en haut
\setlength{\headheight}{13pt} % Haut de page
\setlength{\headsep}{10pt} % Entre le haut de page et le texte
\setlength{\footskip}{27pt} % Bas de page + séparation
\setlength{\textheight}{668pt} % Hauteur de la zone de texte (25cm)

\begin{document}""")

    header = raw"""
\begin{center}
\renewcommand{\arraystretch}{1.4} 
 \begin{tabular}{l"""

    # Name of the subfolder of the result folder (i.e, the resolution methods used)
    folderName = Array{String, 1}()

    # List of all the instances solved by at least one resolution method
    solvedInstances = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file
        
        # If it is a subfolder
        if isdir(path)

            # Add its name to the folder list
            folderName = vcat(folderName, file)
             
            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            # Add all its files in the solvedInstances array
            for file2 in filter(x->occursin(".txt", x), readdir(path))
                solvedInstances = vcat(solvedInstances, file2)
            end 

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # 修复：定义一个函数来提取文件名中的数字部分
    function extractNumberFromFilename(filename::String)
        # 匹配文件名末尾的数字部分，即"_数字.txt"
        m = match(r"_(\d+)\.txt$", filename)
        if m !== nothing
            return parse(Int, m.captures[1])
        else
            # 如果没有找到符合模式的数字，则保持原样
            return filename
        end
    end
    
    # 基于文件名中相同的部分和序号进行排序
    function sortInstances(filename::String)
        # 分离文件名的基本部分和数字部分
        # 例如，从"unequal_n5_d30_7.txt"提取"unequal_n5_d30_"和"7"
        base_end = findlast('_', filename) # 最后一个下划线的位置
        if base_end !== nothing
            base_name = filename[1:base_end]
            num_part = extractNumberFromFilename(filename)
            # 返回排序元组：先按基本名称排序，然后按数字排序
            return (base_name, num_part)
        else
            # 如果文件名不符合预期格式，则按原始文件名排序
            return (filename, 0)
        end
    end

    # 去重
    unique!(solvedInstances)
    
    # 按照新定义的规则排序
    sort!(solvedInstances, by=sortInstances)

    # For each resolution method, add two columns in the array
    for folder in folderName
        header *= "rr"
    end

    header *= "}\n\t\\hline\n"

    # Create the header line which contains the methods name
    for folder in folderName
        header *= " & \\multicolumn{2}{c}{\\textbf{" * folder * "}}"
    end

    header *= "\\\\\n\\textbf{Instance} "

    # Create the second header line with the content of the result columns
    for folder in folderName
        header *= " & \\textbf{Temps (s)} & \\textbf{Optimal ?} "
    end

    header *= "\\\\\\hline\n"

    footer = raw"""\hline\end{tabular}
\end{center}

"""
    println(fout, header)

    # On each page an array will contain at most maxInstancePerPage lines with results
    maxInstancePerPage = 27
    id = 1

    # For each solved files
    for solvedInstance in solvedInstances

        # If we do not start a new array on a new page
        if rem(id, maxInstancePerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end 

        # Replace the potential underscores '_' in file names
        print(fout, replace(solvedInstance, "_" => "\\_"))

        # For each resolution method
        for method in folderName

            path = resultFolder * method * "/" * solvedInstance

            # If the instance has been solved by this method
            if isfile(path)

                include(path)

                println(fout, " & ", round(solveTime, digits=2), " & ")

                if isOptimal
                    println(fout, "\$\\times\$")
                end 
                
            # If the instance has not been solved by this method
            else
                println(fout, " & - & - ")
            end
        end

        println(fout, "\\\\")

        id += 1
    end

    # Print the end of the latex file
    println(fout, footer)

    println(fout, "\\end{document}")

    close(fout)
    
end 

generateDataSet()
solveDataSet()
performanceDiagram("../res/performance_diagram.pdf")
resultsArray("../res/results.tex")