using DataStructures  # For disjoint set

# Directions: up, right, down, left
const DIRECTIONS = [(-1, 0), (0, 1), (1, 0), (0, -1)]

# Boundary check for grid coordinates
function in_bounds(game::PalisadeGame, i, j)
    return 1 <= i <= game.height && 1 <= j <= game.width
end

# ------------ Helper functions ----------------------------

# Convert (i,j) coordinates to a 1D index for the disjoint set
function flatten_coord(game::PalisadeGame, i::Int, j::Int)
    return (i - 1) * game.width + j
end

# Convert a 1D index back to (i,j) coordinates
function unflatten_coord(game::PalisadeGame, index::Int)
    i = div(index - 1, game.width) + 1
    j = mod(index - 1, game.width) + 1
    return (i, j)
end

# Get the region index of a grid
function region(game::PalisadeGame, ds::DisjointSet, i::Int, j::Int)
    index = flatten_coord(game, i, j)
    while ds.parent[index] != index
        index = ds.parent[index]
    end
    return index
end

# Whether two grids are in the same region
function connected(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    return region(game, ds, a[1], a[2]) == region(game, ds, b[1], b[2])
end

# Connect the two regions of the grids, set the parents as the minimum index, and the wall between them to 0
function connect(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    ra = region(game, ds, a[1], a[2])
    rb = region(game, ds, b[1], b[2])

    if ra == rb
        return 0  # No need to connect
    end

    # Boarder walls are set to 0
    count = 0
    for (dir, r, c) in boarder(game, ds, a, b)
        if get_wall(game, dir, r, c) == -1
            set_wall!(game, dir, r, c, Int8(0))
            count += 1
        end
    end

    # Get all the indices of the two regions
    indices_a = [i for i in eachindex(ds.parent) if region(game, ds, unflatten_coord(game, i)...) == ra]
    indices_b = [i for i in eachindex(ds.parent) if region(game, ds, unflatten_coord(game, i)...) == rb]
    combined = vcat(indices_a, indices_b)
    new_rep = minimum(combined)
    new_size = length(combined)

    for i in combined
        ds.parent[i] = new_rep
        ds.size[i] = new_size
    end

    return count  # return the number of walls set
end


# Disconnect two regions, set the walls between them to 1
function disconnect(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    walls = boarder(game, ds, a, b)
    count = 0
    for (dir, r, c) in walls
        if get_wall(game, dir, r, c) == -1
            set_wall!(game, dir, r, c, 1)
            count += 1
        end
    end
    return count
end

# Check if two regions are disconnected
function disconnected(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    return region(game, ds, a[1], a[2]) != region(game, ds, b[1], b[2])
end

# Get all cells in each region
function get_region_cells(game::PalisadeGame, ds::DisjointSet)
    regions = Dict{Int, Vector{Tuple{Int,Int}}}()
    for i in 1:game.height, j in 1:game.width
        rep = region(game, ds, i, j)
        if !haskey(regions, rep)
            regions[rep] = []
        end
        push!(regions[rep], (i, j))
    end
    return regions
end

# Get the status of the walls around a cell, return by [up, right, down, left], 1=wall, 0=no wall, -1=unknown
function wall(game::PalisadeGame, i::Int, j::Int)
    up    = game.h_walls[i, j]
    right = game.v_walls[i, j+1]
    down  = game.h_walls[i+1, j]
    left  = game.v_walls[i, j]
    return (up, right, down, left)
end

# Get the wall between two cells, return (:h or :v, row, col) indicating the wall position
function get_wall_between(game::PalisadeGame, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    di, dj = b[1] - a[1], b[2] - a[2]
    if (di, dj) == (-1, 0)
        return (:h, a[1], a[2])
    elseif (di, dj) == (1, 0)
        return (:h, b[1], b[2])
    elseif (di, dj) == (0, -1)
        return (:v, a[1], a[2])
    elseif (di, dj) == (0, 1)
        return (:v, b[1], b[2])
    else
        error("Not neighboring, cannot get the wall")
    end
end

# Set the wall value (1=wall, 0=no wall)
function set_wall!(game::PalisadeGame, dir::Symbol, r::Int, c::Int, value::Int8)
    if dir == :h
        game.h_walls[r, c] = value
    elseif dir == :v
        game.v_walls[r, c] = value
    else
        error("Unknown wall direction: $dir")
    end
end

# Get the current status of the wall
function get_wall(game::PalisadeGame, dir::Symbol, r::Int, c::Int)
    return dir == :h ? game.h_walls[r, c] : game.v_walls[r, c]
end

# Get the neighboring cells of the current cell (up, right, down, left)
function neighbors(game::PalisadeGame, pos::Tuple{Int,Int})
    i, j = pos
    res = Tuple{Int,Int}[]
    for (di, dj) in DIRECTIONS
        ni, nj = i + di, j + dj
        if in_bounds(game, ni, nj)
            push!(res, (ni, nj))
        end
    end
    return res
end

# Get the representative cells of all "neighboring regions" of the current cell (only one neighbor cell for each different region)
function neighbor(game::PalisadeGame, ds::DisjointSet, pos::Tuple{Int,Int})
    root = region(game, ds, pos[1], pos[2])
    seen_regions = Set{Int}()
    result = Set{Tuple{Int,Int}}()

    # For all the cells in the region of the current cell
    region_cells = get_region_cells(game, ds)[root]
    for cell in region_cells
        for ncell in neighbors(game, cell)
            nroot = region(game, ds, ncell[1], ncell[2])
            if nroot != root && !(nroot in seen_regions)
                push!(seen_regions, nroot)
                push!(result, ncell)  # Push the representative
            end
        end
    end
    return collect(result)
end

# Check if two regions are adjacent and return all walls between them (as a list)
function boarder(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    if region(game, ds, a[1], a[2]) == region(game, ds, b[1], b[2])
        return []
    end
    walls = []
    for cell_a in get_region_cells(game, ds)[region(game, ds, a[1], a[2])]
        for cell_b in neighbors(game, cell_a)
            if region(game, ds, cell_b[1], cell_b[2]) == region(game, ds, b[1], b[2])
                push!(walls, get_wall_between(game, cell_a, cell_b))
            end
        end
    end
    return walls
end

# Check if two cells belong to different regions and are adjacent
# if the wall between them is 1 (already a wall), or if connecting them exceeds the size limit, return false
function canconnect(game::PalisadeGame, ds::DisjointSet, a::Tuple{Int,Int}, b::Tuple{Int,Int})
    if !in_bounds(game, a...) || !in_bounds(game, b...)
        return false
    end
    if connected(game, ds, a, b)
        return false
    end
    walls = boarder(game, ds, a, b)
    for (dir, r, c) in walls
        if get_wall(game, dir, r, c) == 1
            return false
        end
    end
    # Check the size limit
    ra = region(game, ds, a[1], a[2])
    rb = region(game, ds, b[1], b[2])
    total_size = ds.size[ra] + ds.size[rb]
    return total_size <= game.region_size
end

# Get all "connectable neighboring cells" of the current cell in the same region
function extendlist(game::PalisadeGame, ds::DisjointSet, pos::Tuple{Int,Int})
    root = region(game, ds, pos[1], pos[2])
    region_cells = get_region_cells(game, ds)[root]

    seen_roots = Set{Int}()
    result = Set{Tuple{Int,Int}}()

    for cell in region_cells
        for ncell in neighbors(game, cell)
            nroot = region(game, ds, ncell[1], ncell[2])
            if nroot != root && !(nroot in seen_roots)
                if canconnect(game, ds, cell, ncell)
                    push!(seen_roots, nroot)
                    push!(result, ncell)
                end
            end
        end
    end

    return collect(result)
end

# -------------- Objective functions and solvers ----------------------------

# Step 1: Set walls between connected clues to 1 if they exceed the region size
function connected_clues_vs_region(game::PalisadeGame, ds::DisjointSet)
    changes = 0
    seen = Set{Tuple{Tuple{Int,Int}, Tuple{Int,Int}}}()

    for ((i, j), p) in game.clues
        for (ni, nj) in neighbors(game, (i, j))
            if !haskey(game.clues, (ni, nj))
                continue
            end

            q = game.clues[(ni, nj)]
            a, b = (i, j), (ni, nj)
            pair = a < b ? (a, b) : (b, a)  # Avoid duplicate pairs

            if pair in seen
                continue
            end
            push!(seen, pair)

            if (8 - p - q > game.region_size) || (p == 3 && q == 3 && game.region_size > 2)
                # Set the wall between them to 1
                dir, r, c = get_wall_between(game, a, b)
                if get_wall(game, dir, r, c) == -1
                    set_wall!(game, dir, r, c, Int8(1))
                    changes += 1
                end
            end
        end
    end

    return changes  # Number of walls set
end

# Step 2: Set walls to 0 if the number of walls is exhausted
function number_exhausted(game::PalisadeGame, ds::DisjointSet)
    changes = 0

    for ((i, j), clue) in game.clues
        wall_vals = wall(game, i, j)
        count_1 = count(x -> x == 1, wall_vals)
        count_0 = count(x -> x == 0, wall_vals)
        count_neg1 = count(x -> x == -1, wall_vals)

        if count_neg1 == 0
            continue
        end

        if clue - count_1 == 0
            # All the unknown walls must be 0 if clue exhausted
            for k in 1:4
                if wall_vals[k] == -1
                    di, dj = DIRECTIONS[k]
                    ni, nj = i + di, j + dj
                    if in_bounds(game, ni, nj)
                        dir, r, c = get_wall_between(game, (i,j), (ni,nj))
                        if get_wall(game, dir, r, c) == -1
                            connect(game, ds, (i,j), (ni,nj))
                            changes += 1
                        end
                    end
                end
            end
        elseif clue - count_1 == count_neg1
            # All the unknown walls must be 1 if walls "fit" the clue
            for k in 1:4
                if wall_vals[k] == -1
                    di, dj = DIRECTIONS[k]
                    ni, nj = i + di, j + dj
                    if in_bounds(game, ni, nj)
                        dir, r, c = get_wall_between(game, (i,j), (ni,nj))
                        if get_wall(game, dir, r, c) == -1
                            set_wall!(game, dir, r, c, Int8(1))
                            changes += 1
                        end
                    end
                end
            end
        end
    end
    return changes  # Number of walls set
end

# Step 3: Set walls to 1 if the size of the region exceeds the limit
function not_too_big(game::PalisadeGame, ds::DisjointSet)
    regions = get_region_cells(game, ds)
    changes = 0
    seen = Set{Tuple{Int,Int}}()

    for (rep1, cells1) in regions
        size1 = ds.size[rep1]

        # Get all the neighboring regions of the current region
        for cell1 in cells1
            for cell2 in neighbors(game, cell1)
                rep2 = region(game, ds, cell2...)
                if rep1 == rep2
                    continue
                end
                key = rep1 < rep2 ? (rep1, rep2) : (rep2, rep1)
                if key in seen
                    continue
                end
                push!(seen, key)

                size2 = ds.size[rep2]
                if size1 + size2 > game.region_size
                    # Set all the unknown walls between the two regions to 1
                    for (dir, r, c) in boarder(game, ds, cell1, cell2)
                        if get_wall(game, dir, r, c) == -1
                            set_wall!(game, dir, r, c, Int8(1))
                            changes += 1
                        end
                    end
                end
            end
        end
    end
    return changes  # Number of walls set
end

# Step 4: Merge small regions if they are adjacent to only one other region
function not_too_small(game::PalisadeGame, ds::DisjointSet)
    regions = get_region_cells(game, ds)
    changes = 0
    visited = Set{Int}()

    for (rep, cells) in regions
        if rep in visited
            continue
        end

        if ds.size[rep] >= game.region_size
            continue
        end

        extend_targets = Set{Int}()
        candidate = nothing

        # Find all the "extendable" neighboring regions of the current region
        for cell in cells
            for npos in extendlist(game, ds, cell)
                nroot = region(game, ds, npos...)
                if nroot != rep
                    push!(extend_targets, nroot)
                    candidate = (cell, npos)  # Choose one of the neighboring cells to connect
                end
            end
        end

        if length(extend_targets) == 1 && !isnothing(candidate)
            a, b = candidate
            changes += connect(game, ds, a, b)
            # Two regions are merged, so we mark them as visited
            push!(visited, region(game, ds, a...))
            push!(visited, region(game, ds, b...))
        end
    end

    return changes
end

# Step 5: No dangling edges (walls that have an end)
function no_dangling_edges(game::PalisadeGame, ds::DisjointSet)
    changes = 0
    h, w = game.height, game.width

    for i in 1:h
        for j in 1:w
            walls = [
                (:h, i, j, game.h_walls[i, j]),
                (:h, i+1, j, game.h_walls[i+1, j]),
                (:v, i, j, game.v_walls[i, j]),
                (:v, i, j+1, game.v_walls[i, j+1])
            ]

            values = getindex.(walls, 4)
            count_neg1 = count(==( -1), values)
            count_0 = count(==(  0), values)
            count_1 = count(==(  1), values)

            if count_neg1 == 1 && count_0 == 2 && count_1 == 1
                for (dir, r, c, val) in walls
                    if val == -1
                        set_wall!(game, dir, r, c, Int8(1))
                        changes += 1
                        break
                    end
                end
            end
        end
    end
    return changes
end

# Step 6: Equivalent edges (walls that can be set to 1 or 0)
function equivalent_edges(game::PalisadeGame, ds::DisjointSet)
    changes = 0

    for ((i, j), clue) in game.clues
        neigh = neighbors(game, (i, j))
        root_self = region(game, ds, i, j)

        # Calculate the number of walls around the current cell
        wall_vals = wall(game, i, j)
        wall_count = count(x -> x == 1, wall_vals)

        # Generate all pairs of neighboring cells
        for m in 1:length(neigh)-1
            for n in m+1:length(neigh)
                a, b = neigh[m], neigh[n]
                root_a = region(game, ds, a...)
                root_b = region(game, ds, b...)

                # If two neighboring cells belong to the same region and that region is not the current cell's region
                if root_a == root_b && root_a != root_self
                    combined_size = ds.size[root_self] + ds.size[root_a]

                    # Case A: The size of the two regions exceeds the limit, set the wall to 1
                    if combined_size > game.region_size
                        for x in (a, b)
                            dir, r, c = get_wall_between(game, (i,j), x)
                            if get_wall(game, dir, r, c) == -1
                                set_wall!(game, dir, r, c, Int8(1))
                                changes += 1
                            end
                        end

                    # Case B: The walls are enough, cannot add more walls, must connect
                    elseif wall_count == clue
                        for x in (a, b)
                            if canconnect(game, ds, (i,j), x)
                                changes += connect(game, ds, (i,j), x)
                            end
                        end
                    end
                end
            end
        end
    end
    return changes
end

# Step 7: Fill the remaining regions using backtracking
# The attempt number is limited to avoid infinite loops, so it has possibility to fail
function fill_all(game::PalisadeGame, ds::DisjointSet; max_attempts::Int=20000)
    attempts = Ref(0)
    original_game = deepcopy(game)
    original_ds = deepcopy(ds)

    function is_complete()
        for (rep, cells) in get_region_cells(game, ds)
            if length(cells) != game.region_size
                return false
            end
        end
        return true
    end

    function dfs_fill()
        attempts[] += 1
        if attempts[] > max_attempts
            return false
        end

        # Get all the uncomplete regions
        regions = [(rep, cells) for (rep, cells) in get_region_cells(game, ds) if ds.size[rep] < game.region_size]
        if isempty(regions)
            return true
        end

        # Priority queue, first extend the region with the least number of neighbors
        options = []

        for (rep, cells) in regions
            extendables = Set{Int}()
            candidate = nothing
            for cell in cells
                for neighbor in extendlist(game, ds, cell)
                    nroot = region(game, ds, neighbor...)
                    if nroot != rep
                        push!(extendables, nroot)
                        candidate = (cell, neighbor)
                    end
                end
            end
            if !isempty(extendables)
                push!(options, (length(extendables), rep, candidate))
            end
        end

        if isempty(options)
            return false  # No options available
        end

        # Start from the region with the least number of neighbors
        sort!(options, by = x -> x[1])
        _, _, (a, _) = options[1]

        # Connect its neighboring cells randomly
        for b in shuffle(extendlist(game, ds, a))
            if canconnect(game, ds, a, b)
                game_backup = deepcopy(game)
                ds_backup = deepcopy(ds)

                connect(game, ds, a, b)

                if dfs_fill()
                    return true
                end

                # Backtrack
                game = deepcopy(game_backup)
                ds.parent .= ds_backup.parent
                ds.size .= ds_backup.size
            end
        end

        return false
    end

    success = dfs_fill()

    if !success
        println("fill(): Still not able to divide after $max_attempts attempts, return to original state.")
        game = deepcopy(original_game)
        ds.parent .= original_ds.parent
        ds.size .= original_ds.size
    end

    return success
end

# Main solver function
function solve2(game::PalisadeGame, ds::DisjointSet; verbose::Bool=true)
    total_changes = 0

    # Step 1: Only one time
    step1 = connected_clues_vs_region(game, ds)
    total_changes += step1
    verbose && println("Step 1: connected_clues_vs_region → $step1 changes")

    # Step 2~6: Loop until no changes
    iteration = 0
    while true
        iteration += 1
        changes = 0
        c2 = number_exhausted(game, ds)
        c3 = not_too_big(game, ds)
        c4 = not_too_small(game, ds)
        c5 = no_dangling_edges(game, ds)
        c6 = equivalent_edges(game, ds)
        changes = c2 + c3 + c4 + c5 + c6
        total_changes += changes
        verbose && println("Iteration $iteration: Steps 2~6 → $changes changes")
        if changes == 0
            break
        end
    end

    # Step 7: fill (backtracking for the remaining regions)
    verbose && println("Step 7: Attempting fill...")
    success = fill_all(game, ds)

    if success
        verbose && println("solve: Puzzle successfully solved!")
    else
        verbose && println("solve: Could not finish all regions, partial solution kept.")
    end

    # Visualize
    if verbose
        show_result(game, ds)
    end

    return success
end