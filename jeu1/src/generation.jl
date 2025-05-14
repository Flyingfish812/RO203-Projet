# This file contains methods to generate a data set of instances (i.e., sudoku grids)
include("io.jl")
using Random
using Glob

"""
Generate an n*n grid with a given density

Argument
- n: size of the grid
- density: percentage in [0, 1] of initial values in the grid
"""

# check the bounds of the CartesianIndex
inbounds(A::AbstractArray, idx::CartesianIndex) = all(1 .<= Tuple(idx) .<= size(A))

# Divide the grid into several regions of region_size connexes (using heuristic flood-fill)
function divide_grid_into_regions(w::Int, h::Int, region_size::Int)::Array{Int,2}
    region_ids = fill(-1, h, w)
    region_id = 1
    remaining = vec(collect(CartesianIndices((1:h, 1:w))))
    visited = falses(h, w)
    retry_count = 0
    max_retries = 2000

    while !isempty(remaining)
        retry_count += 1
        if retry_count > max_retries
            error("divide_grid_into_regions: surpass the maximum number of retries")
        end

        start = rand(remaining)
        region = Set{CartesianIndex}()
        frontier = [start]

        while !isempty(frontier) && length(region) < region_size
            current = pop!(frontier)
            if current in region || visited[current]
                continue
            end

            push!(region, current)
            visited[current] = true

            for offset in [(0,1), (1,0), (0,-1), (-1,0)]
                nei = current + CartesianIndex(offset)
                if inbounds(region_ids, nei) && !visited[nei] && !(nei in region)
                    push!(frontier, nei)
                end
            end
        end

        if length(region) == region_size
            for idx in region
                region_ids[idx] = region_id
                idx_pos = findfirst(isequal(idx), remaining)
                if idx_pos !== nothing
                    deleteat!(remaining, idx_pos)
                end
            end
            region_id += 1
        else
            # Not big enough, rollback, mark as unvisited
            for idx in region
                visited[idx] = false
            end
        end
    end

    return region_ids
end

# Generate an instance of the game and save it to a file
function generateInstance(m::Int64, n::Int64, k::Int64, density::Float64, idx::Int64)
    @assert (m * n) % k == 0 "Grid area m*n must be divisible by region size k"

    region_ids = divide_grid_into_regions(n, m, k)  # width: n, height: m

    # Initialize the walls
    h_walls = fill(Int8(0), m + 1, n)
    v_walls = fill(Int8(0), m, n + 1)

    h_walls[1, :] .= 1
    h_walls[end, :] .= 1
    v_walls[:, 1] .= 1
    v_walls[:, end] .= 1

    # Internal walls
    for i in 1:m
        for j in 1:n
            if i < m && region_ids[i, j] != region_ids[i+1, j]
                h_walls[i+1, j] = 1
            end
            if j < n && region_ids[i, j] != region_ids[i, j+1]
                v_walls[i, j+1] = 1
            end
        end
    end

    # Calculate the clues
    clues = Dict{Tuple{Int,Int}, Int}()
    for i in 1:m, j in 1:n
        cnt = 0
        cnt += h_walls[i, j] == 1 ? 1 : 0
        cnt += h_walls[i+1, j] == 1 ? 1 : 0
        cnt += v_walls[i, j] == 1 ? 1 : 0
        cnt += v_walls[i, j+1] == 1 ? 1 : 0
        clues[(i, j)] = cnt
    end

    # Keep only a fraction of the clues based on the density
    all_positions = collect(keys(clues))
    keep = rand(all_positions, round(Int, length(all_positions) * density))
    filtered_clues = Dict(k => clues[k] for k in keep)

    # Output to data, with idx
    filename = "data/generated_$(m)x$(n)_k$(k)_d$(Int(round(density*100)))_id$(idx).txt"
    open(filename, "w") do io
        println(io, "width $n")
        println(io, "height $m")
        println(io, "region $k")
        println(io, "clues")
        for ((i,j), v) in sort(collect(filtered_clues))
            println(io, "$i $j $v")
        end
        println(io, "end")
    end

    println("Generated game: $filename")
end

# Full dataset
function generateDataSet(clear_data::Bool=true)
    base_seed = 1234
    Random.seed!(base_seed)  # Fix the seed for reproducibility
    seed_offset = 0

    if clear_data
        for f in Glob.glob("data/generated_*.txt")
            rm(f, force=true)
        end
        println("Cleared all previous data files.")
    end

    sizes = [(5,5,5), (8,6,6), (10,8,8), (15,12,10)]
    densities = [0.2, 0.4, 0.5]
    count = 10

    for (m, n, k) in sizes
        for d in densities
            for idx in 1:count
                success = false
                while !success
                    try
                        Random.seed!(base_seed + seed_offset)
                        generateInstance(m, n, k, d, idx)
                        success = true
                    catch e
                        seed_offset += 1  # Use a different seed
                    end
                end
                seed_offset += 1  # Ensure that the next instance is different
            end
        end
    end
end

# -------------- Test case and functions ----------------------------

function create_test_case()
    width, height, region_size = 5, 5, 5

    # Initialize the clues
    clues = Dict{Tuple{Int,Int}, Int}()
    clues[(1,1)] = 3
    clues[(1,2)] = 3
    clues[(1,3)] = 3
    clues[(5,1)] = 2
    clues[(5,2)] = 1

    # Initialize the walls
    h_walls = fill(Int8(-1), height + 1, width)
    v_walls = fill(Int8(-1), height, width + 1)
    h_walls[1, :] .= 1
    h_walls[end, :] .= 1
    v_walls[:, 1] .= 1
    v_walls[:, end] .= 1

    game = PalisadeGame(width, height, region_size, clues, h_walls, v_walls)
    ds = DisjointSet(width * height)

    return game, ds
end

function test_function(game::PalisadeGame, ds::DisjointSet)
    println("=== parent and size ===")
    println("parent: ", ds.parent)
    println("size:   ", ds.size)

    println("\n=== flatten_coord / unflatten_coord ===")
    pos = (4,2)
    flat = flatten_coord(game, pos...)
    println("flatten_coord(4,2) = ", flat)
    println("unflatten_coord($flat) = ", unflatten_coord(game, flat))

    println("\n=== region ===")
    for pos in [(4,2), (5,1), (5,2), (5,3)]
        println("region(", pos, ") = ", region(game, ds, pos...))
    end

    println("\n=== connected ===")
    println("connected((4,2), (5,1)) = ", connected(game, ds, (4,2), (5,1)))
    println("connected((4,2), (5,3)) = ", connected(game, ds, (4,2), (5,3)))
    println("connected((4,2), (1,1)) = ", connected(game, ds, (4,2), (1,1)))  # not connected

    println("\n=== wall() ===")
    for pos in [(1,1), (1,2), (5,2)]
        println("wall(", pos, ") = ", wall(game, pos...))
    end

    println("\n=== get_wall_between() ===")
    println("get_wall_between((1,1), (1,2)) = ", get_wall_between(game, (1,1), (1,2)))  # with wall
    println("get_wall_between((5,1), (5,2)) = ", get_wall_between(game, (5,1), (5,2)))  # without wall

    println("\n=== get_wall() ===")
    dir, r, c = get_wall_between(game, (1,1), (1,2))
    println("get_wall((1,1)-(1,2)) = ", get_wall(game, dir, r, c))  # 1

    dir, r, c = get_wall_between(game, (5,1), (5,2))
    println("get_wall((5,1)-(5,2)) = ", get_wall(game, dir, r, c))  # 0

    println("\n=== neighbors() ===")
    println("neighbors((1,1)) = ", neighbors(game, (1,1)))
    println("neighbors((3,3)) = ", neighbors(game, (3,3)))  # center

    println("\n=== neighbor() ===")
    println("neighbor((1,1)) = ", neighbor(game, ds, (1,1)))
    println("neighbor((5,2)) = ", neighbor(game, ds, (5,2)))

    println("\n=== boarder() ===")
    println("boarder((1,1), (1,2)) = ", boarder(game, ds, (1,1), (1,2)))  # should have wall
    println("boarder((5,1), (5,2)) = ", boarder(game, ds, (5,1), (5,2)))  # same region without boarder

    println("\n=== canconnect() ===")
    println("canconnect((1,1), (1,2)) = ", canconnect(game, ds, (1,1), (1,2)))  # false, have wall
    println("canconnect((3,3), (4,3)) = ", canconnect(game, ds, (3,3), (4,3)))  # true
    println("canconnect((4,2), (4,3)) = ", canconnect(game, ds, (4,2), (4,3)))  # true 

    println("\n=== extendlist() ===")
    println("extendlist((4,2)) = ", extendlist(game, ds, (4,2)))
    println("extendlist((1,1)) = ", extendlist(game, ds, (1,1)))  # blocked

    println("\n=== disconnected() ===")
    println("disconnected((1,1), (1,2)) = ", disconnected(game, ds, (1,1), (1,2)))
    println("disconnected((5,1), (5,2)) = ", disconnected(game, ds, (5,1), (5,2)))  # false
end

function show_result(game::PalisadeGame, ds::DisjointSet)
    println("=== Display the game status ===")

    # Grids
    if @isdefined displayGrid
        displayGrid(game)
    end

    println("\n=== Disjoint set parent ===")
    for i in 1:game.height
        for j in 1:game.width
            idx = flatten_coord(game, i, j)
            print(rpad(ds.parent[idx], 4))
        end
        println()
    end

    println("\n=== Disjoint set size ===")
    for i in 1:game.height
        for j in 1:game.width
            idx = flatten_coord(game, i, j)
            print(rpad(ds.size[idx], 4))
        end
        println()
    end

    println("\n=== Walls ===")
    # Visualize walls
    function wallchar(val)
        val == 1 && return "1"
        val == 0 && return "0"
        return "-"  # Unknown -1
    end

    h, w = game.height, game.width

    for i in 1:h+1
        # Horizontal walls
        print(" ")
        for j in 1:w
            print(" ", wallchar(game.h_walls[i, j]), "  ")
        end
        println()

        if i <= h
            # Vertical walls
            for j in 1:w+1
                print(wallchar(game.v_walls[i, j]), "   ")
            end
            println()
        end
    end
end