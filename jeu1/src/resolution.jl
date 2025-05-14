# This file contains methods to solve an instance (heuristically or with CPLEX)
using CPLEX
using Dates
using JSON3
using LinearAlgebra

include("generation.jl")
include("solver1.jl")
include("solver2.jl")

TOL = 0.00001

"""
Solve an instance with CPLEX
"""
function cplexSolve(game::PalisadeGame; verbose::Bool = true)
    original_game = deepcopy(game)
    m = build_model(original_game)

    if verbose
        write_to_file(m, "model.lp")
    end

    start = time()
    optimize!(m)
    elapsed = time() - start

    if verbose
        ds = DisjointSet(game.height * game.width)
        extractSolution(original_game, ds, m)
        show_result(original_game, ds)
    end

    return JuMP.primal_status(m) == JuMP.MOI.FEASIBLE_POINT, elapsed
end

"""
Heuristically solve an instance
"""

function heuristicSolve(game::PalisadeGame; verbose::Bool = true)
    original_game = deepcopy(game)
    ds = DisjointSet(game.height * game.width)

    start = time()
    success = solve2(original_game, ds; verbose = verbose)
    elapsed = time() - start

    return success, elapsed
end

function solveOneData(file::String; method::String, verbose::Bool = true)
    # Read one instance
    game = readInputFile(file)
    if verbose
        println("Solving instance: ", file)
    end

    # Call for the solver according to the method
    if method == "cplex"
        isOptimal, solveTime = cplexSolve(game; verbose=verbose)
    elseif method == "heuristique" || method == "heuristic"
        isOptimal, solveTime = heuristicSolve(game; verbose=verbose)
    else
        error("Unknown method: $method")
    end
    if verbose
        println("Results: isOptimal = ", isOptimal, ", solveTime = ", solveTime, " s")
    end

    return isOptimal, solveTime
end

"""
Solve all the instances contained in "../data" through CPLEX and heuristics

The results are written in "../res/cplex" and "../res/heuristic"

Remark: If an instance has previously been solved (either by cplex or the heuristic) it will not be solved again
"""
function solveDataSet()
    dataFolder = "./data/"
    resFolder = "./res/"

    # Resolution method: ["cplex"] or ["cplex", "heuristique"]
    resolutionMethods = ["cplex", "heuristique"]

    # Output folders
    resolutionFolders = [joinpath(resFolder, meth) for meth in resolutionMethods]
    for folder in resolutionFolders
        if !isdir(folder)
            mkdir(folder)
        end
    end

    # For all files in dataFolder ended with ".txt"
    files = filter(f -> endswith(f, ".txt"), readdir(dataFolder))
    for file in files
        dataFilePath = joinpath(dataFolder, file)
        println("===========================================")
        println("Solving: ", file)

        # Use each method
        for (idx, meth) in enumerate(resolutionMethods)
            outputFile = joinpath(resolutionFolders[idx], file)
            # If file exists, skip
            if isfile(outputFile)
                println("File ", outputFile, " already exists, skipping...")
                continue
            end

            println("Method: ", meth)
            isOpt, timeUsed = solveOneData(dataFilePath; method=meth, verbose=false)
            # 这里可以进一步写入输出文件中，比如写入求解时刻、求解耗时以及部分游戏状态信息（例如墙体、区域划分等）
            open(outputFile, "w") do fout
                println(fout, "solveTime = ", timeUsed)
                println(fout, "isOptimal = ", isOpt)
            end
            println(meth, " results have been written to: ", outputFile)
        end
    end
end