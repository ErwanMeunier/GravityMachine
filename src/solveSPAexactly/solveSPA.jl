include("setPartitioning.jl") # parser and other tools 

using Distributed

@everywhere using JuMP
@everywhere import HiGHS
@everywhere import MultiObjectiveAlgorithms as MOA

@everywhere using DataFrames, CSV

const instancesPath = "../../SPA/instances/"
const verbose = false 

# TIMING OUT ------------------------------------------------------------
# Inspired from: https://discourse.julialang.org/t/help-writing-a-timeout-macro/16591/6
function timeout(f, arg, seconds, fail)
    tsk = @task f(arg...)
    schedule(tsk)
    Timer(seconds) do timer
        istaskdone(tsk) || Base.throwto(tsk, InterruptException())
    end
    try
        fetch(tsk)
    catch _;
        fail
    end
end
#------------------------------------------------------------------------

# JuMP Model and solving
# https://voptsolver.github.io/vOptLib/SPA/spa2021.html

function solveBiSPA!(c1::Vector{Int}, c2::Vector{Int}, A::Array{Int,2}, nbsolutionsLimit::Int64=5, timeLimit::Float64=180) #ok
    return @elapsed begin
        m = size(A,1) # number of columns/constraints
        n = size(A,2) # number of rows/variables

        model = JuMP.Model(()->MOA.Optimizer(HiGHS.Optimizer))
        set_time_limit_sec(model, timeLimit)
        set_optimizer_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
        set_optimizer_attribute(model, MOA.SolutionLimit(), nbsolutionsLimit)

        @variable(model, x[1:n], Bin)
        @expression(model, obj1, sum(x[j]*c1[j] for j=1:n))
        @expression(model, obj2, sum(x[j]*c2[j] for j=1:n))
        @objective(model, Min, [obj1,obj2])
        @constraint(model, [i=1:m], (sum((x[j]*A[i,j]) for j in 1:n)) == 1)

        optimize!(model)

        verbose ? solution_summary(model) : nothing
    end
end

function main()
    timeLimit = 500.
    nbsolLimit = 10
    instances::Vector{String} = readdir(instancesPath)[5:end]
    times::Vector{Float64} = fill(-1., length(instances)) # none benchmarked instance are clearly discriminated
    
    for i in eachindex(instances)
        # would be better to put a try catch end around here
        file = instances[i]
        c1, c2, A = loadInstance2SPA(instancesPath*file)
        times[i] = solveBiSPA!(c1,c2,A,nbsolLimit,timeLimit)
        println(times)
        println("Parsing is OK")
    end
    output = DataFrame([instances,times],[:Instances,:Times])
    CSV.write("../results/timeExactSolving.csv", output)
end
