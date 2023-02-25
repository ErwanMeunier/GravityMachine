include("setPartitioning.jl") # parser and other tools 

using Distributed

@everywhere using JuMP
@everywhere import HiGHS
@everywhere import MultiObjectiveAlgorithms as MOA

@everywhere using DataFrames, CSV

const instancesPath = "../../SPA/instances/"
const verbose = false 

# TIMING OUT ------------------------------------------------------------
# https://discourse.julialang.org/t/help-writing-a-timeout-macro/16591/6
macro timeout(seconds, expr, fail)
    quote
        tsk = @task $expr
        schedule(tsk)
        Timer($seconds) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        #=try
            println(fetch(tsk))
        catch InterruptException
            $fail
        end=#
        fetch(tsk)
    end
end
#------------------------------------------------------------------------

# JuMP Model and solving
# https://voptsolver.github.io/vOptLib/SPA/spa2021.html

function solveBiSPA!(c1::Vector{Int}, c2::Vector{Int}, A::Array{Int,2}, nbsolutionsLimit::Int64=5) #ok
    m = size(A,1) # number of columns/constraints
    n = size(A,2) # number of rows/variables

    model = JuMP.Model(()->MOA.Optimizer(HiGHS.Optimizer))
    set_optimizer_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())
    set_optimizer_attribute(model, MOA.SolutionLimit(), nbsolutionsLimit)

    @variable(model, x[1:n], Bin)
    @expression(model, obj1, sum(x[j]*c1[j] for j=1:n))
    @expression(model, obj2, sum(x[j]*c2[j] for j=1:n))
    @objective(model, Min, [obj1,obj2])
    @constraint(model, [i=1:m], (sum((x[j]*A[i,j]) for j in 1:n)) == 1)

    optimize!(model)

    verbose ? solution_summary(model) : nothing
    return true
end

function main()
    timeLimit = 10.
    nbsolLimit = 5
    instances::Vector{String} = readdir(instancesPath)[1:1]
    times::Vector{Float64} = fill(-1., length(instances)) # none benchmarked instance are clearly discriminated
    
    for i in eachindex(instances)
        # would be better to put a try catch end around here
        file = instances[i]
        #try 
        println(pwd())
        c1, c2, A = loadInstance2SPA(instancesPath*file)
        println("parsing OK")
        println(timeLimit)
        solveBiSPA!(c1,c2,A,times,i,nbsolLimit)
        println(state)
        #catch SystemError;
        #    println("WARNING : "*instancesPath*file*" cannot be opened.")
        #finally
        #    System.out.println("Error not handled -> exit(2) !!!")
        #    exit(2)
        #    nothing
        #end
    end
    output = DataFrame([instances,times],[:Instances,:Times])
    CSV.write("../results/timeExactSolving.csv", output)
end

f(x) = begin 2*x ; println(2*x) ; sleep(1.) end

function test()
    
    #println("Testing f: f(3)")
    for i in 1:10
        @timeout 10 f(i) "fail"
        #println(i)
    end
end