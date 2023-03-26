include("GMparsers.jl")
#include("GMprojection.jl")
using PyPlot ; pygui(true)
using HiGHS, JuMP, MathOptInterface, LaTeXStrings

# n is the size of vectors returned
function generateRandomBinSolution(n,nbsol)::Vector{Vector{Int}}
    result::Vector{Vector{Int}} = [Vector{Int}(undef,n) for k in 1:nbsol]
    sizehint!(result,nbsol)
    for k in 1:nbsol
        for i in 1:n
            result[k][i] = rand() <= 1/2 ? 0 : 1
        end
    end
    return result
end

function modelFeasibility(A::Matrix{Int})
    m::Int64 = size(A,1)
    n::Int64 = size(A,2)
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:n], Bin)
    @objective(model, Min, 0)
    @constraint(model, [i=1:m], sum(x[j]*A[i,j] for j in 1:n)==1)

    return model
end

function testFeasibility(model::Model, sol::Vector{Int})::Bool
    n = length(sol)
    @assert num_variables(model)==n "Not the same number of variables"
    variables = all_variables(model)
    for i in 1:n
        fix(variables[i],sol[i];force=true)
    end
    optimize!(model)
    status = termination_status(model)
    # Reset the model removing the constraint
    for i in 1:n
        unfix(variables[i])
    end
    #unregister(model, Symbol("fixedVar"))
    return status == MOI.OPTIMAL
end

#=
pointsDensity must be between 0 and 1. 
The number of points generated equals to ~pointsDensity*2^n with n the length of the vector
=#
function studyDensity(fname::String, nbsol::Int)
    c1, c2, A = loadInstance2SPA(fname)
    n = length(c1)
    println("n=", n)
    XN, YN = loadNDPoints2SPA(fname)
    plot(XN, YN, color="black", linewidth=0.75, marker="+", markersize=1.0, linestyle=":", label = L"y \in Y_N")
    scatter(XN, YN, color="black", marker="+")
    randomBinSolutions::Vector{Vector{Int}} = unique(generateRandomBinSolution(n,nbsol)) # 
    #println(randomBinSolutions)
    model = modelFeasibility(A)
    feasibleSolutions::Vector{Vector{Int}} = [v for v in randomBinSolutions if  testFeasibility(model, v)]
    z1::Vector{Int} = [sum(c1 .* v) for v in feasibleSolutions]
    z2::Vector{Int} = [sum(c2 .* v) for v in feasibleSolutions]
    for v in feasibleSolutions
        scatter(z1, z2)
    end
    show()
end

function studyProjection(fname::String, nbSamples::Int)
    c1, c2, A = loadInstance2SPA(fname)
    n = length(c1)
end

function main()
    studyDensity("sppnw10.txt", 1000)
end

main()