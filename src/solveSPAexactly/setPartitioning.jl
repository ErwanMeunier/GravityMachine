# From: https://github.com/vOptSolver/vOptGeneric.jl/blob/master/examples/setPartitionning.jl
# Bi-objective set partitionning problem (2-SPA)


# ---- Packages to use
using JuMP, GLPK


# ---- Parser reading an instance of 2-SPA (format of instances compliant with vOptLib)
function loadInstance2SPA(fname::String)

    f = open(fname)
    nbctr, nbvar = parse.(Int, split(readline(f))) # number of constraints , number of variables
    A = zeros(Int, nbctr, nbvar)                   # matrices of constraints
    c1 = zeros(Int, nbvar)                         # 1st vector of costs
    c2 = zeros(Int, nbvar)                         # 2nd vector of costs
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c1[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                c2[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                A[j,i] = 1
            end
        end
    end
    close(f)
    return c1, c2, A
end

function loadInstance2SPA(fname::String)

    f = open("../../SPA/instances/"*fname)
    nbctr, nbvar = parse.(Int, split(readline(f))) # nombre de contraintes , nombre de variables
    A = zeros(Int, nbctr, nbvar)                   # matrice des contraintes
    c1 = zeros(Int, nbvar)                         # vecteur des couts
    c2 = zeros(Int, nbvar)                         # deuxième vecteur des couts
    nb = zeros(Int, nbvar)
    for i in 1:nbvar
        flag = 1
        for valeur in split(readline(f))
            if flag == 1
                c1[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 2
                c2[i] = parse(Int, valeur)
                flag +=1
            elseif flag == 3
                nb[i] = parse(Int, valeur)
                flag +=1
            else
                j = parse(Int, valeur)
                A[j,i] = 1
            end
        end
    end
    close(f)
    return c1, c2, A
end

function loadNDPoints2SPA(fname::String)
    fname = "../../SPA/Y/Y_N_"*fname
    f = open(fname)
    temps = parse.(Float64, readline(f))
    len = parse.(Int, readline(f))
    xN = Vector{Int64}(undef, len)
    yN = Vector{Int64}(undef, len)
    for i in 1:len
        x, y =  parse.(Float64, split(readline(f)))
        xN[i],yN[i] = round(Int,x),round(Int,y)
    end
    close(f)
    return xN, yN
end

#=
# ---- Compute the set of non-dominated points Y_N of a 2SPA with vOptGeneric
function computeYNfor2SPA(  nbvar::Int,
                            nbctr::Int,
                            A::Array{Int,2},
                            c1::Array{Int,1},
                            c2::Array{Int,1}
                         )

    # ---- setting the model
    model = vModel( GLPK.Optimizer )
    JuMP.set_silent(model)

    @variable(model, x[1:nbvar], Bin)
    @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    @addobjective(model, Min, sum(c1[i]*x[i] for i in 1:nbvar))
    @addobjective(model, Min, sum(c2[i]*x[i] for i in 1:nbvar))

    # ---- Invoking the solver (epsilon constraint method)
    vSolve(model, method=:epsilon, step = 0.5)

    # ---- Querying the results
    start = time()
    Y_N = getY_N(model)
    elapsedTime = time() - start
    return Y_N, elapsedTime
end


# ---- Save an instance of 2-SPA (format of instances compliant with vOptLib)
function saveInstance2SPA(fname::String, Y_N, elapsedTime::Float64)
    way = "Y_N_"*fname
    open(way, "w") do f
        write(f, "$elapsedTime \n")
        write(f, "$(length(Y_N)) \n")
        for i in 1:length(Y_N)
            write(f, "$(Y_N[i][1]) $(Y_N[i][2]) \n")
        end
    end
end


# ---- Entry point
function main(fname::String)

    print("Compute Y_N with vOPtGeneric for a set partitionning problem with 2 objectives (2-SPA) \n\n")

    # load a numerical instance of 2SPA ----------------------------------------
    c1, c2, A = loadInstance2SPA(fname)
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    # compute the set of non-dominated points Y_N ------------------------------
    Y_N, elapsedTime = computeYNfor2SPA(nbvar, nbctr, A, c1, c2)

    # save the result on a file
    saveInstance2SPA(fname, Y_N, elapsedTime)
end
=#