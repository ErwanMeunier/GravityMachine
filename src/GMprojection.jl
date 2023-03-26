const MINIMAL_NUMBER_OF_ARGUMENTS_PROJ::Int = 2 # must be fixed by the programmer
# Tools
function computeCλ(k::Int,λ1::Vector{Float64},λ2::Float64,c1::Float64,c2::Float64, guided::Bool)::Vector{Float64}
    guided ? (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2 : ones(length(c1))
end

# ============================Some tools================================================
# Self-explanatory
function computeLocalNadirs(L::Vector{tSolution{Float64}})::Vector{tPoint}
    @assert length(L)>2 "vg must be of size greater or equal than 3"
    nadirs::Vector{tPoint} = Vector{tPoint}(undef, length(L))
    nadirs[1] = tPoint(L[end].y[1],L[1].y[2]) # global nadir
    nadirs[end] = nadirs[1]
    for k = 2:(length(L)-1)
        nadirs[k] = tPoint(L[k+1].sRel.y[1], L[k-1].sRel.y[2])
    end
    return nadirs
end

# ==============================PROJECTION METHODS======================= 
#=
Any new projection method must fill the following specifications:
- Their arguments follow the order of the method with the larger number of arguments i.e. (A::Array{Int,2}, xTilde::Array{Int,1}, c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64}, etc)
- The method is referenced in the mapping "configurationProjection::Dict{Int,Tuple{Function,Int}}"
- The field of integer type is dedicated to the NUMBER OF ARGUMENTS needed for the added method
- TODO
=#
# =======================================================================

# Projete xTilde sur le polyedre X du SPA avec norme-L1
# version FP 2005

function Δ2SPA(A::Array{Int,2}, xTilde::Array{Int,1})

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
    @objective(proj, Min, sum(x[i] for i in idxTilde0) + sum((1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

# ==============================================================================
# Projete xTilde sur le polyedre X du SPA avec norme-L1
# version avec somme ponderee donnant la direction vers le generateur k
# α ∈ [0,1] -> Module l'utilisation de la projection guidée

function Δ2SPAbis(A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64,
                  c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2
    #cλ = fill(1.,length(xTilde))
    #cλ = 1.0 .+ α.*(λ1[k].*c1 + λ2[k].*c2  .-1)
    #println("cλ :", cλ)
    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
#    @objective(proj, Min, sum(λ1[k]*x[i] for i in idxTilde0) + sum(λ2[k]*(1-x[i]) for i in idxTilde1) )
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    optimize!(proj)
    return objective_value(proj), value.(x)
end

#TODO: factorizing Δ2SPAbisCone arguments
function Δ2SPAbisCone(L::Vector{tSolution{Float64}},A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64,
                  c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, nadirs::Vector{tPoint}, vg::Vector{tGenerateur}, α::Float64=1., β::Float64=0.5)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)
    gN = tPoint(L[end].y[1],L[1].y[2]) # global nadir point
    lN = tPoint(nadirs[k].x+ β*(gN.x-nadirs[k].x),nadirs[k].y+ β*(gN.x-nadirs[k].y)) # local nadir point 
    cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2
    #cλ = 1.0 .+ α.*(λ1[k].*c1 + λ2[k].*c2  .-1)
    #println("cλ :", cλ)
    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0 )
#   @objective(proj, Min, sum(λ1[k]*x[i] for i in idxTilde0) + sum(λ2[k]*(1-x[i]) for i in idxTilde1) )
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1) )
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    # Staying in the cone
    if (k!=1) & (k!=length(vg)) 
        @expression(proj, z1, sum(c1[i]*x[i] for i in 1:nbvar))
        @expression(proj, z2, sum(c2[i]*x[i] for i in 1:nbvar))
        # vg are in decreasing order considering their second coordinate
        # Hence, lines delimiting the cone are computed using:
        # The k-1th and the k+1th generators
        u = tPoint(vg[k-1].sRel.y[1],vg[k-1].sRel.y[2])
        v = tPoint(vg[k+1].sRel.y[1],vg[k+1].sRel.y[2])
        #--- Equation delimiting the current cone of research:
        # first line equation
        a1 = (lN.y - u.y)/(lN.x - u.x)
        b1 = u.y- a1*u.x
        # second line equation
        a2 = (lN.y - v.y)/(lN.x - v.x)
        b2 = v.y- a2*v.x
        #--- 
        println("k=",k,"/",length(vg))
        println("First equation : ", a1, "*x+",b1)
        println("Second equation : ", a2, "*x+",b2)
        # building up the cone
        if (k!=1) & (k!=length(vg))
            if (β==0.) # lines which delimits the current cone are respectively vertical and horizontal
                println("β=0 => Vertical and horizontal constraints")
                @constraint(proj, z2<= u.y)
                @constraint(proj, z1<= v.x)
            elseif (k-1)==1 # upmost cone
                println("UPMOST CONE")
                @constraint(proj, z2 <= u.y)
                @constraint(proj, z2 >= a2*z1 + b2)
            elseif (k+1)==length(vg) # lowest cone
                println("LOWEST CONE")
                @constraint(proj, z2 <= a1*z1 + b1)
                @constraint(proj, z2 >= v.y)
            else # inner cone
                println("INNER CONE")
                @constraint(proj, z2 <= a1*z1 + b1)
                @constraint(proj, z2 >= a2*z1 + b2)
            end   
        end
    end
    # ---
    optimize!(proj)
    return objective_value(proj), value.(x)
end

function Δ2SPABelgique(A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64,
    c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64=1.)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2

    proj::Model = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0)
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1)) 
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)    
    optimize!(proj)

    xOutput = value.(x)
    [set_binary(proj[:x][i]) for i=1:length(xTilde) if !(isapprox(xOutput[i],0,atol=10^-3)||isapprox(xOutput[i],1,atol=10^-3))]
    optimize!(proj)

    return objective_value(proj), value.(x)
end


#=
knownSol : Already known solution
nbCycles : Number of times where the projection reach an already known solution

=#
function Δ2SPABelgique2(A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64,
    c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64, solutionsHist::Set{Vector{Int}}, knownSol::Vector{Int}, nbCycles::Int)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2

    proj::Model = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0)
    if nbCycles > 0
        idx0attrac, idx1attrac = split01(knownSol)
        @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1)
                                +nbCycles*(sum(cλ[i]*x[i] for i in idx1attrac)+sum(cλ[i]*(1-x[i]) for i in idx0attrac))
                  )
    else
        @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1)) 
    end
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)    
    # diff constraint 
    for sol in solutionsHist
        N0, N1 = split01(sol)
        @constraint(proj, sum(x[j] for j in N0) + sum(1-x[j] for j in N1) >= 1) # diff from knownSol does not work
    end
    # end - diff constraint
    optimize!(proj)

    xOutput = value.(x)
    [set_binary(proj[:x][i]) for i=1:length(xTilde) if !(isapprox(xOutput[i],0,atol=10^-3)||isapprox(xOutput[i],1,atol=10^-3))]
    
    #=
    if nbCycles > 0 
        println("Modifying the model used in the projection in order to avoid reaching the known solution twice")
        #xOutput = value.(x)
        idxTilde0KnownSol, idxTilde1KnownSol = split01(knownSol)
        #[set_binary(proj[:x][i]) for i in 1:nbvar]
        @constraint(proj, diffToKnownSolution, sum(x[i] for i in idxTilde0KnownSol)+sum(x[i] for i in idxTilde1KnownSol)<=nbvar-1)
    end
    =#
    optimize!(proj)
    

    return objective_value(proj), value.(x)
end


# ==============================================================================

# projecte la solution entiere correspondant au generateur k et test d'admissibilite
# (A::Array{Int,2}, xTilde::Array{Int,1}, 
# c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, protectedIndexOfInt::Vector{Int64}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)
function projectingSolution!(A::Array{Int,2}, vg::Vector{tGenerateur}, k::Int, c1::Vector{Int}, c2::Vector{Int}, d::tListDisplay, args::Vararg{Any}=Tuple{Any}())
    println("---Projection solution--- | generator ", k)
    # --------------------------------------------------------------------------
    # Projete la solution entiere sur le polytope X 
    # generalNadir = fill(tPoint(L[end].y[1],L[1].y[2]),length(vg))
    # nadirs = computeLocalNadirs(vg,L)
    fPrj, vg[k].sPrj.x = interface_projection!(A,vg[k].sInt.x,k,c1,c2,args...) #: Δ2SPAbis(A,vg[k].sInt.x,c1,c2,k,λ1,λ2,α)
    
    # Nettoyage de la valeur de vg[k].sPrj.x et calcul du point bi-objectif
    # reconditionne les valeurs 0 et 1 et arrondi les autres valeurs
    nettoyageSolution!(vg[k].sPrj.x)
#    verbose ? @printf("  %2dP : fPrj = %8.2f  ",k, round(fPrj, digits=2)) : nothing

    # recalcule la solution au regard des 2 objectifs
    vg[k].sPrj.y[1], vg[k].sPrj.y[2] = evaluerSolution(vg[k].sPrj.x, c1, c2)
    verbose ? @printf("  %2dP : [ %8.2f , %8.2f ] ",k, vg[k].sPrj.y[1], vg[k].sPrj.y[2]) : nothing

    # archive le point obtenu pour les besoins d'affichage
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XProj, vg[k].sPrj.y[1])
        push!(d.YProj, vg[k].sPrj.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XProj, vg[k].sPrj.y[1])
        push!(d.YProj, vg[k].sPrj.y[2])
    end            

    # ----------------------------------------------------------------
    # Teste si la projection est admissible

    if estAdmissible(vg[k].sPrj.x)
        # sauvegarde de la solution entiere admissible obtenue
        # vg[k].sInt.x = deepcopy(vg[k].sPrj.x)
        # vg[k].sInt.y[1] = vg[k].sPrj.y[1]
        # vg[k].sInt.y[2] = vg[k].sPrj.y[2]
        vg[k].sFea = true
        @printf("→ Admissible "); print("                       ")
        # ajouterXtilde!(vg, k, convert.(Int, round.(vg[k].sPrj.x)), convert.(Int, round.(vg[k].sPrj.y))) # TODO: GET RID OF THIS LINE 
        vg[k].sPrj.x = round.(vg[k].sPrj.x)
        vg[k].sPrj.y = round.(vg[k].sPrj.y)
        # archive le point obtenu pour les besoins d'affichage
        if generateurVisualise == -1 
            # archivage pour tous les generateurs
            push!(d.XFeas, copy(vg[k].sPrj.y[1]))
            push!(d.YFeas, copy(vg[k].sPrj.y[2]))
        elseif generateurVisualise == k
            # archivage seulement pour le generateur k
            push!(d.XFeas, copy(vg[k].sPrj.y[1]))
            push!(d.YFeas, copy(vg[k].sPrj.y[2]))
        end  
    else

        vg[k].sFea = false
        @printf("→ x          "); print("                       ")
        # prepare pour l'iteration suivante
#        vg[k].xRlx = deepcopy(vg[k].sPrj.x) !!!!!!!!!!!!!
    end
end

# =================================================PROJECTION INTERFACE=======================================================================

const configurationProjection::Dict{Int,Tuple{Function,Int}} = Dict{Int,Tuple{Function,Int}}( # signature of each methods input format
                                                                1 => (Δ2SPA,2), # (A::Array{Int,2}, xTilde::Array{Int,1})
                                                                2 => (Δ2SPAbis,8), # (A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)
                                                                3 => (Δ2SPAbisCone,11), # (A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64, nadirs::Vector{tPoint}, vg::Vector{tGenerateur}, β::Float64)
                                                                4 => (Δ2SPABelgique,8), # (A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)
                                                                5 => (Δ2SPABelgique2,11) # (A::Array{Int,2}, xTilde::Array{Int,1}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64, solutionsHist::Set{Vector{Int}, knownSol, nbCycles})
                                                            )

# args are in fact additional arguments
function interface_projection!(A::Array{Int,2},xTilde::Vector{Int},args::Vararg{Any}=Typle{Any}();CHOICE::Int=CHOICE_PROJECTION) 
    @assert configurationProjection[CHOICE][2] <= (MINIMAL_NUMBER_OF_ARGUMENTS_PROJ + length(args)) "Bad number of arguments for such a projection"
    return configurationProjection[CHOICE][1](A,xTilde,args[1:(configurationProjection[CHOICE][2]-MINIMAL_NUMBER_OF_ARGUMENTS_PROJ)]...)
end