# Tools
# ==============================================================================
# Self-explanatory
function computeLocalNadirs(vg::Vector{tGenerateur}, L::Vector{tSolution{Float64}})::Vector{tPoint}
    @assert length(vg)>2 "vg must be of size greater or equal than 3"
    nadirs::Vector{tPoint} = Vector{tPoint}(undef, length(vg))
    nadirs[1] = tPoint(L[end].y[1],L[1].y[2]) # global nadir
    nadirs[end] = nadirs[1]
    for k = 2:(length(vg)-1)
        nadirs[k] = tPoint(max(vg[k-1].sRel.y[1],vg[k+1].sRel.y[1]),max(vg[k-1].sRel.y[2],vg[k+1].sRel.y[2]))
    end
    return nadirs
end

# ==============================================================================
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
function Δ2SPAbis(A::Array{Int,2}, xTilde::Array{Int,1}, 
                  c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    #cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2
    cλ = fill(1.,length(xTilde))
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

function Δ2SPAbisCone(L::Vector{tSolution{Float64}} ,A::Array{Int,2}, xTilde::Array{Int,1}, 
                  c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64}, nadirs::Vector{tPoint}, vg::Vector{tGenerateur}, α::Float64=1., β::Float64=0.5)

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

function Δ2SPAbisKeepInt(A::Array{Int,2}, xTilde::Array{Int,1}, 
    c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, protectedIndexOfInt::Vector{Int64}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)
    # α = 1
    #cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2
    #cλ = 1.0 .+ α.*(λ1[k].*c1 + λ2[k].*c2  .-1)
    #println("cλ :", cλ)
    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0)
    #    @objective(proj, Min, sum(λ1[k]*x[i] for i in idxTilde0) + sum(λ2[k]*(1-x[i]) for i in idxTilde1) )
    #@objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1))
    @objective(proj, Min, sum(x[i] for i in idxTilde0) + sum((1-x[i]) for i in idxTilde1)) 
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
    # Protecting components which are still integer
    @constraint(proj, [i in protectedIndexOfInt], x[i]==xTilde[i])
    
    #
    optimize!(proj)
    return objective_value(proj), value.(x)
end

function Δ2SPABelgique(A::Array{Int,2}, xTilde::Array{Int,1}, 
    c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64=1.)

    nbctr = size(A,1)
    nbvar = size(A,2)
    idxTilde0, idxTilde1 = split01(xTilde)

    cλ = (1+α*(λ1[k]-1)).*c1 + (1+α*(λ2[k]-1)).*c2

    proj = Model(GLPK.Optimizer)
    @variable(proj, 0.0 <= x[1:length(xTilde)] <= 1.0)
    @objective(proj, Min, sum(cλ[i]*x[i] for i in idxTilde0) + sum(cλ[i]*(1-x[i]) for i in idxTilde1)) 
    @constraint(proj, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)    
    optimize!(proj)

    xOutput = value.(x)
    [set_binary(projInt[:x][i]) for i=1:length(xTilde) if !(isapprox(xOutput[i],0,atol=10^-3)||isapprox(xOutput[i],1,atol=10^-3))]
    optimize!(proj)

    return objective_value(proj), value.(x)
end

# ==============================================================================
# projecte la solution entiere correspondant au generateur k et test d'admissibilite
# (A::Array{Int,2}, xTilde::Array{Int,1}, 
# c1::Array{Int,1}, c2::Array{Int,1}, k::Int64, protectedIndexOfInt::Vector{Int64}, λ1::Vector{Float64}, λ2::Vector{Float64}, α::Float64)
function projectingSolution!(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur}, k::Int64, 
                             A::Array{Int,2}, c1::Array{Int,1}, c2::Array{Int,1},
                             d::tListDisplay, α::Float64=1.,β::Float64=0.5, MIPproj::Bool=true)

    # --------------------------------------------------------------------------
    # Projete la solution entiere sur le polytope X 
    #generalNadir = fill(tPoint(L[end].y[1],L[1].y[2]),length(vg))
    nadirs = computeLocalNadirs(vg,L)
    λ1,λ2 = calculerDirections2(L,vg)
    #fPrj, vg[k].sPrj.x = Δ2SPA(A,vg[k].sInt.x)
    #fPrj, vg[k].sPrj.x = Δ2SPAbisCone(L,A,vg[k].sInt.x,c1,c2,k,λ1,λ2,nadirs,vg,α,β)

    fPrj, vg[k].sPrj.x = Δ2SPABelgique(A,vg[k].sInt.x,c1,c2,k,λ1,λ2,α) #: Δ2SPAbis(A,vg[k].sInt.x,c1,c2,k,λ1,λ2,α)
    
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
        vg[k].sInt.x = deepcopy(vg[k].sPrj.x)
        vg[k].sInt.y[1] = vg[k].sPrj.y[1]
        vg[k].sInt.y[2] = vg[k].sPrj.y[2]
        vg[k].sFea = true
        @printf("→ Admissible "); print("                       ")

        # archive le point obtenu pour les besoins d'affichage
        if generateurVisualise == -1 
            # archivage pour tous les generateurs
            push!(d.XFeas, vg[k].sPrj.y[1])
            push!(d.YFeas, vg[k].sPrj.y[2])
        elseif generateurVisualise == k
            # archivage seulement pour le generateur k
            push!(d.XFeas, vg[k].sPrj.y[1])
            push!(d.YFeas, vg[k].sPrj.y[2])
        end  

    else

        vg[k].sFea = false
        @printf("→ x          "); print("                       ")
        # prepare pour l'iteration suivante
#        vg[k].xRlx = deepcopy(vg[k].sPrj.x) !!!!!!!!!!!!!
    end
end
