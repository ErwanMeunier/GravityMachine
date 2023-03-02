# ==============================================================================
# The gravity machine (Man of Steel) -> to terraform the world

println("""\nAlgorithme "Gravity machine" --------------------------------\n""")

const verbose = true
const graphic = true

# Display dinamically (step-by-step) the inner operations of gravity machine 
const slowexec = false
const slowtime = 1
# ---

const plotGenerators = true
global generateurVisualise = -1

verbose ? println("-) Active les packages requis\n") : nothing
using JuMP, GLPK, PyPlot, Printf, Random
verbose ? println("  Fait \n") : nothing


graphic ? (println("-) Mise en place de l'affichage\n") ; pygui(true); println(" Fait \n")) : nothing


# ==============================================================================

include("GMdatastructures.jl") # types, datastructures and global variables specially defined for GM
include("GMparsers.jl")        # parsers of instances and non-dominated points
include("GMgenerators.jl")     # compute the generators giving the L bound set
include("GMjumpModels.jl")     # JuMP models for computing relaxed optima of the SPA
include("GMrounding.jl")       # Strategies for rounding a LP-solution to a 01-solution
include("GMdirection.jl")      # Strategies for computing the direction of the projection process
include("GMprojection.jl")     # JuMP models for computing the projection on the polytope of the SPA
include("GMmopPrimitives.jl")  # usuals algorithms in multiobjective optimization
include("GMperturbation.jl")   # routines dealing with the perturbation of a solution when a cycle is detected
include("GMquality.jl")        # quality indicator of the bound set U generated
include("TriangleTools.jl") # Some tools to sample uniformly integral points in simplex

# ==============================================================================
# Ajout d'une solution relachee initiale a un generateur

function ajouterX0!(vg::Vector{tGenerateur}, k::Int64, s::tSolution{Float64})

    vg[k].sRel = deepcopy(s) # met en place le generateur \bar{x}^k
    vg[k].sPrj = deepcopy(s) # le generateur est la premiere projection \bar{x}^{k,0}
    return nothing
end


# ==============================================================================
# Ajout d'une solution entiere (arrondie ou perturbee) a un generateur

function ajouterXtilde!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Int64}, y::Vector{Int64})

    vg[k].sInt.x = copy(x)
    vg[k].sInt.y = copy(y)
    return nothing
end


# ==============================================================================
# Ajout d'une solution fractionnaire (projetee) a un generateur

function ajouterXbar!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Float64}, y::Vector{Float64})

    vg[k].sPrj.x = copy(x)
    vg[k].sPrj.y = copy(y)
    return nothing
end


# ==============================================================================
# Elabore 2 ensembles d'indices selon que xTilde[i] vaut 0 ou 1

function split01(xTilde::Array{Int,1})

   indices0 = Vector{Int64}()
   indices1 = Vector{Int64}()

   sizehint!(indices0,length(xTilde))
   sizehint!(indices1,length(xTilde))

   for i=1:length(xTilde)
       if xTilde[i] == 0
           push!(indices0,i)
       else
           push!(indices1,i)
       end
    end

   return indices0, indices1
end

# ==============================================================================
# Returns two vectors where: 
# The first one stores the index of integer values of x whereas the seconde one does the same for the floating point values of x
function splitByType(x::Vector{Float64})::Tuple{Vector{Int64},Vector{Int64}}
    # Vectors of indexes corresponding to each type: integer or floating point number
    xInt::Vector{Int64} = Vector{Int64}()
    xFloat::Vector{Int64} = Vector{Int64}()
    # Allocates enough memory to store the indexes
    sizehint!(xInt,length(x))
    sizehint!(xFloat,length(x))

    for i in 1:length(x)
        # atol parameter should be set carefully and be consistent with other occurences of isapprox function
        if isapprox(x[i],0,atol=10^-3) || isapprox(x[i],1,atol=10^-3)
            push!(xInt,i)
        else
            push!(xFloat,i)
        end
    end

    return xInt, xFloat
end

# ==============================================================================
# test si une solution est admissible en verifiant si sa relaxation lineaire
# conduit a une solution entiere

function estAdmissible(x::Vector{Float64})

    admissible = true
    i=1
    while admissible && i<=length(x)
        if round(x[i], digits=3)!=0.0 && round(x[i], digits=3)!=1.0
            admissible = false
        end
        i+=1
    end
    return admissible
end


# ==============================================================================
# calcule la performance z d'une solution x sur les 2 objectifs

function evaluerSolution(x::Vector{Float64}, c1::Array{Int,1}, c2::Array{Int,1})

    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return round(z1, digits=2), round(z2, digits=2)
end

# ==============================================================================
# predicat : verifie si une solution entiere est realisable
function isFeasible(vg::Vector{tGenerateur}, k::Int64)
    #verbose && vg[k].sFea == true ? println("   feasible") : nothing
    return (vg[k].sFea == true)
end


# ==============================================================================
# predicat : verifie si le nombre d'essai maximum a ete tente
function isFinished(trial::Int64, maxTrial::Int64)
#    verbose && trial > maxTrial ? println("   maxTrial") : nothing
    return (trial > maxTrial)
end


# ==============================================================================
# predicat : verifie si le budget de calcul maximum a ete consomme
function isTimeout(temps, maxTime)
#    verbose && time()- temps > maxTime ? println("   maxTime") : nothing
    return (time()- temps > maxTime)
end


# ==============================================================================
# elabore pC le pointeur du cone ouvert vers L

function elaborePointConeOuvertversL(vg::Vector{tGenerateur}, k::Int64, pB::tPoint, pA::tPoint)

    # recupere les coordonnees du point projete
    pC=tPoint(vg[k].sPrj.y[1], vg[k].sPrj.y[2])

    # etablit le point nadir pN au depart des points pA et pB adjacents au generateur k
    pN = tPoint( pA.x , pB.y )

#    print("Coordonnees du cone 2 : ")
#    @show pC, pN

    # retient pN si pC domine pN (afin d'ouvrir le cone)
    if (pC.x < pN.x)  &&  (pC.y < pN.y)
        # remplace pC par pN
        pC=tPoint( pA.x , pB.y )
    end

    return pC
end


# ==============================================================================
#= Retourne un booléen indiquant si un point se trouve dans un secteur défini dans
  le sens de rotation trigonométrique (repère X de gauche à droite, Y du haut vers
  le bas).
  https://www.stashofcode.fr/presence-dun-point-dans-un-secteur-angulaire/#more-328
  M    Point dont la position est à tester (point resultant a tester)
  O    Point sommet du secteur (point generateur)
  A    Point de départ du secteur (point adjacent inferieur)
  B    Point d'arrivée du secteur (point adjacent superieur)
  sortie : Booléen indiquant si le point est dans le secteur ou non.

  Exemple :

  B=point(2.0,1.0)
  O=point(2.5,2.5)
  A=point(5.0,5.0)

  M=point(5.0,4.0)
  inSector(M, O, A, B)
=#

function inSector(M, O, A, B)

    cpAB = (A.y - O.y) * (B.x - O.x) - (A.x - O.x) * (B.y - O.y)
    cpAM = (A.y - O.y) * (M.x - O.x) - (A.x - O.x) * (M.y - O.y)
    cpBM = (B.y - O.y) * (M.x - O.x) - (B.x - O.x) * (M.y - O.y)

    if (cpAB > 0)
        if ((cpAM > 0) && (cpBM < 0))
            return true
        else
            return false
        end
    else
        if (!((cpAM < 0) && (cpBM > 0)))
            return true
        else
            return false
        end
    end
end

function inCone(pOrg, pDeb, pFin, pCur)
    # pOrg : point origine du cone (la ou il est pointe)
    # pDeb : point depart du cone (point du rayon [pOrg,pDeb])
    # pFin : point final du cone (point du rayon [pOrg,pFin])
    # pCur : point courant a tester
    # retourne VRAI si pCur est dans le cone pDeb-pFin-pOrg, FAUX sinon

    cp_pDeb_pFin = (pDeb.x - pOrg.x) * (pFin.y - pOrg.y) - (pDeb.y - pOrg.y) * (pFin.x - pOrg.x)
    cp_pDeb_pCur = (pDeb.x - pOrg.x) * (pCur.y - pOrg.y) - (pDeb.y - pOrg.y) * (pCur.x - pOrg.x)
    cp_pFin_pCur = (pFin.x - pOrg.x) * (pCur.y - pOrg.y) - (pFin.y - pOrg.y) * (pCur.x - pOrg.x)

    if (cp_pDeb_pFin > 0)
        if ((cp_pDeb_pCur >= 0) && (cp_pFin_pCur <= 0))
            return true
        else
            return false
        end
    else
        if (!((cp_pDeb_pCur < 0) && (cp_pFin_pCur > 0)))
            return true
        else
            return false
        end
    end
end

function inCone1VersZ(pOrg, pDeb, pFin, pCur)
    return inCone(pOrg, pDeb, pFin, pCur)
end

function inCone2Vers0(pOrg, pDeb, pFin, pCur)
    return !inCone(pOrg, pDeb, pFin, pCur)
end


# ==============================================================================
# Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
function selectionPoints(vg::Vector{tGenerateur}, k::Int64)
    nbgen = size(vg,1)
    if k==1
        # premier generateur (point predecesseur fictif)
        pPrec = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2]+1.0)
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    elseif k==nbgen
        # dernier generateur (point suivant fictif)
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k].sRel.y[1]+1.0, vg[k].sRel.y[2])
    else
        # generateur non extreme
        pPrec = tPoint(vg[k-1].sRel.y[1], vg[k-1].sRel.y[2])
        pCour = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2])
        pSuiv = tPoint(vg[k+1].sRel.y[1], vg[k+1].sRel.y[2])
    end
#    print("Coordonnees du cone 1 : ")
#    @show pPrec, pCour, pSuiv
    return pPrec, pCour, pSuiv
end

# ==============================================================================
# Forces non-integers variables to be integer. Integer variables may become

function transformLowerBoundedSet!(vg::Vector{tGenerateur}, A::Array{Int,2}, λ1::Vector{Float64}, λ2::Vector{Float64}, c1::Vector{Int}, c2::Vector{Int})

    nbvar::Int = size(A,2)
    nbctr::Int = size(A,1)

    for k in eachindex(vg)
        cλ::Vector{Float64} = (1+(λ1[k]-1)).*c1 + (1+(λ2[k]-1)).*c2

        model::Model = Model(GLPK.Optimizer)
        @variable(model, 0<=x[1:nbvar]<=1)
        @objective(model, Min, sum(cλ[j]*x[j] for j in 1:nbvar))
        @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)
        # The new solution must be dominated by the "former generator"
        @expression(model, obj1, sum(c1[j]*x[j] for j in 1:nbvar)) # first objective 
        @expression(model, obj2, sum(c2[j]*x[j] for j in 1:nbvar)) # second objective
        @constraint(model, obj1 >= vg[k].sRel.y[1]) # make a freehand figure and trust your intuition... TODO
        @constraint(model, obj2 >= vg[k].sRel.y[2])
        #---
        [set_binary(model[:x][i]) for i=1:nbvar if !(isapprox(vg[k].sRel.x[i],0,atol=10^-3)||isapprox(vg[k].sRel.x[i],1,atol=10^-3))]
        
        optimize!(model)

        vg[k].sRel.x = value.(x)
        
        arrowBaseX = vg[k].sRel.y[1] # graphic
        arrowBaseY = vg[k].sRel.y[2] # graphic

        vg[k].sRel.y[1], vg[k].sRel.y[2] = evaluerSolution(vg[k].sRel.x,c1,c2)
        
        dX = vg[k].sRel.y[1] - arrowBaseX # graphic
        dY = vg[k].sRel.y[2] - arrowBaseY # graphic  
        graphic ? arrow(arrowBaseX, arrowBaseY, dX, dY, color="fuchsia") : nothing # graphic
    end
end

# ==============================================================================
# point d'entree principal

function GM( fname::String,
             tailleSampling::Int64,
             maxTrial::Int64,
             maxTime::Int64
           )

    @assert tailleSampling>=3 "Erreur : Au moins 3 sont requis"

    @printf("0) instance et parametres \n\n")
    true ? println("  instance = $fname | tailleSampling = $tailleSampling | maxTrial = $maxTrial | maxTime = $maxTime\n\n") : nothing

    # chargement de l'instance numerique ---------------------------------------
    c1, c2, A = loadInstance2SPA(fname) # instance numerique de SPA
    nbctr = size(A,1)
    nbvar = size(A,2)
    nbobj = 2

    # structure pour les points qui apparaitront dans l'affichage graphique [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    d::tListDisplay = tListDisplay([Vector{tSolution}() for k in 1:16]...)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("1) calcule les etendues de valeurs sur les 2 objectifs\n\n")

    # calcule la valeur optimale relachee de f1 seule et le point (z1,z2) correspondant
    f1RL, xf1RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 1) # opt fct 1
    minf1RL, maxf2RL = evaluerSolution(xf1RL, c1, c2)

    # calcule la valeur optimale relachee de f2 seule et le point (z1,z2) correspondant
    f2RL, xf2RL = computeLinearRelax2SPA(nbvar, nbctr, A, c1, c2, typemax(Int), 2) # opt fct 2
    maxf1RL, minf2RL = evaluerSolution(xf2RL, c1, c2)

    verbose ? @printf("  f1_min=%8.2f ↔ f1_max=%8.2f (Δ=%.2f) \n", minf1RL, maxf1RL, maxf1RL-minf1RL) : nothing
    verbose ? @printf("  f2_min=%8.2f ↔ f2_max=%8.2f (Δ=%.2f) \n\n", minf2RL, maxf2RL, maxf2RL-minf2RL) : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("2) calcule les generateurs par e-contrainte alternant minimiser z1 et z2\n\n")

    nbgen, L = calculGenerateurs(A, c1, c2, tailleSampling, minf1RL, maxf2RL, maxf1RL, minf2RL, d)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # allocation de memoire pour la structure de donnees -----------------------

    vg::Vector{tGenerateur} = allocateDatastructure(nbgen, nbvar, nbobj)

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    @printf("3) place L dans structure et verifie l'admissibilite de chaque generateur\n\n")
    
    # TEMPORARY BENCHMARK
    nbcyclestotal = 0
    nbcyclesMax = 0

    for k=1:nbgen

        verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, L[k].y[1], L[k].y[2]) : nothing

        # copie de l'ensemble bornant inferieur dans la stru de donnees iterative ---
        ajouterX0!(vg, k, L[k])
        generateurVisualise = plotGenerators ? k : -1
        # test d'admissibilite et marquage de la solution le cas echeant -------
        if estAdmissible(vg[k].sRel.x)
            ajouterXtilde!(vg, k, convert.(Int, vg[k].sRel.x), convert.(Int, L[k].y))
            vg[k].sFea   = true
            verbose ? @printf("→ Admissible \n") : nothing
            # archive le point obtenu pour les besoins d'affichage    
            if generateurVisualise == -1 
                # archivage pour tous les generateurs
                push!(d.XFeas,vg[k].sInt.y[1])
                push!(d.YFeas,vg[k].sInt.y[2])
            elseif generateurVisualise == k
                # archivage seulement pour le generateur k
                push!(d.XFeas,vg[k].sInt.y[1])
                push!(d.YFeas,vg[k].sInt.y[2])
            end 
        else
            vg[k].sFea   = false
            verbose ? @printf("→ x          \n") : nothing
        end

    end
    verbose ? println("") : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Sortie graphique

    graphic ? figure("Gravity Machine",figsize=(6.5,5)) : nothing
    #xlim(25000,45000)
    #ylim(20000,40000)
    graphic ? xlabel(L"z^1(x)") : nothing
    graphic ? ylabel(L"z^2(x)") : nothing
    graphic ? PyPlot.title("Cone | 1 rounding | 2-$fname") : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # calcule les directions (λ1,λ2) pour chaque generateur a utiliser lors des projections
    λ1,λ2 = calculerDirections2(L,vg)
    # ==========================================================================
    println("3)bis Préparation pour 4) -> tentative d'amélioration des générateurs ")

    transformLowerBoundedSet!(vg,A,λ1,λ2,c1,c2)

    #d.xLf1Improved = [vg[k].sRel.x[1] for k in eachindex(vg)] ;  d.yLf1Improved # liste des points (x,y) relaches améliorés #Recent improvement
    #d.xLf2Improved = [] ;  d.yLf2Improved # liste des points (x,y) relaches améliorés #Recent improvement 
    d.xLImproved = [g.sRel.y[1] for g in vg]; d.yLImproved = [g.sRel.y[2] for g in vg]    # liste des points (x,y) relaches améliorés #Recent improvement

    improvedNadir::tPoint = tPoint(vg[end].sRel.y[1],vg[1].sRel.y[2]) 

    @printf("4) terraformation generateur par generateur \n\n")
    labelInt = 1 # graphical purpose
    #--- Number of trials allowed

    globalNadir = tPoint(L[end].y[1],L[1].y[2])
    #budgetMaxTrials = [Int64(ceil(maxTrial*nbgen/countLP(tPoint(vg[1].sRel.y[1],vg[1].sRel.y[2]),tPoint(vg[3].sRel.y[1],vg[3].sRel.y[2]),globalNadir)))]
    #append!(budgetMaxTrials,[Int64(ceil(maxTrial*nbgen/countLP(tPoint(vg[k-1].sRel.y[1],vg[k-1].sRel.y[2]),tPoint(vg[k+1].sRel.y[1],vg[k+1].sRel.y[2]),globalNadir))) for k in 2:(nbgen-1)])
    #append!(budgetMaxTrials,Int64(ceil(maxTrial*nbgen/countLP(tPoint(vg[nbgen-2].sRel.y[1],vg[nbgen-2].sRel.y[2]),tPoint(vg[nbgen].sRel.y[1],vg[nbgen].sRel.y[2]),globalNadir))))
    #--- End of "Number of trials allowed"
    #println("Budget par générateur : ", budgetMaxTrials)
    nbFeasible = 0
    nbMaxTrials = 0
    nbMaxTime = 0
    for k in [i for i in 1:nbgen if !isFeasible(vg,i)] # ORIGINAL: for k in [i for i in 1:nbgen if !isFeasible(vg,i)]
        temps = time()
        trial = 0
        H = Set{Vector{Int64}}()

        #perturbSolution30!(vg,k,c1,c2,d)

        # rounding solution : met a jour sInt dans vg --------------------------
        #roundingSolution!(vg,k,c1,c2,d)  # un cone
        #roundingSolutionnew24!(vg,k,c1,c2,d) # deux cones
        #protectedIndexOfInt::Vector{Int64}, xFloat::Vector{Int64} = splitByType(vg[k].sRel.x) # TODO: xFloat doesn't matter for now

        arrowBaseX = vg[k].sRel.y[1] # graphic
        arrowBaseY = vg[k].sRel.y[2] # graphic
        roundingSolutionNew23!(vg,k,c1,c2,d) # un cone et LS sur generateur
        labelInt += 1
        dX = vg[k].sInt.y[1] - arrowBaseX
        dY = vg[k].sInt.y[2] - arrowBaseY
        graphic ? plt.arrow(arrowBaseX, arrowBaseY, dX, dY, color="blue") : nothing
        slowexec ? sleep(slowtime) : nothing
        # => Only floating point value are modified so splitByType does have style a sense
        push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])
        verbose ? println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4)) : nothing
        nbcycles = 0
        while !(t1=isFeasible(vg,k)) && !(t2=isFinished(trial, maxTrial)) && !(t3=isTimeout(temps, maxTime))
            trial+=1
            α = 1.# 1/(2^(trial-1))
            β = 0.#0.4 + 0.6*trial/maxTrial
            γ = 1.
            nbcyclesMax = max(nbcyclesMax,nbcycles)
            println("   α = ", α)
            println("   β = ", β)
            # projecting solution : met a jour sPrj, sInt, sFea dans vg --------
            arrowBaseX = vg[k].sInt.y[1] # graphic
            arrowBaseY = vg[k].sInt.y[2] # graphic
            projectingSolution!(L,vg,k,A,c1,c2,d,α,β,trial==1) # first projection uses the integrity constraint 
            labelInt += 1
            dX = vg[k].sPrj.y[1] - arrowBaseX
            dY = vg[k].sPrj.y[2] - arrowBaseY
            graphic ? plt.arrow(arrowBaseX, arrowBaseY, dX, dY, color="red",label=string(labelInt)) : nothing
            slowexec ? sleep(slowtime) : nothing
            println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

            if !isFeasible(vg,k)

                # rounding solution : met a jour sInt dans vg --------------------------
                #roundingSolution!(vg,k,c1,c2,d)
                #roundingSolutionnew24!(vg,k,c1,c2,d)
                arrowBaseX = vg[k].sPrj.y[1] # graphic
                arrowBaseY = vg[k].sPrj.y[2]
                roundingSolutionNew23!(vg,k,c1,c2,d) # graphic
                labelInt+=1
                dX = vg[k].sInt.y[1] - arrowBaseX
                dY = vg[k].sInt.y[2] - arrowBaseY 
                graphic ? plt.arrow(arrowBaseX, arrowBaseY, dX, dY, color="orange",label=string(labelInt)) : nothing
                slowexec ? sleep(slowtime) : nothing
                println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

                # test detection cycle sur solutions entieres ------------------
                cycle = [vg[k].sInt.y[1],vg[k].sInt.y[2]] in H
                if (cycle == true)
                    println("CYCLE!!!!!!!!!!!!!!!")
                    nbcycles += 1
                    # perturb solution
                    arrowBaseX = vg[k].sInt.y[1] # graphic
                    arrowBaseY = vg[k].sInt.y[2]
                    perturbSolution30!(vg,k,c1,c2,d)
                    labelInt+=1
                    dX = vg[k].sInt.y[1] - arrowBaseX
                    dY = vg[k].sInt.y[2] - arrowBaseY
                    graphic ? plt.arrow(arrowBaseX, arrowBaseY, dX, dY, shape="right", color="pink",label=string(labelInt)) : nothing
                    slowexec ? sleep(slowtime) : nothing
                    #perturbSolution40!(vg,k,c1,c2,d,λ1,λ2,γ)
                end
                push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])
            end
        end
        nbcyclestotal += nbcycles

        if t1
            println("   feasible \n")
            nbFeasible+=1
        elseif t2
            println("   maxTrial \n")
            nbMaxTrials+=1
        elseif t3
            println("   maxTime \n")
            nbMaxTime+=1
        end


    end

    println("");

    # ==========================================================================

    @printf("5) Extraction des resultats\n\n")


    for k=1:nbgen
        verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, vg[k].sInt.y[1],vg[k].sInt.y[2]) : nothing
        # test d'admissibilite et marquage de la solution le cas echeant -------
        if vg[k].sFea
            verbose ? @printf("→ Admissible \n") : nothing
        else
            verbose ? @printf("→ x          \n") : nothing
        end
    end

    # allocation de memoire pour les ensembles bornants ------------------------
    U = Vector{tSolution{Int64}}(undef,nbgen)
    for j = 1:nbgen
        U[j] = tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj))
    end
    #--> TODO : stocker l'EBP dans U proprement


    # ==========================================================================
    @printf("6) Edition des resultats \n\n")

    graphic ? figure("Gravity Machine",figsize=(6.5,5)) : nothing
    #xlim(25000,45000)
    #ylim(20000,40000)
    graphic ? xlabel(L"z^1(x)") : nothing
    graphic ? ylabel(L"z^2(x)") : nothing
    # Donne les points relaches initiaux ---------------------------------------
    
    #graphic ? scatter(d.xLf1,d.yLf1,color="green", marker="p") : nothing # generators found considering an ϵ-constraint on f1
    #graphic ? scatter(d.xLf2,d.yLf2,color="pink", marker="p") : nothing # generators found considering an ϵ-constraint on f2
    graphic ? scatter(d.xL,d.yL,color="blue", marker="+", label = L"y \in L") : nothing
    graphic ? scatter(d.xLImproved,d.yLImproved, color="fuchsia", marker=".",label = "Improved generators", lw = 2) : nothing
    # Donne le nadir global  ---------------------------------------------------
    # Des générateurs non améliorés
    println(d.xL)
    println(d.yL)
    println(d.xLImproved)
    println(d.yLImproved)
    graphic ? scatter([L[end].y[1]],[L[1].y[2]],color="blue",marker="*", label = "Non-improved Nadir", lw = 2) : nothing
    # Des générateurs améliorés
    graphic ? scatter([improvedNadir.x],[improvedNadir.y],color="fuchsia",marker="*", label = "Improved Nadir", lw = 2) : nothing
    
    println("Both nadirs points are "*(globalNadir==improvedNadir ? "equal" : "different"))

    # Donne les points entiers -------------------------------------------------
    graphic ? scatter(d.XInt,d.YInt,color="orange", marker="s", label = L"y"*" rounded") : nothing
    graphic ? (@show d.XInt) : nothing
    graphic ? (@show d.YInt) : nothing

    # Donne les points apres projection Δ(x,x̃) ---------------------------------
    graphic ? scatter(d.XProj,d.YProj, color="red", marker="x", label = L"y"*" projected") : nothing
    graphic ? (@show d.XProj) : nothing
    graphic ? (@show d.YProj) : nothing

    # Donne les points admissibles ---------------------------------------------
    graphic ? scatter(d.XFeas,d.YFeas, color="green", marker="o", label = L"y \in F") : nothing
    graphic ? (@show d.XFeas) : nothing
    graphic ? (@show d.YFeas) : nothing

    # Donne l'ensemble bornant primal obtenu + la frontiere correspondante -----
    #--> TODO : stocker l'EBP dans U proprement
    X_EBP_frontiere, Y_EBP_frontiere, X_EBP, Y_EBP = ExtractEBP(d.XFeas, d.YFeas)
    graphic ? plot(X_EBP_frontiere, Y_EBP_frontiere, color="green", markersize=3.0, marker="x") : nothing
    graphic ? scatter(X_EBP, Y_EBP, color="green", s = 150, alpha = 0.3, label = L"y \in U") : nothing
    graphic ? (@show X_EBP) : nothing
    graphic ? (@show Y_EBP) : nothing

    # Donne les points qui ont fait l'objet d'une perturbation -----------------
    graphic ? scatter(d.XPert,d.YPert, color="magenta", marker="s", label ="pertub") : nothing

    # Donne les points non-domines exacts de cette instance --------------------
    XN,YN = loadNDPoints2SPA(fname)
    graphic ? plot(XN, YN, color="black", linewidth=0.75, marker="+", markersize=1.0, linestyle=":", label = L"y \in Y_N") : nothing
    graphic ? scatter(XN, YN, color="black", marker="+") : nothing
    graphic ? (@show XN) : nothing
    graphic ? (@show YN) : nothing

    # Affiche le cadre avec les legendes des differents traces -----------------
    graphic ? legend(bbox_to_anchor=[1,1], loc=0, borderaxespad=0, fontsize = "x-small") : nothing
    #PyPlot.title("Cone | 1 rounding | 2-$fname")

    # Compute the quality indicator of the bound set U generated ---------------
    # Need at least 2 points in EBP to compute the quality indicator
    quality = 0. # TEMPORARY BENCHMARK
    if length(X_EBP) > 1
        quality = qualityMeasure(XN,YN, X_EBP,Y_EBP)
        @printf("Quality measure: %5.2f %%\n", quality*100)
    end
    
    # TEMPORARY TO BENCHMARK
    return quality*100, nbcyclestotal, nbcyclesMax, length(XN), length(X_EBP), nbFeasible, nbMaxTime, nbMaxTrials
end

# ==============================================================================

#@time GM("sppaa02.txt", 6, 20, 20)
#@time GM("sppnw03.txt", 6, 20, 20) #pb glpk
#@time GM("sppnw10.txt", 6, 20, 20)
#@time GM("sppnw16.txt", 6, 20, 20)
#@time GM("sppnw31.txt", 6, 20, 20)
#@time GM("sppnw30.txt", 6, 20, 20)
@time GM("sppnw40.txt", 6, 20, 20)
#@time GM("didactic5.txt", 5, 5, 10)
#@time GM("sppnw29.txt", 6, 30, 20)
#nothing