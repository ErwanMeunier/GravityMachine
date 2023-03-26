using Revise # allowing the redefinition of constant variables while compiling and executing 
# ==============================================================================
# The gravity machine (Man of Steel) -> to terraform the world

println("""\nAlgorithme "Gravity machine" --------------------------------\n""")

const verbose = true

# Figures
const graphic = false
const savegraphic = false # set to false by default
const savingDir = "./results/figuresSol/"

if savegraphic
    try
        mkdir(savingDir)
    catch; Base.IOError # the directory already exists
        println("WARNING: The directory already exists -> Figures inside the directory will still be updated") 
    end
end

# Display dinamically (step-by-step) the inner operations of gravity machine 
const slowexec = false
const slowtime = 0.
# ---

const plotGenerators = true
global generateurVisualise = -1

global CHOICE_ROUNDING = 2 # FROM 1 TO 3
global CHOICE_PROJECTION = 5 # FROM 1 TO 5
global CHOICE_COMPUTEDIRECTIONS = 2 # FROM 1 TO 4
global CHOICE_PERTUBATION = 2 # FROM 1 TO 3
global CONES_CONSTRAINED_IMPROVE_GENERATORS = false

global maxRatioBinaryVariables::Float64 = 1.# Must be between 0 and 1 meaning 0% to 100%

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

function split01(xTilde::Vector{Int})

   indices0::Vector{Int64} = Vector{Int64}()
   indices1::Vector{Int64} = Vector{Int64}()

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

function estAdmissible(x::Vector{Float64})::Bool

    admissible::Bool = true
    i::Int=1
    while admissible && i<=length(x)
        if !isapprox(x[i], 0.0; atol=10^-3) && !isapprox(x[i],1.0;atol=10^-3)
            admissible = false
        end
        i+=1
    end
    return admissible
end

# ==============================================================================
# Tests whether the given solution is feasible and not previously already starting from a previous generator, or not
# Leads to a integer solution

function estAdmissibleDiff(x::Vector{Float64}, d::tListDisplay)::Bool
    #TODO
end

# ==============================================================================
# calcule la performance z d'une solution x sur les 2 objectifs

function evaluerSolution(x::Vector{Float64}, c1::Array{Int,1}, c2::Array{Int,1})

    z1 = 0.0; z2 = 0.0
    for i in 1:length(x)
        z1 += x[i] * c1[i]
        z2 += x[i] * c2[i]
    end
    return round(z1, digits=3), round(z2, digits=3)
end

# ==============================================================================
# Nettoyage des valeurs des variables d'une solution x relachee sur [0,1]

function nettoyageSolution!(x::Vector{Float64})
    # TODO : using isapprox function could be better
    nbvar::Int = length(x)
    for i in 1:nbvar
        if     isapprox(x[i],0.,atol=10^-3)
                   x[i] = 0.0
        elseif isapprox(x[i],1.,atol=10^-3)
                   x[i] = 1.0
        else
                   x[i] = round(x[i], digits=3)
        end
    end
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
# ===============================================================================
# from : https://discourse.julialang.org/t/help-writing-a-timeout-macro/16591/7
macro timeout(seconds, expr, fail)
    quote
        tsk = @task $esc(expr)
        schedule(tsk)
        Timer($(esc(seconds))) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $(esc(fail))
        end
    end
end

# ==============================MACROS FOR PLOTTING ARROWS============================

macro makearrow(expr, xbefore, ybefore, xafter, yafter, color)
    quote
        if $(esc(graphic)) 
            $(esc(slowexec)) ? sleep($(esc(slowtime))) : nothing # slowing the display
            arrowBaseX::Float64 = $(esc(xbefore))
            arrowBaseY::Float64 = $(esc(ybefore))
            $(esc(expr))
            dX::Float64 = $(esc(xafter)) - arrowBaseX
            dY::Float64 = $(esc(yafter)) - arrowBaseY
            println(dX,dY)
            arrow(arrowBaseX, arrowBaseY, dX, dY, color=$(color))
        else 
            $(esc(expr))
        end
    end
end

# ==============================================================================

# ==============================================================================
# Forces non-integers variables to be integer. Integer variables may become

function transformLowerBoundedSet!(vg::Vector{tGenerateur}, A::Array{Int,2}, L::Vector{tSolution{Float64}}, λ1::Vector{Float64}, λ2::Vector{Float64}, c1::Vector{Int}, c2::Vector{Int}, d::tListDisplay)::Vector{tSolution{Float64}}
    nbvar::Int = size(A,2)
    nbctr::Int = size(A,1)
    println("NBVAR :", nbvar)
    println("TRANSFORMING LOWER BOUNDED SET")
    nadir = tPoint(L[end].y[1],L[1].y[2])
    nbgen = length(L)
    for k in eachindex(vg)
        #cλ::Vector{Float64} = λ1[k].*c2 + λ2[k].*c1 
        println("k=",k)
        cλ::Vector{Float64} = λ1[k].*c1 + λ2[k].*c2 # 

        model::Model = Model(GLPK.Optimizer)
        @variable(model, 0<=x[1:nbvar]<=1)
        
        @objective(model, Min, sum(cλ[j]*x[j] for j in 1:nbvar))
        @constraint(model, [i=1:nbctr],(sum((x[j]*A[i,j]) for j in 1:nbvar)) == 1)

        @expression(model, z1, sum(c1[j]*x[j] for j in 1:nbvar))
        @expression(model, z2, sum(c2[j]*x[j] for j in 1:nbvar))
        
        if CONES_CONSTRAINED_IMPROVE_GENERATORS && (1<k) && (k<nbgen)
            u = tPoint(L[k-1].y[1],L[k-1].y[2])
            v = tPoint(L[k+1].y[1],L[k+1].y[2])
            if k==2
                α2 = (nadir.y - v.y)/(nadir.x - v.x)
                β2 = nadir.y - α2*nadir.x
                @constraint(model, z2 <= nadir.y) # TODO: can be lower bounded
                @constraint(model, α2*z1+β2 <= z2)
            elseif k==(nbgen-1)
                α1 = (nadir.y - u.y)/(nadir.x - u.x)
                β1 = nadir.y - α1*nadir.x
                @constraint(model, α1*z1+β1 >= z2)
                @constraint(model, z1 <= nadir.x) # TODO: can be upper bounded
            else 
                α1 = (nadir.y - u.y)/(nadir.x - u.x)
                β1 = nadir.y - α1*nadir.x
                α2 = (nadir.y - v.y)/(nadir.x - v.x)
                β2 = nadir.y - α2*nadir.x
                @constraint(model, α1*z1+β1 >= z2)
                @constraint(model, α2*z1+β2 <= z2)
            end
        end
        idx = [l for l in 1:(max(1,k-1)) if vg[l].sFea]
        #println("INDEX : ", idx)

        for q in idx
            N0, N1 = split01(vg[q].sInt.x)
            println("Generator : ", q)
            println("y1=", vg[q].sInt.y[1])
            println("y2=", vg[q].sInt.y[2])
            #=
            @expression(model, obj1[q], sum(c1[j]*x[j] for j in 1:nbvar))
            @expression(model, obj2[q], sum(c2[j]*x[j] for j in 1:nbvar))

            @constraint(model, vg[q].sInt.y[1] <= obj1)
            @constraint(model, vg[q].sInt.y[2] <= obj2)
            =#
            @constraint(model, sum(x[j] for j in N0) + sum(1-x[j] for j in N1) >= 1)
            #rightMember = 0.001*sum(cλ)
            #sum(cλ[j]*x[j] for j in N0) + sum(cλ[j] for j in N1)
            #@constraint(model, sum(cλ[j]*x[j] for j in N0) + sum(cλ[j]*(1-x[j]) for j in N1) >= 1)
        end
        
        #maxRatioBinaryVariables
        
        idxNonInt::Vector{Int} = [i for i=1:nbvar if !(isapprox(vg[k].sRel.x[i],0.,atol=10^-3)||isapprox(vg[k].sRel.x[i],1.,atol=10^-3))]
        println("Number of non integral variables:", length(idxNonInt))
        idxNonInt = sort(idxNonInt, by=i->abs(vg[k].sRel.x[i]-1/2)) # The more x is close to 1/2 the more it is lucky to become a binary variable

        nbBinVar = Int(ceil((maxRatioBinaryVariables * length(idxNonInt))))
        println("Number of variables set integral: ", nbBinVar)
        for i in idxNonInt[1:min(end,nbBinVar)]
            set_binary(model[:x][i])
        end
        
        optimize!(model)

        vg[k].sRel.x = value.(x)

        @makearrow begin vg[k].sRel.y[1], vg[k].sRel.y[2] = evaluerSolution(vg[k].sRel.x,c1,c2) end vg[k].sRel.y[1] vg[k].sRel.y[2] vg[k].sRel.y[1] vg[k].sRel.y[2] "fuchsia"

        if estAdmissible(vg[k].sRel.x);
            println("Admissibilité du générateur amélioré ", k)
            vg[k].sFea = true
            ajouterXtilde!(vg, k, convert.(Int, round.(vg[k].sRel.x)), convert.(Int, round.(vg[k].sRel.y)))
            if generateurVisualise == -1 
                # archivage pour tous les generateurs
                push!(d.XFeas,round(vg[k].sRel.y[1]))
                push!(d.YFeas,round(vg[k].sRel.y[2]))
            elseif generateurVisualise == k
                # archivage seulement pour le generateur k
                push!(d.XFeas,round(vg[k].sRel.y[1]))
                push!(d.YFeas,round(vg[k].sRel.y[2]))
            end 
        else
            println("Le générateur amélioré ", k, " n'est pas admissible")
        end
    end
    
    #=println("Number of diff values: ", 
        length([v for v in redundantVar([Float64.(vg[l].sInt.x) for l in eachindex(vg) if vg[l].sFea])[1] if v==0])
        )=#
    println("FIN DE L AMELIORATION")
    Limproved = [tSolution(deepcopy(vg[k].sRel.x),deepcopy(vg[k].sRel.y)) for k in eachindex(vg)]
    #Lrounded = [tSolution(round.(g.x)) for g in Limproved]
    for k in eachindex(vg)
        println("[",Limproved[k].y[1],";",Limproved[k].y[2],"]")
        g = Limproved[k].x
        restrictedLimproved = [convert.(Bool, round.(Limproved[i].x)) for i in eachindex(vg) if i!=k]
        #println("Coordonnées à 1: ", findall(convert.(Bool, round.(g))))
        println("Générateur courant a déjà une solution soeur: ", convert.(Bool, round.(g)) in restrictedLimproved)
    end
    #break # TO BE REMOVED
    return Limproved #  (new) improved Lower Bound Set
end

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

    # structure pour les points qui apparaitront dans l'affichage graphique
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
    
    # ==========================================================================
    # ANALYSIS OF THE FIRST (NON-IMPROVED) GENERATORS
    for k=1:nbgen

        verbose ? @printf("  %2d  : [ %8.2f , %8.2f ] ", k, L[k].y[1], L[k].y[2]) : nothing

        # copie de l'ensemble bornant inferieur dans la stru de donnees iterative ---
        ajouterX0!(vg, k, L[k])

        generateurVisualise = plotGenerators ? k : -1

        #function ajouterXtilde!(vg::Vector{tGenerateur}, k::Int64, x::Vector{Int64}, y::Vector{Int64})

        # test d'admissibilite et marquage de la solution le cas echeant -------
        if estAdmissible(vg[k].sRel.x)
            ajouterXtilde!(vg, k, convert.(Int, vg[k].sRel.x), convert.(Int, L[k].y))
            vg[k].sFea   = true
            verbose ? @printf("→ Admissible \n") : nothing
            # archive le point obtenu pour les besoins d'affichage    
            if generateurVisualise == -1 
                # archivage pour tous les generateurs
                push!(d.XFeas,vg[k].sInt.y[1]) # instead of sInt
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

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # Sortie graphique

    graphic ? figure("Gravity Machine",figsize=(6.5,5)) : nothing
    #xlim(25000,45000)
    #ylim(20000,40000)
    graphic ? xlabel(L"z^1(x)") : nothing
    graphic ? ylabel(L"z^2(x)") : nothing
    graphic ? title("Cone | 1 rounding | 2-$fname") : nothing

    # --------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # calcule les directions (λ1,λ2) pour chaque generateur a utiliser lors des projections
    

    #d.xLf1Improved = [vg[k].sRel.x[1] for k in eachindex(vg)] ;  d.yLf1Improved # liste des points (x,y) relaches améliorés #Recent improvement
    #d.xLf2Improved = [] ;  d.yLf2Improved # liste des points (x,y) relaches améliorés #Recent improvement 
    
    println("3)bis Préparation pour 4) -> tentative d'amélioration des générateurs ")
    # ANALYSIS OF THE SECOND (IMPROVED) GENERATORS
    
    λ1::Vector{Float64}, λ2::Vector{Float64} = interface_computeDirections(L,vg)
    verbose ? println("---> Directions computed") : nothing
    
    Limproved::Vector{tSolution{Float64}} = transformLowerBoundedSet!(vg,A,L,λ1,λ2,c1,c2,d)
    verbose ? println("---> Generators improved") : nothing

    λ1, λ2 = interface_computeDirections(L,vg) # TODO: with new generators
    verbose ? println("---> Directions updated according to the improved generators") : nothing

    #d.xLImproved = [g.sRel.y[1] for g in vg]; d.yLImproved = [g.sRel.y[2] for g in vg]    # liste des points (x,y) relaches améliorés #Recent improvement
    d.xLImproved = [g.sRel.y[1] for g in vg]; d.yLImproved = [g.sRel.y[2] for g in vg]    # liste des points (x,y) relaches améliorés #Recent improvement

    improvedNadir::tPoint = tPoint(Limproved[end].y[1],Limproved[1].y[2]) 

    @printf("4) terraformation generateur par generateur \n\n")
    labelInt = 1 # graphical purpose
    #--- Number of trials allowed

    globalNadir = tPoint(L[end].y[1],L[1].y[2])
    
    #println("Budget par générateur : ", budgetMaxTrials)
    nbFeasible = 0
    nbMaxTrials = 0
    nbMaxTime = 0

    idxFeasibleGenerators::Vector{Int} = [k for k in 1:nbgen if isFeasible(vg,k)]
    solutionsHist::Set{Vector{Int}} = Set{Vector{Int}}([vg[k].sInt.x for k in idxFeasibleGenerators]) # tabu memory
    
    antecedantPoint::Dict{Vector{Int},Vector{Vector{Float64}}} = Dict{Vector{Int},Vector{Vector{Float64}}}([Pair(copy(vg[k].sInt.x),[]) for k in idxFeasibleGenerators]...)
    #=
    Contain : 
    Feasible solution and their antecedants which may be:
    -> Some improved generators
    -> Some perturbed solutions
    =#

    # archiving of improved generators which lead to feasible solution
    for k in idxFeasibleGenerators
        push!(antecedantPoint[vg[k].sInt.x],Limproved[k].x)
    end

    for k in [i for i in 1:nbgen if !isFeasible(vg,i)]
        temps = time()
        trial = 0
    
        # rounding solution : met a jour sInt dans vg --------------------------
        #roundingSolution!(vg,k,c1,c2,d)  # un cone
        #roundingSolutionnew24!(vg,k,c1,c2,d) # deux cones

        @makearrow(interface_roundingSolution!(vg,k,c1,c2,d),vg[k].sRel.y[1],vg[k].sRel.y[2],vg[k].sInt.y[1],vg[k].sInt.y[2],"orange")

        slowexec ? sleep(slowtime) : nothing
        # => Only floating point value are modified so splitByType does have style a sense
        verbose ? println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4)) : nothing
        nbcycles = 0
        nbcyclesSameSol = 0
        
        while !(t1=isFeasible(vg,k)) && !(t2=isFinished(trial, maxTrial)) && !(t3=isTimeout(temps, maxTime))
            H::Vector{Vector{Int}} = Vector{Vector{Int}}()
            trial+=1
            α = 1.#1/(2^(trial-1)) # TODO
            β = 0.#0.4 + 0.6*trial/maxTrial
            γ = 1.
            nbcyclesMax = max(nbcyclesMax,nbcycles)
            println("   α = ", α)
            println("   β = ", β)
            # projecting solution : met a jour sPrj, sInt, sFea dans vg --------
            temphist = vg[k].sInt 
            if nbcyclesSameSol > 0 
                @makearrow projectingSolution!(A,vg,k,c1,c2,d,λ1,λ2,α,solutionsHist,vg[k].sPrj.y[1],nbcyclesSameSol) vg[k].sPrj.y[1] vg[k].sPrj.y[2] vg[k].sPrj.y[1] vg[k].sPrj.y[2] "red"    
            else
                @makearrow projectingSolution!(A,vg,k,c1,c2,d,λ1,λ2,α,solutionsHist,fill(0,length(vg[k].sRel.x)),nbcyclesSameSol) vg[k].sInt.y[1] vg[k].sInt.y[2] vg[k].sPrj.y[1] vg[k].sPrj.y[2] "red"
            end
            #slowexec ? sleep(slowtime) : nothing
            println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))
            #if isFeasible(vg,k) && 
            if isFeasible(vg,k)
                if Int.(vg[k].sPrj.x) in solutionsHist
                    push!(antecedantPoint[Int.(vg[k].sPrj.x)],copy(vg[k].sInt.x))
                    vg[k].sFea = false
                    #@makearrow perturbSolutionInt!(vg,k,c1,c2,d,antecedentPoint[vg[k].sInt.x]) vg[k].sPrj.y[1] vg[k].sPrj.y[2] vg[k].sInt.y[1] vg[k].sInt.y[2] "pink"
                    @makearrow perturbSolution45!(vg,k,c1,c2,d,λ1,λ2,nbcyclesSameSol,antecedantPoint[Int.(vg[k].sPrj.x)],1.) vg[k].sPrj.y[1] vg[k].sPrj.y[2] vg[k].sPrj.y[1] vg[k].sPrj.y[2] "cyan"
                    nbcyclesSameSol += 1
                else # the new solution is added to the tabulist
                    #antecedantPoint[vg[k].sPrj.x] = [Int.(vg[k].sInt.x)]
                    push!(solutionsHist,Int.(vg[k].sPrj.x))
                    antecedantPoint[Int.(vg[k].sPrj.x)] = [Int.(vg[k].sInt.x)] # initializer
                end
            else
                # rounding solution : met a jour sInt dans vg --------------------------

                @makearrow interface_roundingSolution!(vg,k,c1,c2,d) vg[k].sPrj.y[1] vg[k].sPrj.y[2] vg[k].sInt.y[1] vg[k].sInt.y[2] "orange"
                slowexec ? sleep(slowtime) : nothing

                push!(H,[vg[k].sInt.y[1],vg[k].sInt.y[2]])

                println("   t=",trial,"  |  Tps=", round(time()- temps, digits=4))

                # test detection cycle sur solutions entieres ------------------
                cycle = (vg[k].sInt.y[1],vg[k].sInt.y[2]) in H
                if (cycle == true)
                    println("CYCLE!!!!!!!!!!!!!!!")
                    nbcycles += 1
                    # perturb solution
                    @makearrow perturbSolution30!(vg,k,c1,c2,d) vg[k].sInt.y[1] vg[k].sInt.y[2] vg[k].sInt.y[1] vg[k].sInt.y[2] "pink"
                    slowexec ? sleep(slowtime) : nothing
                    #perturbSolution40!(vg,k,c1,c2,d,λ1,λ2,γ)
                end
            end
        end
        nbcyclestotal += nbcycles

        if t1
            println("   feasible \n")
            push!(solutionsHist,vg[k].sInt.x)
            nbFeasible+=1
        elseif t2
            println("   maxTrial \n")
            nbMaxTrials+=1
        elseif t3
            println("   maxTime \n")
            nbMaxTime+=1
        end


    end
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
    
    if savegraphic
        savefig(savingDir*fname*".png")
        close()
    end
    # TEMPORARY TO BENCHMARK
    return quality*100, nbcyclestotal, nbcyclesMax, length(XN), length(X_EBP), nbFeasible, nbMaxTime, nbMaxTrials
end

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

# ==============================================================================
#=for name in readdir("../SPA/instances/")
    @timeout 30 GM(name[4:end], 6, 20, 20) nothing
end=#
#[GM(name[4:end], 6, 20, 20) for name in readdir("../SPA/instances")[20:end]]
#@time GM("sppaa02.txt", 6, 20, 20)
#@time GM("sppnw03.txt", 6, 20, 20) #pb glpk
#@time GM("sppnw01.txt", 6, 20, 20)
#@time GM("sppnw06.txt", 6, 20, 20)
@time GM("sppnw31.txt", 6, 20, 20)
#@time GM("sppnw30.txt", 6, 20, 20)
#@time GM("sppnw40.txt", 6, 20, 20)
#@time GM("didactic5.txt", 5, 5, 10)
#@time GM("sppnw29.txt", 6, 30, 20)
#nothing