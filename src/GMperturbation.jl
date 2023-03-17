# ==============================================================================
# applique une perturbation sur la solution entiere faisant l'objet d'un cycle

function perturbSolution!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)
    # nombre de variables maximum a considerer
    T = 100
    # nombre effectif de variables a flipper
    TT = rand(T/2:3*T/2)

    # liste des candidats (valeur, indice) et tri decroissant
    nbvar = length(vg[k].sInt.x)
    candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar]
    println(length([true for i=1:nbvar if vg[k].sPrj.x[i]==0.0]))
    println("Candidats : ", length(candidats))
    sort!(candidats, rev=true, by = x -> x[1])
#    sort!(candidats,  by = x -> x[1])

#@show vg[k].sPrj.x
#@show vg[k].sInt.x
#@show candidats
#@show TT
#@show nbvar

    i = 1
    while (i<= nbvar) && (i<=TT)
        j=candidats[i][2]
#        @show candidats[i][2]
        if vg[k].sInt.x[j] == 0
            vg[k].sInt.x[j] = 1
            vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] + c2[j]
        else
            vg[k].sInt.x[j] = 0
            vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[j]
            vg[k].sInt.y[2] = vg[k].sInt.y[2] - c2[j]
        end
        i+=1
    end
    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage    
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    end 
#    @show vg[k].sInt.x
    return nothing
end


# ==============================================================================
# applique une perturbation sur la solution entiere faisant l'objet d'un cycle

function perturbSolution30!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)
    T::Int64 = 10
    # liste des candidats (valeur, indice) et tri decroissant
    nbvar = length(vg[k].sInt.x)
    #idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

    #candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]<1]

    candidats=[( vg[k].sPrj.x[i] , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0. && vg[k].sPrj.x[i]<1.]
#    sort!(candidats, rev=true, by = x -> x[1])
    println("Candidats : ", candidats)
    println("Nombre de candidats : ", length(candidats))

    #@show vg[k].sPrj.x
    #@show vg[k].sInt.x
    #@show candidats
    #@show nbvar

    seq = randperm(length(candidats)) # melange les candidats afin d'avoir une composante variee
    println("Seq : ", seq)
    println("Etat : ", )
    etat = vg[k].sInt.x[ candidats[seq[1]][2] ] # etat 0 ou 1 de la premiere variable candidate
    for i = 1:length(candidats)
        j=candidats[seq[i]][2]
        #@show candidats[seq[i]][2]
        if etat == 0
            if vg[k].sInt.x[j] == 0
                vg[k].sInt.x[j] = 1
                vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[j]
                vg[k].sInt.y[2] = vg[k].sInt.y[2] + c2[j]
            end
        else
            if vg[k].sInt.x[j] == 1
                vg[k].sInt.x[j] = 0
                vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[j]
                vg[k].sInt.y[2] = vg[k].sInt.y[2] - c2[j]
            end
        end
        etat=(etat+1)%2
    end

    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage    
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    end     
#    @show vg[k].sInt.x

    return nothing
end


function perturbSolution40!(vg::Vector{tGenerateur}, k::Int64, c1::Vector{Int64}, c2::Vector{Int64}, d::tListDisplay,λ1::Vector{Float64},λ2::Vector{Float64},nbcyclesSameSol::Int,γ::Float64=0.5)
    nbvar::Int64 = length(vg[k].sInt.x)
    T::Float64 = nbcyclesSameSol
    limitPerturb::Int64 = Int64(ceil(T*nbvar/sum(c1)))
    # liste des candidats (valeur, indice) et tri decroissant
    
    #idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

    #candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]

    utilities::Vector{Float64} = [λ1[k]*c1[i]+λ2[k]*c2[i] for i in 1:nbvar]
    maxU::Float64 = maximum(utilities)
    map!(x->x/(maxU+1),utilities,utilities) # normalized utilities
    proba::Vector{Float64} = [(vg[k].sInt.x[i]==1 ? utilities[i] : 1-utilities[i]) for i in 1:nbvar]
    #println("Probabilities : ", proba)
    #seq::Vector{Int64} = sortperm(proba,rev=true)
    seq = sortperm(1:nbvar,rev=true)
    nbflipped::Int64 = 0
    
    for i in seq
        if rand() <= (1/2 + γ*( proba[i] -1/2)) # γ is very important by fixing the quantity of "randomness" in the "flipping process"
            nbflipped += 1
            if vg[k].sPrj.x[i] == 1
                vg[k].sPrj.x[i] = 0
                vg[k].sPrj.y[1] -= c1[i]
                vg[k].sPrj.y[2] -= c2[i]
            else
                vg[k].sPrj.x[i] = 1
                vg[k].sPrj.y[1] += c1[i]
                vg[k].sPrj.y[2] += c2[i]
            end            
        end
        nbflipped > limitPerturb ? (println("BREAK"); break) : nothing 
    end

    #sort!(candidats, rev=true, by = x -> x[1])
    #@show vg[k].sPrj.x
    #@show vg[k].sInt.x
    #@show candidats
    #@show nbvar

    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage    
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XPert,vg[k].sPrj.y[1])
        push!(d.YPert,vg[k].sPrj.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XPert,vg[k].sPrj.y[1])
        push!(d.YPert,vg[k].sPrj.y[2])
    end     
#    @show vg[k].sInt.x

    return nothing
end

# TODO the returned vector will be very sparse...
# Returns the number of 
function redundantVar(solutions::Vector{Vector{Float64}})::Tuple{Vector{Int64},Vector{Int64}}
    n::Int = length(solutions[1]) # number of variables
    m::Int = length(solutions) # number of solutions
    allequal::Vector{Bool} = trues(n)
    nbEqual::Vector{Int} = zeros(n)
    for i in 1:n
        k::Int = 2
        while allequal[i] && k<=m
            allequal[i] &= (solutions[1][i] == solutions[k][i]) && (isapprox(solutions[k][i],0;atol=10^-3) || isapprox(solutions[k][i],1;atol=10^-3))
            nbEqual[i] += 1
            k+=1
        end
        !allequal[i] ? nbEqual[i]=0 : nothing
    end
    return nbEqual, convert.(Int,round.(solutions[1]))
end

function redundantVar2(solutions::Vector{Vector{Float64}})::Tuple{Vector{Int64},Vector{Int64}}
    n::Int = length(solutions[1]) # number of variables
    m::Int = length(solutions) # number of solutions
    allequal::Vector{Bool} = trues(n)
    nbEqual::Vector{Int} = zeros(n)
    for i in 1:n
        k::Int = 2
        while allequal[i] && k<=m
            allequal[i] &= (solutions[1][i] == solutions[k][i]) && (isapprox(solutions[k][i],0;atol=10^-3) || isapprox(solutions[k][i],1;atol=10^-3))
            nbEqual[i] += 1
            k+=1
        end
        !allequal[i] ? nbEqual[i]=0 : nothing
    end
    return nbEqual, convert.(Int,round.(solutions[1]))
end


function perturbSolution45!(vg::Vector{tGenerateur}, k::Int64, c1::Vector{Int64}, c2::Vector{Int64}, d::tListDisplay,λ1::Vector{Float64},λ2::Vector{Float64},nbcyclesSameSol::Int,antecedents::Vector{Vector{Float64}},γ::Float64=0.5)
    nbvar::Int64 = length(vg[k].sInt.x)
    T::Float64 = nbcyclesSameSol/2
    limitPerturb::Int64 = Int64(ceil(T*nbvar/sum(c1)))
    # liste des candidats (valeur, indice) et tri decroissant
    
    #idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

    #candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]
    redundanceVar::Vector{Int}, correspondingValues::Vector{Int} = redundantVar(antecedents)
    
    utilities::Vector{Float64} = [λ1[k]*c1[i]+λ2[k]*c2[i] for i in 1:nbvar]
    maxU::Float64 = maximum(utilities)
    map!(x->x/(maxU+1),utilities,utilities) # normalized utilities
    proba::Vector{Float64} = [(vg[k].sInt.x[i]==1 ? utilities[i] : 1-utilities[i]) for i in 1:nbvar]
    #println("Probabilities : ", proba)
    #seq::Vector{Int64} = sortperm(proba,rev=true)
    #seq = sortperm(1:nbvar,rev=true)
    nbflipped::Int64 = 0
    #println(redundanceVar)
    seq = sortperm(redundanceVar,rev=true)
    
    for i in seq
        if rand() <= (1/2 + γ*( proba[i] -1/2)) # γ is very important by fixing the quantity of "randomness" in the "flipping process"
            nbflipped += 1
            if vg[k].sPrj.x[i] == 1
                vg[k].sPrj.x[i] = 0
                vg[k].sPrj.y[1] -= c1[i]
                vg[k].sPrj.y[2] -= c2[i]
            else
                vg[k].sPrj.x[i] = 1
                vg[k].sPrj.y[1] += c1[i]
                vg[k].sPrj.y[2] += c2[i]
            end            
        end
        nbflipped > limitPerturb ? (println("BREAK"); break) : nothing 
    end

    #sort!(candidats, rev=true, by = x -> x[1])
    #@show vg[k].sPrj.x
    #@show vg[k].sInt.x
    #@show candidats
    #@show nbvar

    @printf("  %2dC : [ %8.2f , %8.2f ] \n", k, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage    
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XPert,vg[k].sPrj.y[1])
        push!(d.YPert,vg[k].sPrj.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XPert,vg[k].sPrj.y[1])
        push!(d.YPert,vg[k].sPrj.y[2])
    end     
#    @show vg[k].sInt.x

    return nothing
end

function perturbSolutionInt!(vg::Vector{tGenerateur}, k::Int64, c1::Vector{Int64}, c2::Vector{Int64}, d::tListDisplay, antecedentPoint::Vector{Vector{Int}})
    redundantVar::Vector{Int}, associatedValues::Vector{Int} = redundantVar(antecedentPoint)
    maxRV = maximum(redundantVar)
    normalizedRedudantVar::Vector{Float64} = collect(map(x->(x+1)/(maxRV+1),redundantVar))

    for i in eachindex(vg[k].sInt.x)
        if rand() <= normalizedRedudantVar[i]
            vg[k].sInt.x[i] = 1-vg[k].sInt.x[i] 
            if vg[k].sInt.x[i]==1 # => previously set to zero
                vg[k].sInt.y[1] += c1[i]
                vg[k].sInt.y[2] += c2[i]
            else # => previously set to one
                vg[k].sInt.y[1] -= c1[i]
                vg[k].sInt.y[2] -= c2[i]
            end 
        end
    end

    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XPert,vg[k].sInt.y[1])
        push!(d.YPert,vg[k].sInt.y[2])
    end     
end