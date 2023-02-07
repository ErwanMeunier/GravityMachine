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
    idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

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


function perturbSolution40!(vg::Vector{tGenerateur}, k::Int64, c1::Vector{Int64}, c2::Vector{Int64}, d::tListDisplay,λ1::Vector{Float64},λ2::Vector{Float64},γ::Float64)
    nbvar::Int64 = length(vg[k].sInt.x)
    T::Float64 = 1.
    limitPerturb::Int64 = Int64(ceil(T*nbvar))
    # liste des candidats (valeur, indice) et tri decroissant
    
    #idxTilde0, idxTilde1 = split01(vg[k].sInt.x)

    #candidats=[( abs( vg[k].sPrj.x[i] - vg[k].sInt.x[i] ) , i ) for i=1:nbvar if vg[k].sPrj.x[i]>0 && vg[k].sPrj.x[i]<1]
    
    utilities::Vector{Float64} = [c1[i]+c2[i] for i in 1:nbvar]
    maxU::Float64 = maximum(utilities)
    map!(x->x/maxU,utilities,utilities) # normalized utilities
    proba::Vector{Float64} = [(vg[k].sInt.x[i]==1 ? utilities[i] : 1-utilities[i]) for i in 1:nbvar]
    #println("Probabilities : ", proba)
    seq::Vector{Int64} = sortperm(proba,rev=true)
    nbflipped::Int64 = 0
    # independant random draw
    
    fig1, ax1 = plt.subplots()
    plt.plot(1:nbvar,[sum([proba[i]*(vg[k].sInt.x[i]==1 ? c1[i] : -(c1[i])) for j=1:i]) for i=1:nbvar])
    plt.show()
    plt.savefig(string(k)*"c1.png")
    plt.close(fig1)
    fig2, ax2 = plt.subplots()
    plt.plot(1:nbvar,[sum([proba[i]*(vg[k].sInt.x[i]==1 ? c2[i] : -(c2[i])) for j=1:i]) for i=1:nbvar])
    plt.show()
    plt.savefig(string(k)*"c2.png")
    plt.close(fig2)
    
    for i in seq
        if rand() <= (1/2 + γ*( 1/2 + proba[i] )) # γ is very important determine the quantity of "randomness" in the "flipping process"
            nbflipped += 1
            if vg[k].sInt.x[i] == 1
                vg[k].sInt.x[i] = 0
                vg[k].sInt.y[1] = vg[k].sInt.y[1] - c1[i]
                vg[k].sInt.y[2] = vg[k].sInt.y[2] - c2[i]
            else
                vg[k].sInt.x[i] = 1
                vg[k].sInt.y[1] = vg[k].sInt.y[1] + c1[i]
                vg[k].sInt.y[2] = vg[k].sInt.y[2] + c2[i]
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