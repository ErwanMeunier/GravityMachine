# ===================Tools for rounding methods===============================
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

# =============================Rounding methods====================================

# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec cone inferieur seulement
function roundingSolution!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
        else
            vg[k].sInt.x[i] = -1
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

    print("Depart a arrondir : ")
    @show pPrec, pCour, pSuiv, pM
    if length(vg[k].sPrj.x) ≤ 20 (@show vg[k].sPrj.x) end
    if length(vg[k].sInt.x) ≤ 20 (@show vg[k].sInt.x) end

    nbVarNonEntiere = 0
    for i in 1:nbvar
        if vg[k].sInt.x[i] == -1

            # la variable i est non entiere
            nbVarNonEntiere += 1

            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[i]=0
            pM = tPoint( z1 - vg[k].sPrj.x[i] * c1[i] , z2 - vg[k].sPrj.x[i] * c2[i])

            if inSector(pM, pCour, pSuiv, pPrec)
                # le point pM obtenu est hors du cone => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
                @printf("  x[%2d]=1 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            else
                # le point pM obtenu est dans le cone => valide x[i]=0
                vg[k].sInt.x[i] = 0
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            end
        end

        # calcule la performance de la solution entiere
        vg[k].sInt.y[1] += vg[k].sInt.x[i] * c1[i]
        vg[k].sInt.y[2] += vg[k].sInt.x[i] * c2[i]
    end

    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    end  

end


# ==============================================================================
# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec cone inferieur et superieur
function roundingSolutionnew24!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
        else
            vg[k].sInt.x[i] = -1
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

#    print("Depart a arrondir : ")
#    @show pPrec, pCour, pSuiv, pM
#    if length(vg[k].sPrj.x) ≤ 20 @show vg[k].sPrj.x end
#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end

    verbose ? @printf("  %2dR : [ %8.2f , %8.2f ] ", k, z1, z2) : nothing

    nbVarNonEntiere = 0
    for i in 1:nbvar
        if vg[k].sInt.x[i] == -1

            # la variable i est non entiere
            nbVarNonEntiere += 1

            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[i]=0
            pM0 = tPoint( z1 - vg[k].sPrj.x[i] * c1[i] , z2 - vg[k].sPrj.x[i] * c2[i])
            pM1 = tPoint( z1 - vg[k].sPrj.x[i] * c1[i]+c1[i] , z2 - vg[k].sPrj.x[i] * c2[i]+c2[i])

            if !inCone1VersZ(pCour, pSuiv, pPrec, pM0)
                # le point pM0 obtenu est hors du cone inferieur => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
#                @printf("  x[%2d]=1 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            elseif inCone2Vers0(pC, pSuiv, pPrec, pM0)
                # le point pM0 obtenu est dans le cone => valide x[i]=0
                vg[k].sInt.x[i] = 0
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
#                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            elseif inCone2Vers0(pC, pSuiv, pPrec, pM1)
                # le point pM1 obtenu est dans le cone => valide x[i]=1
                vg[k].sInt.x[i] = 1
                z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
#                @printf("  x[%2d]=0 |  : [ %12.5f , %12.5f ] \n",i, z1, z2)
            else
#                println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                # choisir entre pM0 et pM1 selon la distance a la droite [O;OO]
                O = tPoint(vg[k].sRel.y[1], vg[k].sRel.y[2]) # generateur
                OO = tPoint(vg[k].sPrj.y[1], vg[k].sPrj.y[2]) # projection
                # les points pM0 et pM1 obtenus sont hors du cone => choix de 0 ou 1 avec la plus courte distance à O et OO
                distM0 = sqrt((O.x-pM0.x)^2+(O.y-pM0.y)^2)+sqrt((OO.x-pM0.x)^2+(OO.y-pM0.y)^2)
                distM1 = sqrt((O.x-pM1.x)^2+(O.y-pM1.y)^2)+sqrt((OO.x-pM1.x)^2+(OO.y-pM1.y)^2)
                if distM0 < distM1
                    vg[k].sInt.x[i] = 0
                    z1 = z1 - (vg[k].sPrj.x[i] * c1[i])
                    z2 = z2 - (vg[k].sPrj.x[i] * c2[i])
                else
                    vg[k].sInt.x[i] = 1
                    z1 = z1 - (vg[k].sPrj.x[i] * c1[i]) + c1[i]
                    z2 = z2 - (vg[k].sPrj.x[i] * c2[i]) + c2[i]
                end
            end
        end

        # calcule la performance de la solution entiere
        vg[k].sInt.y[1] += vg[k].sInt.x[i] * c1[i]
        vg[k].sInt.y[2] += vg[k].sInt.x[i] * c2[i]
    end

    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    end      

end


# ==============================================================================
# arrondi la solution correspondant au generateur (pas d'historique donc)
# version avec voisinage et selection d'un voisin selon distance L1 avec generateur
function roundingSolutionNew23!(vg::Vector{tGenerateur}, k::Int64, c1::Array{Int,1}, c2::Array{Int,1}, d::tListDisplay)

    nbvar = length(vg[k].sInt.x)
    nbgen = size(vg,1)
    lstIdxFrac =(Int64)[]

    vg[k].sInt.y[1] = 0
    vg[k].sInt.y[2] = 0

    # --------------------------------------------------------------------------
    # identifie les variables fractionnaires et marquage par valeur sentinelle -1
    for i in 1:nbvar
        if  isapprox(vg[k].sPrj.x[i] , 0.0, atol=1e-3)
            vg[k].sInt.x[i] = 0
        elseif isapprox(vg[k].sPrj.x[i] , 1.0, atol=1e-3)
            vg[k].sInt.x[i] = 1
            vg[k].sInt.y[1] += c1[i] # Comptabilise la contribution des var a 1
            vg[k].sInt.y[2] += c2[i]
        else
            vg[k].sInt.x[i] = -1
            push!(lstIdxFrac,i) # marquage par valeur sentinelle
        end
    end

    # --------------------------------------------------------------------------
    # Selectionne les points pour le cone pointe sur le generateur k (pCour) et ouvert vers Y
    pPrec, pCour, pSuiv = selectionPoints(vg, k)

    # --------------------------------------------------------------------------
    # elabore le point pointeur du cone ouvert vers L
    pC = elaborePointConeOuvertversL(vg, k, pPrec, pSuiv)

    # --------------------------------------------------------------------------
    # Arrondi les valeurs non-entieres d'une solution fractionnaire
    # 1) applique un arrondi qui fixe 1 variable (la premiere venue) par iteration

    # calcule en differentiel par rapport a l'image de la solution fractionnaire
    z1 = vg[k].sPrj.y[1]
    z2 = vg[k].sPrj.y[2]
    pM = tPoint( z1 , z2 )

#    if length(vg[k].sPrj.x) ≤ 20 @show vg[k].sPrj.x end
#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end

    nbVarNonEntiere = length(lstIdxFrac)
    # traite toutes les variables non entieres

#=
    @printf("\n\n")
    @show c1
    @show c2
    @printf("VARIABLES FIXEES : \n")
    @show vg[k].sPrj.x
    @show vg[k].sInt.x
    @show vg[k].sInt.y
    @show pM
    @show lstIdxFrac
=#
    verbose ? @printf("  %2dR : [ %8.2f , %8.2f ] ", k, z1, z2) : nothing

    # traitement de toutes les variables fractionnaires une a une
    for i= 1:nbVarNonEntiere
#        @printf("\n")

        distMin = Inf; ind = -1; val = -1
        # voisinage a 1 variable sur toutes les variables
        for iFrac in [i for i=1:nbvar if vg[k].sInt.x[i] == -1]

            # Evalue x[iFrac]=0
            # defalque la contribution fractionnaire de la variable aux objectifs => pose x[iFrac]=0
            pM = tPoint( z1 - c1[iFrac] * vg[k].sPrj.x[iFrac] , z2 - c2[iFrac] * vg[k].sPrj.x[iFrac])
            if inCone1VersZ(pCour, pSuiv, pPrec, pM)
                distL1 = abs(vg[k].sRel.y[1]-pM.x) + abs(vg[k].sRel.y[2]-pM.y)
                if distL1 < distMin
                    distMin = distL1; ind = iFrac; val = 0
                end
            end

            # Evalue x[iFrac]=1
            # ajoute la contribution entiere de la variable aux objectifs => pose x[iFrac]=1
            pM = tPoint( z1 - c1[iFrac] * vg[k].sPrj.x[iFrac] + c1[iFrac], z2 - c2[iFrac] * vg[k].sPrj.x[iFrac] + c2[iFrac])
            distL1 = abs(vg[k].sRel.y[1]-pM.x) + abs(vg[k].sRel.y[2]-pM.y)
            if distL1 < distMin
                distMin = distL1; ind = iFrac; val = 1
            end
        end

        # applique le voisin le plus interessant trouve
        vg[k].sInt.x[ind] = val
        vg[k].sInt.y[1] += vg[k].sInt.x[ind] * c1[ind] * val
        vg[k].sInt.y[2] += vg[k].sInt.x[ind] * c2[ind] * val
        z1 = z1 - (vg[k].sPrj.x[ind] * c1[ind]) + c1[ind] * val
        z2 = z2 - (vg[k].sPrj.x[ind] * c2[ind]) + c2[ind] * val
    #    @printf("  x[%2d]=%2d |  : [ %12.5f , %12.5f ] \n",ind, val, z1, z2)

    end

#    if length(vg[k].sInt.x) ≤ 20 @show vg[k].sInt.x end
    @printf("→ #round : %4d → [ %5d , %5d ] ", nbVarNonEntiere, vg[k].sInt.y[1], vg[k].sInt.y[2])

    # archive le point obtenu pour les besoins d'affichage
    if generateurVisualise == -1 
        # archivage pour tous les generateurs
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    elseif generateurVisualise == k
        # archivage seulement pour le generateur k
        push!(d.XInt,vg[k].sInt.y[1])
        push!(d.YInt,vg[k].sInt.y[2])
    end      

end

# ==================================================================Rounding Interface===========================================================
# Interface intended in manipulating and testing differents configurations of GM more easily
const configurationRounding::Dict{Int,Function} = Dict{Int,Function}(
                                                                1 => roundingSolution!,
                                                                2 => roundingSolutionNew23!,
                                                                3 => roundingSolutionnew24!
                                                            )

function interface_roundingSolution!(vg::Vector{tGenerateur},k::Int64,c1::Vector{Int},c2::Vector{Int},d::tListDisplay;CHOICE::Int=CHOICE_ROUNDING)
    return configurationRounding[CHOICE](vg,k,c1,c2,d)
end
