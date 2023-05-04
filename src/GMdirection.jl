# ==============================================================================
# Calcule la direction d'interet du nadir vers le milieu de segment reliant deux points generateurs

# Main functions
function calculerDirections(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})::Tuple{Vector{Float64},Vector{Float64}}
    # function calculerDirections(L, vg::Vector{tGenerateur})
 
     nbgen = size(vg,1)
     for k in 2:nbgen
 
         n1 = L[end].y[1]
         n2 = L[1].y[2]
 
         x1,y1 = vg[k-1].sRel.y[1], vg[k-1].sRel.y[2]
         x2,y2 = vg[k].sRel.y[1], vg[k].sRel.y[2]
         xm=(x1+x2)/2.0
         ym=(y1+y2)/2.0
         Δx = abs(n1-xm)
         Δy = abs(n2-ym)
         λ1 =  1 - Δx / (Δx+Δy)
         λ2 =  1 - Δy / (Δx+Δy)
         @printf("  x1= %7.2f   y1= %7.2f \n",x1,y1)
         @printf("  x2= %7.2f   y2= %7.2f \n",x2,y2)
         @printf("  Δx= %7.2f    Δy= %7.2f \n",Δx,Δy)
         @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1,λ2)
         #=plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+")
         annotate("",
                  xy=[xm;ym],# Arrow tip
                  xytext=[n1;n2], # Text offset from tip
                  arrowprops=Dict("arrowstyle"=>"->"))
         println("")=#
     end
     return λ1, λ2
 end
 
 # ==============================================================================
 # Calcule la direction d'interet du nadir vers un point generateur
 function calculerDirections2(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})::Tuple{Vector{Float64},Vector{Float64}} # OK 
     #function calculerDirections2(L, vg::Vector{tGenerateur})
     #println(L)
     nbgen = size(vg,1)
     λ1=Vector{Float64}(undef, nbgen)
     λ2=Vector{Float64}(undef, nbgen)
     #println("VG : ", vg)
     for k in 1:nbgen
 
         # Calcul du nadir
         n1 = L[end].y[1] # max de la valeur du premier objectif
         n2 = L[1].y[2] # max de la valeur du deuxième objectif
 
         xm=vg[k].sRel.y[1] # coordonnée z1 du générateur
         ym=vg[k].sRel.y[2] # coordonnée z2 du générateur
         Δy1 = n1-xm # abs(n1-xm)
         Δy2 = n2-ym # abs(n2-ym)
         λ2[k] =  1 - Δy1 / (Δy1+Δy2) # = Δy/(Δx + Δy)
         λ1[k] =  1 - Δy2 / (Δy1+Δy2) # = Δx/(Δx + Δy) = 1 - λ1[k]
         @printf("  k= %3d   ",k)
         @printf("  xm= %7.2f   ym= %7.2f ",xm,ym)
         @printf("  Δx= %8.2f    Δy= %8.2f ",Δy1,Δy2)
         @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1[k],λ2[k])
         #=if generateurVisualise == -1 
             # affichage pour tous les generateurs
            graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
            #=graphic ? annotate("",k,
                      xy=[xm;ym],# Arrow tip
                      xytext=[n1;n2], # Text offset from tip
                      arrowprops=Dict("arrowstyle"=>"->")) : nothing=#
         elseif generateurVisualise == k
             # affichage seulement pour le generateur k
             graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
             graphic ? annotate("",
                      xy=[xm;ym],# Arrow tip
                      xytext=[n1;n2], # Text offset from tip
                      arrowprops=Dict("arrowstyle"=>"->")) : nothing  
         end
         =#
         #println("")
     end
     return λ1, λ2
 end
 
 # ERWAN
 # TODO
 #= Calcul aléatoirement une direction basée sur une solution entière (admissible ou non) dans le cône de source nadir 
 pointant vers vg[k-1] et vg[k+1] 
 =# 
 
 function calculerDirections3(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})::Tuple{Vector{Float64},Vector{Float64}}
     nbgen::Int64 = size(vg,1)
     λ1::Vector{Float64} = Vector{Float64}(undef, nbgen)
     λ2::Vector{Float64} = Vector{Float64}(undef, nbgen)
 
     n1 = L[end].y[1] # max value found of the first objective
     n2 = L[1].y[2]  # max value found of the second objective
     nadir = tPoint(n1,n2)
     
     for k in 1:nbgen
         if (k==1) | (k==nbgen) # Edge generator
             xm=vg[k].sRel.y[1] # coordonnée z1 du générateur
             ym=vg[k].sRel.y[2] # coordonnée z2 du générateur
             Δx = abs(n1-xm)
             Δy = abs(n2-ym)
             λ1[k] =  1 - Δx / (Δx+Δy)
             λ2[k] =  1 - Δy / (Δx+Δy)
         else # Inner generators
             xuk = vg[k-1].sRel.y[1]
             yuk = vg[k-1].sRel.y[2]
             upperk = tPoint(xuk,yuk)
 
             xlk = vg[k+1].sRel.y[1]
             ylk = vg[k+1].sRel.y[2]
             lowerk = tPoint(xlk,ylk)
 
             xm=vg[k].sRel.y[1] # first coordinate of the k-th generator
             ym=vg[k].sRel.y[2] # second coordinate of the k-th generator
             Δx = abs(n1-xm)
             Δy = abs(n2-ym)
 
             sampling::Vector{tPoint} = sampleInSimplex([nadir,lowerk,upperk])
 
             if !isempty(sampling) # An integral point random driven vector is computed
                 m = sampling[1]
                 Δx = abs(n1-m.x)
                 Δy = abs(n2-m.y)
             end
         end
 
         λ1[k] = 1 - Δx / (Δx+Δy)
         λ2[k] = 1 - Δy / (Δx+Δy)
 
         @printf("  k= %3d   ",k)
         @printf("  xm= %7.2f   ym= %7.2f ",xm,ym)
         @printf("  Δx= %8.2f    Δy= %8.2f ",Δx,Δy)
         @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1[k],λ2[k])
         #=if generateurVisualise == -1 
             # affichage pour tous les generateurs
             graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
             graphic ? annotate("",
                      xy=[xm;ym],# Arrow tip
                      xytext=[n1;n2], # Text offset from tip
                      arrowprops=Dict("arrowstyle"=>"->")) : nothing      
         elseif generateurVisualise == k
             # affichage seulement pour le generateur k
             graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
             graphic ? annotate("",
                      xy=[xm;ym],# Arrow tip
                      xytext=[n1;n2], # Text offset from tip
                      arrowprops=Dict("arrowstyle"=>"->")) : nothing
         end =#
     end
     return λ1, λ2
 end
 
function calculerDirections4(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur})
    nbgen::Int64 = size(vg,1)
    λ1::Vector{Float64} = Vector{Float64}(undef, nbgen)
    λ2::Vector{Float64} = Vector{Float64}(undef, nbgen)

    globalNadir1 = L[end].y[1] # max value found of the first objective
    globalNadir2 = L[1].y[2]  # max value found of the second objective
    nadirs::Vector{tPoint} = computeLocalNadirs(vg,L)
    
    nbgen = length(vg)
    λ1=Vector{Float64}(undef, nbgen)
    λ2=Vector{Float64}(undef, nbgen)
    for k in 1:nbgen
        n1 = nadirs[k].x
        n2 = nadirs[k].y

        xm=vg[k].sRel.y[1] # coordonnée z1 du générateur
        ym=vg[k].sRel.y[2] # coordonnée z2 du générateur
        Δx = abs(n1-xm)
        Δy = abs(n2-ym)
        λ1[k] =  1 - Δx / (Δx+Δy)
        λ2[k] =  1 - Δy / (Δx+Δy)
        @printf("  k= %3d   ",k)
        @printf("  xm= %7.2f   ym= %7.2f ",xm,ym)
        @printf("  Δx= %8.2f    Δy= %8.2f ",Δx,Δy)
        @printf("  λ1= %6.5f    λ2= %6.5f \n",λ1[k],λ2[k])
        #=if generateurVisualise == -1 
            # affichage pour tous les generateurs
            graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
            graphic ? annotate("",k,
                    xy=[xm;ym],# Arrow tip
                    xytext=[n1;n2], # Text offset from tip
                    arrowprops=Dict("arrowstyle"=>"->")) : nothing
        elseif generateurVisualise == k
             # affichage seulement pour le generateur k
            graphic ? plot(n1, n2, xm, ym, linestyle="-", color="blue", marker="+") : nothing
            graphic ? annotate("",
                    xy=[xm;ym],# Arrow tip
                    xytext=[n1;n2], # Text offset from tip
                    arrowprops=Dict("arrowstyle"=>"->")) : nothing  
        end=#
         #println("")
    end
    return λ1, λ2
end

const configurationComputeDirections::Dict{Int,Function} = Dict{Int,Function}(
                                                                1 => calculerDirections,
                                                                2 => calculerDirections2,
                                                                3 => calculerDirections3,
                                                                4 => calculerDirections4                                                                
                                                            )

function interface_computeDirections(L::Vector{tSolution{Float64}}, vg::Vector{tGenerateur};CHOICE::Int=CHOICE_COMPUTEDIRECTIONS)
    return configurationComputeDirections[CHOICE](L,vg)
end

function printlambdas(v1,v2)
    λ1 = fill(0.,length(v1))
    λ2 = fill(0.,length(v2))
    n1 = maximum(v1)
    n2 = maximum(v2)
    for k in 1:length(v1)
        xm=v1[k] # coordonnée z1 du générateur
        ym=v2[k] # coordonnée z2 du générateur
        Δy1 = n1-xm # abs(n1-xm)
        Δy2 = n2-ym # abs(n2-ym)
        λ2[k] =  1 - Δy1 / (Δy1+Δy2) # = Δy/(Δx + Δy)
        λ1[k] =  1 - Δy2 / (Δy1+Δy2) # = Δx/(Δx + Δy) = 1 - λ1[k]
    end
    result = [(round(2*λ1[k]+v1[k],digits=2), round(2*λ2[k]+v2[k],digits=2)) for k in 1:length(v1)]
    for k in 1:length(v1)
        println(round(v1[k],digits=2),"/",round(v2[k],digits=2),"/",result[k][1],"/",result[k][2])
    end
end