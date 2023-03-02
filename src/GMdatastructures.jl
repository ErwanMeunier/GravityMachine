import Base.== # to allow the extention of this method for tPoint
# types ------------------------------------------------------------------------

# type corresponding to a solution
mutable struct tSolution{T}
    x :: Vector{T}                # vector variables x (1..n)
    y :: Vector{T}                # vector outcomes  y (1..p)
end

# type corresponding to a point generator
mutable struct tGenerateur
    sRel :: tSolution{Float64}    # initial relaxation
    sInt :: tSolution{Int64}      # integer solution
    sPrj :: tSolution{Float64}    # projected solution
    sFea :: Bool                  # indicate if sInt is feasible or not
end

# type of a point (x,y) in the objective space
mutable struct tPoint
    x :: Float64
    y :: Float64
end

function ==(a::tPoint, b::tPoint) # Extending equality relation for the tPoint
    return (a.x == b.x) && (a.y == b.y)
end

# type grouping the lists of points for dysplaying purposes
mutable struct tListDisplay
    # generators found considering an ϵ-constraint on f1
    xLf1  :: Vector{Float64};  yLf1  :: Vector{Float64} # liste des points (x,y) relaches -> Sur f1 pendant la génération par ϵ-contraintes
    # generators found considering an ϵ-constraint on f2
    xLf2  :: Vector{Float64};  yLf2  :: Vector{Float64} # liste des points (x,y) relaches -> Sur f2 pendant la génération par ϵ-contraintes
    # All the generators
    xL    :: Vector{Float64};  yL    :: Vector{Float64} # liste des points (x,y) relaches
    #xLf1Improved :: Vector{Float64};  yLf1Improved  :: Vector{Float64} # liste des points (x,y) relaches améliorés #Recent improvement
    #xLf2Improved  :: Vector{Float64};  yLf2Improved  :: Vector{Float64} # liste des points (x,y) relaches améliorés #Recent improvement 
    xLImproved   :: Vector{Float64};  yLImproved    :: Vector{Float64} # liste des points (x,y) relaches améliorés #Recent improvement
    XInt  :: Vector{Int64};    YInt  :: Vector{Int64}   # liste des points (x,y) entiers
    XProj :: Vector{Float64};  YProj :: Vector{Float64} # liste des points (x,y) projetes
    XFeas :: Vector{Int64};    YFeas :: Vector{Int64}   # liste des points (x,y) admissibles
    XPert :: Vector{Int64};    YPert :: Vector{Int64}   # liste des points (x,y) perturbes
end


# ==============================================================================
# Initialisation structure donnees contenant tous les generateurs

function allocateDatastructure(nbgen::Int64, nbvar::Int64, nbobj::Int64)

    verbose ? println("\n  → Allocation memoire pour ",nbgen," generateurs\n") : nothing

    vg::Vector{tGenerateur} = Vector{tGenerateur}(undef, nbgen)
    for k = 1:nbgen
        vg[k] = tGenerateur(
                            tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                            tSolution{Int64}(zeros(Int64,nbvar),zeros(Int64,nbobj)),
                            tSolution{Float64}(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                            false
                            )
    end
    return vg
end