# Erwan Meunier - Nantes Université - January 2023 - erwan.meunier@etu.univ-nantes.fr

using SparseArrays
using Random
using HypothesisTests # Useful in benchmarking sampleInSimplex
include("GMdatastructures.jl")

# Overriding isequal and hash functions for tPoint datastructure using in Set
Base.isequal(a::tPoint,b::tPoint) = (a.x == b.x) & (a.y == b.y)
Base.hash(a::tPoint, h::UInt) = hash(a.y, hash(a.x, hash(tPoint, h)))
#Base.==(a::tPoint, b::tPoint) = Base.isequal(a,b)

###
Random.seed!(314)

#=
If horizontal = true then the ranges stores the vertical ranges of integrals points along the horizontal axis.
=#
struct integralTriangle
    ranges::Dict{Int64,Tuple{Int64,Int64}}
    tranposed::Bool
    nbintegrals::Int64
end


#= Computes the area of the triangle delimited by a, b and c using the Heron's formula.
-vertices must be of length 3
https://en.wikipedia.org/wiki/Heron%27s_formula
=#
function areaTriangle(a::tPoint,b::tPoint,c::tPoint)::Float64 # OK
    ab::Float64 = sqrt((a.x - b.x)^2 + (a.y - b.y)^2)
    bc::Float64 = sqrt((b.x - c.x)^2 + (b.y - c.y)^2)
    ca::Float64 = sqrt((c.x - a.x)^2 + (c.y - a.y)^2)
    s::Float64 = (ab+bc+ca)/2

    return sqrt(abs(s*(s-ab)*(s-bc)*(s-ca)))
end

#=
Compares the sum of areas of each sub-triangle with the area of the triangle delimited by vertices.
If both areas are equals, then the d is in the considered triangle.
=#
function inTriangle(vertices::Vector{tPoint},d::tPoint)::Bool # OK
    A::Float64 = areaTriangle(vertices[1],vertices[2],vertices[3])
    # Area of each sub-triangles
    A1::Float64 = areaTriangle(vertices[1],vertices[2],d)
    A2::Float64 = areaTriangle(vertices[2],vertices[3],d)
    A3::Float64 = areaTriangle(vertices[3],vertices[1],d)
    return isapprox(A,A1+A2+A3,atol=10^-6)
end

# vertices must be of length 3
function extractIntegralTriangle(vertices::Vector{tPoint})::integralTriangle #OK
    # Centering the triangle
    minX::Int64 = Int64(floor(minimum(map(a->a.x,vertices))))
    minY::Int64 = Int64(floor(minimum(map(a->a.y,vertices))))
    newvertices::Vector{tPoint} = collect(map(a->tPoint(a.x-minX,a.y-minY),vertices))
    l::Int64 = Int64(ceil(maximum(map(a->a.x, newvertices))))
    h::Int64 = Int64(ceil(maximum(map(a->a.y, newvertices))))
    transposed::Bool = false
    sizeX::Int64 = l
    sizeY::Int64 = h
    if h <= l # ranges will be represented as horizontal ranges for each integer between 0 and ceil(l)
        sizeX = h
        sizeY = l
        transposed = true
    end
    ranges::Dict{Int64,Tuple{Int64,Int64}} = Dict{Int64,Tuple{Int64,Int64}}()
    sizehint!(ranges,sizeX)
    nbint::Int64 = 0
    # Transposition function
    tr(point::tPoint) = transposed ? tPoint(point.y,point.x) : point
    # Each integral point candidate is tested for its belonging to the triangle
    for i = 0:sizeX
        belonging::Vector{Bool} = [inTriangle(newvertices,tr(tPoint(i,j))) for j = 0:sizeY]
        if reduce(|,belonging) # at least one integral point to create an interval
            lb::Int64 = findfirst(belonging)-1 # because the index starts from 0
            ub::Int64 = lb + 1
            while inTriangle(newvertices,tr(tPoint(i,ub)))
                ub += 1
            end
            ranges[(transposed ? i+minY : i+minX)]=(lb,ub-1)
            nbint += (ub-1)-lb+1
        end
    end
    for k in keys(ranges) # translating interval
        ranges[k]= transposed ? (ranges[k][1]+minX, ranges[k][2]+minX) : (ranges[k][1]+minY,ranges[k][2]+minY)
    end
    nbint==0 ? println("WARNING: This triangle : ", vertices, " does not admit any integral point") : nothing
    return integralTriangle(ranges,transposed,nbint)
end

# Returns one and only one integral point which belongs to the given triangle, at random (uniformly).
function sampleInIntegralTriangle(triangle::integralTriangle)::tPoint
    keyset::Vector{Int64} = [i for i in keys(triangle.ranges)]
    #println(keyset)
    draw::Int64 = rand(1:triangle.nbintegrals)
    #println(draw)
    p::Int64 = 1
    k::Int64 = 1
    notFound::Bool = true
    x::Int64, y::Int64 = 0, 0 # sampled point to be returned
    while notFound # no need to check for the bounding of k because the drawing is "number of intergral points"-wise
        itr::Int64 = 0
        lb::Int64 = triangle.ranges[keyset[k]][1] # lower bound
        ub::Int64 = triangle.ranges[keyset[k]][2] # upper bound
        while notFound & ((lb+itr)<=ub)
            if (p+itr)==draw
                x, y = k, itr+lb
                notFound = false
            end
            itr+=1
        end
        p += itr
        k+=1
    end
    return triangle.tranposed ? tPoint(y,x) : tPoint(x,y)
end

#=
Return nbp all different integral points at Random, in the simplex who is defined according to its vertices.
The sampling follow the Uniform Law in the number of integral points inside the simplex.
If nbp > number of integral points, the maximum number of integral points is returned.
ATTENTION: 
- For now, this function only handles triangles.
- This method uses a discounted draw. Then, nbp MUST BE small (less than the halved number of integral points).
=#
function sampleInSimplex(vertices::Vector{tPoint}, nbp::Int64=1)::Vector{tPoint}
    @assert length(vertices)==3 "ERROR: This method only handles triangles -> Sampling aborted" # For now...
    triangle::integralTriangle = extractIntegralTriangle(vertices)
    #@assert triangle.nbintegrals>0 "ERROR: This triangle has not any integral point -> Sampling aborted"
    result::Vector{tPoint} = Vector{tPoint}()
    if nbp >= triangle.nbintegrals # All integral points are returned
        println("WARNING: The number of points to be sampled is greater or equal than the number of integral points in the simplex")
        result = [tPoint(x,y) for x in keys(triangle.ranges) for y in triangle.ranges[x][1]:triangle.ranges[x][2]]
        # Above: resp. lower bound and upper bound <=> triangle.ranges[x][1]:triangle.ranges[x][2]
    else
        sampledPoints::Set{tPoint} = Set{tPoint}()
        while length(sampledPoints)<nbp
            #println(sampledPoints)
            push!(sampledPoints,sampleInIntegralTriangle(triangle))
        end
        result = [tp for tp in sampledPoints]
    end
    # If needed -> Transposition of points
    if triangle.tranposed
        map!(tp->tPoint(tp.y,tp.x),result,result)
    end
    return result
end

function testUniformity(x::tPoint,y::tPoint,z::tPoint)
    triangle = extractIntegralTriangle([x,y,z])
    test = Dict{tPoint,Int64}()
    nbsamp = 100000
    for i = 1:nbsamp
        samp = sampleInSimplex([x,y,z])[1]
        if haskey(test,samp)
            test[samp] += 1
        else
            test[samp] = 1
        end
    end
    post = [test[k] for k in keys(test)]
    println("Frequencies: ", post ./ nbsamp)
    ChisqTest([test[k] for k in keys(test)])
end

function main()
    x = tPoint(1.,1.)
    y = tPoint(0.,6.)
    z = tPoint(6.,3.)
    #= Instance with warnings
    x = tPoint(0.,0.1)
    y = tPoint(0.1,0.2)
    z = tPoint(0.3,0.4)
    =#
    sampleInSimplex([x,y,z])
    #testUniformity(x,y,z)
end

