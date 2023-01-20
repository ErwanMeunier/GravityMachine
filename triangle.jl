using SparseArrays
using Random
include("./src/GMdatastructures.jl")

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
    return sqrt(s*(s-ab)*(s-bc)*(s-ca))
end

#=
Compares the sum of areas of each sub-triangle with the area of the triangle delimited by vertices.
If both areas are equals, then the d is in the considered triangle.
=#
function inTriangle(vertices::Vector{tPoint},d::tPoint)::Bool # OK
    @assert length(vertices)==3 "Triangle has 3 vertices in 2-Dimensions"
    A::Float64 = areaTriangle(vertices[1],vertices[2],vertices[3])
    # Area of each sub-triangles
    A1::Float64 = areaTriangle(vertices[1],vertices[2],d)
    A2::Float64 = areaTriangle(vertices[2],vertices[3],d)
    A3::Float64 = areaTriangle(vertices[3],vertices[1],d)
    return isapprox(A,A1+A2+A3,atol=10^-6)
end

# vertices must be of length 3
function extractIntegralTriangle(vertices::Vector{tPoint})::integralTriangle
    @assert length(vertices)==3 "Triangle has 3 vertices in 2-Dimensions"
    # Centering the triangle
    minX::Int64 = Int64(floor(minimum(map(a->a.x,vertices))))
    minY::Int64 = Int64(floor(minimum(map(a->a.y,vertices))))
    newvertices::Vector{tPoint} = collect(map(a->tPoint(a.x-minX,a.y-minY),vertices))
    println(newvertices)
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
        lb::Int64 = findfirst([inTriangle(newvertices,tr(tPoint(i,j))) for j = 0:sizeY])-1 # because the index starts from 0
        ub::Int64 = lb + 1
        while inTriangle(newvertices,tr(tPoint(i,ub)))
            ub += 1
        end
        ranges[i]=(lb,ub-1)
        nbint += (ub-1)-lb+1
    end
    #TODO UPDATE INTERVAL
    for k in keys(ranges)
        ranges[k]= transposed ? (ranges[k][1]+minY+1, ranges[k][2]+minX) : (ranges[k][1]+minX+1,ranges[k][2]+minY)
    end
    return integralTriangle(ranges,transposed,nbint)
end

function sampleInTriangle(triangle::integralTriangle)::tPoint
    keyset = keys(triangle.ranges)
    draw::Int64 = rand(1:triangle.nbintegrals)
    p::Int64 = 1
    k::Int64 = keyset[1]
    sampledPoint::tPoint = tPoint()
    while p!=draw
        itr::Int64 = 0
        lb::Int64 = triangle.ranges[k][1] # lower bound
        ub::Int64 = triangle.ranges[k][2]
        while ((base+itr)!=draw) && (itr+base)<=ub
            itr+=1
        end
        p += itr
    end
    return 
end

function main()
    x = tPoint(1.5,1.)
    y = tPoint(4.,3.)
    z = tPoint(6.,0.)
    d = tPoint(6.,0.)
    println(inTriangle([x,y,z],d))
    println(extractIntegralTriangle([x,y,z]))
end