using SparseArrays
include("./src/GMdatastructures.jl")


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
    size::Float64 = l
    transposed = false
    sizeX = l
    sizeY = h
    if h <= l # ranges will be represented as horizontal ranges for each integer between 0 and ceil(l)
        sizeX = h
        sizeY = l
        transposed = true
    end
    ranges::Dict{Int64,Tuple{Int64,Int64}} = Dict{Int64,Tuple{Int64,Int64}}()
    sizehint!(ranges,sizeX)
    nbint::Int64 = 0
    #TODO TRANSPOSITION
    for i = 0:sizeX
        println("testLB : ", [inTriangle(newvertices,tPoint(i,j)) for j = 0:sizeY])
        lb::Int64 = findfirst([inTriangle(newvertices,tPoint(i,j)) for j = 0:sizeY])-1 # because the index starts from 0
        ub::Int64 = lb + 1
        while inTriangle(newvertices,tPoint(i,ub))
            ub += 1
        end
        ranges[i]=(lb,ub-1)
        nbint += (ub-1)-lb+1
    end
    #TODO UPDATE INTEVRAL
    for k in keys(ranges)
        ranges[k]=tPoint(ranges[k][+minX,ranges[k].y+minY)
    end
    return integralTriangle(ranges,transposed,nbint)
end

function main()
    x = tPoint(1.,1.5)
    y = tPoint(3.,4.)
    z = tPoint(0.,6.)
    d = tPoint(6.,0.)
    println(inTriangle([x,y,z],d))
    println(extractIntegralTriangle([x,y,z]))
end