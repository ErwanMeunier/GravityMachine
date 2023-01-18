using SparseMatrixCSC
include("./src/GMdatastructures.jl")

struct integralTriangle
    range::Vector{Tuple{Int64,Int64}}
end


# vertices must be of length 3
# With the area method
function inTriangle(vertices::Vector{tPoint},d::tPoint)::Bool
    minX::Float64 = minimum(map(w->w[1],vertices))
    minY::Float64 = minimum(map(w->w[2],vertices))
    A = 
    A1 = 
    A2 = 
    A3 =
    h = 
end

# vertices must be of length 3
function extractIntegralTriangle(vertices::Vector{tPoint})::integralTriangle
    @assert length(vertices==3) "Triangle has 3 vertices in 2-Dimensions"
    minX::Int64 = Int64(floor(minimum(map(w->w[1],vertices))))
    minY::Int64 = Int64(floor(minimum(map(w->w[2],vertices))))

end