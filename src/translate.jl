using CoordinateTransformations
using LinearAlgebra
using Optim


const XYPosition = SArray{Tuple{2},Float64,1,2}

"""
    translate(oldPts::AbstractVector{XYPosition}, newPts::AbstractVector{XYPosition})::AffineMap

Determines the optimal AffineMap to transform the X-Y coordinates in `oldPts` into the X-Y coordinates in `newPts`.
The coordinates in `oldPts` and `newPts` are assumed to represent the same features in the same order and thus
must contain the same number of X-Y coordinates.
"""
function translate(oldPts::AbstractVector{XYPosition}, newPts::AbstractVector{XYPosition})::AffineMap
    @assert length(oldPts)==length(newPts) "The number of old and new points must match."
    if length(oldPts)==1
        offset, scale, rotatio = SA[ newPts[1] .- oldPts[1] ], SA[1.0, 1.0], 0.0
    else # Compute a first estimate from points 1 and 2
        offset = 0.5*(newPts[2] .+ newPts[1]) .- 0.5*(oldPts[2] .+ oldPts[1])
        dold, dnew = oldPts[2] .- oldPts[1], newPts[2] .- newPts[1]
        scale = SA[ norm(dnew)/norm(dold), norm(dnew)/norm(dold) ]
        rotation = acos((dold⋅dnew)/(norm(dold)*norm(dnew)))
    end
    res = AffineMap(SA[ cos(rotation) -sin(rotation) ; sin(rotation) cos(rotation) ] * SA[ scale[1] 0.0; 0.0 scale[2]], offset)
    if length(oldPts)>2
        # Simplex optimize to discover the minimum offset, scale and rotation over all point pairs
        x0 = [offset..., scale..., rotation]
        function f(x)
            res = AffineMap(SA[ cos(x[5]) -sin(x[5]) ; sin(x[5]) cos(x[5]) ] * SA[ x[3] 0.0; 0.0 x[4] ], SA[x[1], x[2]])
            return sum(norm(res(oldPt) .- newPt) for (oldPt, newPt) in zip(oldPts, newPts))
        end
        xp = Optim.minimizer(optimize(f, x0, Optim.NelderMead()))
        res = AffineMap(SA[ cos(xp[5]) -sin(xp[5]) ; sin(xp[5]) cos(xp[5]) ] * SA[ xp[3] 0.0; 0.0 xp[4] ], SA[xp[1], xp[2]])
    end
    return res
end
