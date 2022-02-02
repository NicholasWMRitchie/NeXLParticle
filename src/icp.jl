using CoordinateTransformations, Rotations
using NearestNeighbors
using LinearAlgebra


centerofmass(pts::AbstractVector{T}) where { T <: AbstractVector{<:AbstractFloat} } = #
    [mean(p -> p[i], pts) for i in 1:length(pts[1])]

"""
    icp_iteration(pts1::AbstractVector{<:AbstractVector{2, T}}, pts2::AbstractVector{<:AbstractVector{2, T}}; orthonormal=false) where { T<: AbstractFloat }

Perform a single iteration of the iterative closest point algorithm using the orthogonal Procrustes algorithm.

Returns the `pts2` rotated and translated to most closely match `pts1`. 
"""
function icp_iteration(pts1::AbstractVector{T}, pts2::AbstractVector{T}) where { T <: AbstractVector{<:AbstractFloat}}
    # Center both pts1 and pts2 on the origin
    com1, com2 = centerofmass(pts1), centerofmass(pts2)
    pts1c, pts2c = Translation(-com1).(pts1), Translation(-com2).(pts2)
    # Find the nearest point in pts1c for each point in pts2c
    pts1cr = map(i->pts1c[i[1]], knn(KDTree(pts1c), pts2c, 1)[1])
    @assert length(pts1cr)==length(pts2c)
    # Turn Vector{Vector{<:AbstractFloat}} into Matrix{<:AbstractFloat}
    m2(vs) = transpose(reshape(collect(Iterators.flatten(vs)), length(vs[1]), length(vs)))
    # Solve the orthogonal Procrustes problem
    f = svd(m2(pts1cr) * transpose(m2(pts2c)))
    Ω = f.U * f.Vt
    # Now translate and rotated the reordered data towards `pts1`
    return Translation(com1).(Ω * pts2c)
end

"""
    icp(pts1::AbstractVector{T}, pts2::AbstractVector{T}; maxiter=10, tol=0.9)

Implements the iterative closest point algorithm.  Transforms `pts2` through rotations and 
translations to come as close as possible to `pts1` using a least-squares metric.

The intention is that `pts1` is the super-set of points in `pts2`.  When aligned most of
the points in `pts2` will match up with points in `pts1`.
"""
function icp(pts1::AbstractVector{T}, pts2::AbstractVector{T}; maxiter=10, tol=0.99) where { T <: AbstractVector{<:AbstractFloat}}
    initialerror = icperror(pts1, pts2)
    next, res, minerror, pos = pts2, pts2, initialerror, 0
    @show initialerror
    for _ in 1:maxiter
        next = icp_iteration(pts1, next)
        nexterror = icperror(pts1, next)
        if nexterror < minerror
            # Meets tolerance
            (nexterror > tol * minerror) && return next
            res, minerror, pos = next, nexterror, 0
        else
            # Increasing error
            (minerror!=initialerror) && ((pos+=1)==2) && return res
        end
    end
    minerror == initialerror && @warn "No improvement after $maxiter steps in icp(...). Returning the initial points."
    return res
end
"""
    icperror(pts1::AbstractVector{T}, pts2::AbstractVector{T}) where { T <: AbstractVector{<:AbstractFloat}}

Measures the sum distance between the points in `pts2` and the point closest to the point in `pts1`.
"""
icperror(pts1::AbstractVector{T}, pts2::AbstractVector{T}) where { T <: AbstractVector{<:AbstractFloat}} =
    sum(i->i[1], knn(KDTree(pts1), pts2, 1)[2])