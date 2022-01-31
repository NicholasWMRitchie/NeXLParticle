
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
    # Find the nearest neighbors for pts2c in pts1c
    tr = KDTree(pts2c)
    i_knn = collect(Iterators.flatten(knn(tr, pts1c, 1)[1]))
    # Reorder pts2c to best correspond to pts1's order
    pts2cr = pts2c[i_knn]
    # Turn Vector{Vector{<:AbstractFloat}} into Matrix{<:AbstractFloat}
    m2(vs) = transpose(reshape(collect(Iterators.flatten(vs)), length(vs[1]), length(vs)))
    # Solve the orthogonal Procrustes problem
    f = svd(m2(pts1c) * transpose(m2(pts2cr)))
    Ω = f.U * f.Vt
    # Now translate and rotated the reordered data towards `pts1`
    return Translation(com1).(Ω * pts2cr)
end

"""
    icp(pts1::AbstractVector{T}, pts2::AbstractVector{T}; maxiter=10, tol=0.9)

Implements the iterative closest point algorithm.  Transforms `pts2` through rotations and 
translations to come as close as possible to `pts1` using a least-squares metric.

The intention is that `pts1` is the super-set of points in `pts2`.  When aligned most of
the points in `pts2` will match up with points in `pts1`.
"""
function icp(pts1::AbstractVector{T}, pts2::AbstractVector{T}; maxiter=10, tol=0.9) where { T <: AbstractVector{<:AbstractFloat}}
    function err(pts1, pts2)
        # Associate each point in pts2 with the closest point in pts1 
        pts2r = pts2[collect(Iterators.flatten(knn(KDTree(pts2), pts1, 1)[1]))]
        # Compute the sum Euclidean distance
        sum(norm(p1 - p2) for (p1, p2) in zip(pts1, pts2r))
    end
    preverror, next = err(pts1, pts2), pts2
    for _ in 1:maxiter
        next = icp_iteration(pts1, next)
        nexterror = err(pts1, next)
        if nexterror > tol*preverror
            return next
        end
        preverror = nexterror
    end
    @warn "Not converging after $maxiter steps in icp(...)."
    return next
end

end