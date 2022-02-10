# Given two particle data sets collected on the same sample but at different
# positions and rotations, find the offset and rotation which best aligns
# datasets.
using Random
using GeometryBasics: Rect2
using StaticArrays
using Rotations
using NearestNeighbors
using LinearAlgebra
using StatsBase
using Distributions
using CoordinateTransformations

"""
generate_ground_truth(rnd::Random, bounds::Rect2, npart::Integer, inside::Function=p->p in bounds)

Generate a simulated ground-truth particle data set.  This represents the "true" location of `npart` particles that are randomly 
distributed within `bounds` and for which `inside(pt)` evaluates true.
"""
function generate_ground_truth(rnd::AbstractRNG, bounds::Rect2, npart::Integer, inside::Function=p->p in bounds)
    @assert minimum(bounds.widths)>0.0 "Rectangle widths must be positive."
    sc = maximum(bounds.widths)
    res = map(Base.OneTo(npart)) do _ 
        pt = SVector{2}(sc .* rand(rnd, Float64, 2) .+ bounds.origin)
        while !inside(pt)
            pt = SVector{2}(sc .* rand(rnd, Float64, 2) .+ bounds.origin)
        end
        pt
    end            
    return res
end

"""
    measure_particles(rnd::AbstractRNG, pd::Vector{<:SVector{2}}, offset::SVector{2}, θ::Float64, frac::Float64, eps::Float64)

Sample a ground-truth particle data set in `pd`.  A randomized fraction `frac` of the data set `pd` is sampled and then the
positions are transformed according to an `offset` and a rotation by `θ`. In addition, a measurement jitter the scale of 
which is determiend by `eps` is added to each measured coordinate.
"""
function measure_particles(rnd::AbstractRNG, pd::Vector{<:SVector{2}}, offset::SVector{2}, θ::Float64, frac::Float64, eps::Float64)
    @assert frac > 0.0 "The fraction measured `frac` must be larger than 0.0."
    @assert frac <= 1.0 "The fraction measured `frac` must be less than or equal to 1.0."
    rot = RotMatrix{2}(θ)
    n = Distributions.Normal(0.0, eps)
    return shuffle!(map(filter(f->rand(rnd)<frac, pd)) do part
        rot*part + offset + rand(rnd, n, 2)
    end)
end

"""
    rough_align(ps1::AbstractVector{StaticVector{2, T}}, ps2::AbstractVector{StaticVector{2, T}}; nnsize=3, tol=0.001, corrtol=10.0) where { T<: AbstractFloat }

Perform a rough alignment of two sets of measured particle coordinates.  The assumption is that the coordinates come from two difference measurements
of the same particle data set which may vary in sample orientation and offset.  The task is to find the affine transformations that
will transform ps1 and ps2 to overlap at the origin.

The particle data sets should include "most" of the same particles or at least a significant overlapping region.  The distance metric should be the 
same in both x & y and also between data sets.

The algorithm uses a nearest neighbor method to identify triplets (or nnsize-lets) of particles.  The triplets are matched from `ps1` with the closest
match in `ps2`.  Then the distances between matched triplets is compared to identify triplets that appear similar in particle spacing and in spacing
between center-of-gravities.  These pairs are then used to compute the rotation angle and the center of gravity for each data set.

  * nnsize is the number of nearest neighbors to consider (nominally 2)
  * tol is the similarity tolerance for separation edge length (<<1.0)
  * corrtol is the correspondance metric multiplier (>1.0)

Example:

    julia> (ct1, ct2) = rough_align(ps1, ps2)
    julia> ct1.(ps1) # Translates ps1 towards the origin (no rotation)
    julia> ct2.(ps2) # Translates ps2 towards the origin and rotates it to overlay `ct1.(ps1)`
    julia> (inv(ct1)∘ct2).(ps2) # Transforms ps2 to overlay ps1
    julia> (inv(ct2)∘ct1).(ps1) # Transforms ps1 to overlay ps2
"""
function rough_align(
    ps1::AbstractVector{<:StaticVector{2, T}}, #
    ps2::AbstractVector{<:StaticVector{2, T}}; #
    nnsize=3, #
    tol=0.001, #
    corrtol=10.0 #
)::NTuple{2, AffineMap} where { T<: AbstractFloat }
    # For each data set find the two nearest neighbors for each particle
    dist1 = collect(zip(knn(KDTree(ps1), ps1, nnsize, true)...))
    dist2 = collect(zip(knn(KDTree(ps2), ps2, nnsize, true)...))
    # Compute the separations between all pairs of nearest neighbors
    function compute_separations(dist, pts)
        function permute2(idx)
            ( ( idx[i], idx[j] ) for i in eachindex(idx) for j in i+1:length(idx) )
        end
        # add a sign to indicate which clockwise/counter-clockwise (additional information to work with!)
        function signme(i1, i2, i, j)
            (i != i1) || (j != i2) ? sign(LinearAlgebra.cross(pts[i1]-pts[i2], pts[i]-pts[j])) : 1.0
        end
        seps= map(dist) do (idx, _)
            SVector{nnsize*(nnsize-1) ÷ 2, Float64}(signme(idx[1], idx[2], i, j) * norm(pts[i] - pts[j]) for (i,j) in permute2(idx) )
        end
        return ( map(d->d[1], dist), seps)
    end
    # Transform the particle data into "nearest neighbor separation space"
    si1 = compute_separations(dist1, ps1)
    si2 = compute_separations(dist2, ps2)
    # Match particle groupings between ps1 and ps2 by matching the ones that are most similar in shape and size
    correspondences = collect(zip(eachindex(si1[2]), nn(KDTree(si2[2]), si1[2])...))
    # Order the correspondences by similarity
    sort!(correspondences, lt=(a,b)->a[3]<b[3])
    # Determine the last index of the "best correspondences" (take at least 10)
    cmin = correspondences[findfirst(c-> c[3] > 0.0, correspondences)][3]
    last = max(min(length(correspondences), 10), findfirst(c->c[3] > corrtol*cmin, correspondences))
    # These functions return the equivalent i-th triplet of particle coordinates in either sp1 or sp2
    triplet1(i::Integer) = ps1[si1[1][correspondences[i][1]]]
    triplet2(i::Integer) = ps2[si2[1][correspondences[i][2]]]
    # Compute the "center-of-gravity" for a series of coordinates
    cog(coords) = mean(coords)
    # Now perform a second filter.  Look at the center-of-gravity for the matching triplets.
    # Build a list of pairs for which the distance between the cog for triplet1(i) and triplet2(i) are similar.
    res = Tuple{Int,Int,Float64}[]
    for ii in 1:last, jj in ii+1:last
        ndc1, ndc2 = norm(cog(triplet1(ii))-cog(triplet1(jj))), norm(cog(triplet2(ii))-cog(triplet2(jj)))
        if abs(ndc1-ndc2) < nnsize*(nnsize-1)*tol && (ndc1 > nnsize*(nnsize-1)*tol)
            push!(res, (ii, jj, abs(ndc1 - ndc2)))
        end
    end
    sort!(res, lt=(a,b)->a[3]<b[3])
    # res now contains the indices into ps1 & ps2 for triplets that match in shape and separation
    rmin = res[findfirst(r->r[3] > 0.0, res)][3]
    lastres = max(min(10,length(res)), findfirst(r->r[3]>10.0*rmin, res))   
    rots = map(res[1:lastres]) do (ii, jj)
        dc1, dc2 = cog(triplet1(ii))-cog(triplet1(jj)), cog(triplet2(ii))-cog(triplet2(jj))
        acos(dot(dc1,dc2)/(norm(dc1)*norm(dc2)))
    end
    # Construct a histogram to determine the most probable angle...
    bins = 0.0:π/180.0:2π
    h = fit(Histogram, rots, bins)
    θi = findmax(h.weights)[2] # most probably rotation angle
    function between_angles(ϕ, ϕmin, ϕmax)
        @assert ϕmin <= ϕmax
        ϕ, ϕmin, ϕmax = mod(ϕ, 2π), mod(ϕmin, 2π), mod(ϕmax, 2π)
        return (ϕmin != ϕmax) && (ϕmin < ϕmax ? (ϕ > ϕmin) && (ϕ < ϕmax) : (ϕ > ϕmin) || (ϕ < ϕmax))
    end    
    # Determine the indices of groupings close to this angle
    idxs = filter(i->between_angles(rots[i], bins[θi]-π/180, bins[θi+1]+π/180), eachindex(rots))
    # Use these indices to compute the best estimates of the rotation and center-of-gravity for ps1 and ps2
    θ = mean(rots[idxs])
    good=unique(mapreduce(x->collect(x[1:2]), append!, res[idxs]))
    cog1 = cog(map(ii->cog(triplet1(ii)), good))
    cog2 = cog(map(ii->cog(triplet2(ii)), good))
    # Convert this into two affine transformations.
    return ( 
        LinearMap(one(RotMatrix{2, Float64}))∘inv(Translation(cog1)),
        LinearMap(RotMatrix{2}(θ))∘inv(Translation(cog2))
    )
end

"""
    correspondences(ps1::AbstractVector{<:StaticVector{2, T}}, ps2::AbstractVector{<:StaticVector{2, T}}; tol=0.01, invert=false) where { T<: AbstractFloat }

Returns a pair of indices `(idx1, idx2)` into `ps1` and `ps2` respectively such that `ps1[idx1]` are likely the same particle as `ps2[idx2]`.

`invert=true` inverts the meaning so returns particles which don't match well with another.
"""
function correspondences(ps1::AbstractVector{<:StaticVector{2, T}}, ps2::AbstractVector{<:StaticVector{2, T}}; tol=0.01, invert=false) where { T<: AbstractFloat }
    (idx1, dist1) = nn(KDTree(ps2), ps1) 
    (idx2, _) = nn(KDTree(ps1), ps2)
    keep = map(eachindex(idx1)) do i
        # Reciprocal matches within tolerance
        b = (i==idx2[idx1[i]]) && (dist1[i] < tol)
        invert ? !b : b
    end
    m1 = idx1[keep]
    m2 = idx2[m1]
    return m2, m1
end

"""
    refined_alignment(ps1::AbstractVector{<:StaticVector{2, T}}, ps2::AbstractVector{<:StaticVector{2, T}}; tol=0.01) where {T <: AbstractFloat}1

Takes two rough aligned data sets `ps1` and `ps2`.  It finds the corresponding particles in each and then uses these
to refine the positions of `ps2` through rotation and translation to best match `ps1`.

Returns a transformed subset of `ps2` that well matches `ps1`.
"""
function refined_alignment(ps1::AbstractVector{<:StaticVector{2, T}}, ps2::AbstractVector{<:StaticVector{2, T}}; tol=0.01) where {T <: AbstractFloat}
    c1, c2 = correspondences2(ps1, ps2; tol=tol, invert=false)
    cps1, cps2 = ps1[c1], ps2[c2]
    com1, com2 = mean(cps1), mean(cps2)
    cpst1, cpst2 = Translation(-com1).(cps1), Translation(-com2).(cps2)
    m2(vs) = transpose(reshape(collect(Iterators.flatten(vs)), length(vs[1]), length(vs)))
    # Solve the orthogonal Procrustes problem
    f = svd(m2(cpst1) * transpose(m2(cpst2)))
    Ω = f.U * f.Vt
    # Now translate and rotated the reordered data towards `pts1`
    return cps1, Translation(com1).(Ω * cpst2)
end
