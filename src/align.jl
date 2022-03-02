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
using LsqFit
using InvertedIndices: Not

"""
generate_ground_truth(rnd::Random, bounds::Rect2, npart::Integer, inside::Function=p->p in bounds)

Generate a simulated ground-truth particle data set.  This represents the "true" location of `npart` particles that are randomly 
distributed within `bounds` and for which `inside(pt)` evaluates true.
"""
function generate_ground_truth(rnd::AbstractRNG, bounds::Rect2, npart::Integer, inside::Function = p -> p in bounds)
    @assert minimum(bounds.widths) > 0.0 "Rectangle widths must be positive."
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
function measure_particles(rnd::AbstractRNG, pd::AbstractVector{<:SVector{2}}, offset::SVector{2}, θ::Float64, frac::Float64, eps::Float64)
    @assert frac > 0.0 "The fraction measured `frac` must be larger than 0.0."
    @assert frac <= 1.0 "The fraction measured `frac` must be less than or equal to 1.0."
    rot = RotMatrix{2}(θ)
    n = Distributions.Normal(0.0, eps)
    return shuffle!(rnd, map(filter(f -> rand(rnd) < frac, pd)) do part
        rot * part + offset + rand(rnd, n, 2)
    end)
end

"""
    rough_align(ps1::AbstractVector{StaticVector{2, T}}, ps2::AbstractVector{StaticVector{2, T}}; groupsize=3, tol=0.001, corrtol=10.0) where { T<: AbstractFloat }

Perform a rough alignment of two sets of measured particle coordinates.  The assumption is that the coordinates come from two difference measurements
of the same particle data set which may vary in sample orientation and offset.  The task is to find the affine transformations that
will transform ps1 and ps2 to overlap at the origin.

The particle data sets should include "most" of the same particles or at least a significant overlapping region.  The distance metric should be the 
same in both x & y and also between data sets.

The algorithm uses a nearest neighbor method to identify triplets (or groupsize-lets) of particles.  The triplets are matched from `ps1` with the closest
match in `ps2`.  Then the distances between matched triplets is compared to identify triplets that appear similar in particle spacing and in spacing
between center-of-gravities.  These pairs are then used to compute the rotation angle and the center of gravity for each data set.

  * `groupsize` is the number of nearest neighbors to consider (nominally 3)
  * `tol` is the similarity tolerance for separation edge length (<<1.0). `tol` should be approximately the uncertainty in each
  component of the measured positions.

A description of the algorithm is available [here](https://docs.google.com/presentation/d/1XvkeB1AY75s0X-f992S1Qf2H895_9MoI/edit?usp=sharing&ouid=113647574378052266554&rtpof=true&sd=true).

Example:

    julia> (ct1, ct2) = rough_align(ps1, ps2)
    julia> ct1.(ps1) # Translates ps1 towards the origin (no rotation)
    julia> ct2.(ps2) # Translates ps2 towards the origin and rotates it to overlay `ct1.(ps1)`
    julia> (inv(ct1)∘ct2).(ps2) # Transforms ps2 to overlay ps1
    julia> (inv(ct2)∘ct1).(ps1) # Transforms ps1 to overlay ps2
"""
function rough_align(
    ps1::AbstractVector{<:StaticVector{2,T}}, #
    ps2::AbstractVector{<:StaticVector{2,T}}; #
    groupsize = 3, #
    tol = 0.001 #
)::NTuple{2,AffineMap} where {T<:AbstractFloat}
    # For each data set find the (groupsize-1) nearest neighbors for each particle
    idx1, _ = knn(KDTree(ps1), ps1, groupsize, true)
    idx2, _ = knn(KDTree(ps2), ps2, groupsize, true)
    # Compute the separations between all pairs of nearest neighbors
    function compute_separations(idxs, pts)
        # add a sign to indicate clockwise vs counter-clockwise
        signednorm(a, b)::T = (a[1] * b[2] > a[2] * b[1] ? -one(T) : one(T)) * norm(b)
        map(idxs) do idx
            # DOF = groupsize particles times 2 dimensions - (rotate, x_off, y_off)
            dnn12 = pts[idx[1]] - pts[idx[2]]
            SA[collect(signednorm(dnn12, pts[idx[i]] - pts[idx[j]]) for i in Base.OneTo(groupsize) for j in i+1:groupsize)...]
        end # separations 
    end
    # Transform the particle data into "nearest neighbor separation space"
    si1 = compute_separations(idx1, ps1)
    si2 = compute_separations(idx2, ps2)
    # Match particle groupings between ps1 and ps2 by matching the ones that are most similar in shape and size
    correspondences = collect(zip(eachindex(idx1), nn(KDTree(si2), si1)...)) # index in idx1, index in idx2, distance
    # Order the correspondences by similarity
    sort!(correspondences, lt = (a, b) -> a[3] < b[3])
    # These functions return the equivalent i-th triplet of particle coordinates in either sp1 or sp2
    triplet1(i::Integer) = ps1[idx1[correspondences[i][1]]]
    triplet2(i::Integer) = ps2[idx2[correspondences[i][2]]]
    # Compute the "center-of-gravity" for a series of coordinates
    cog(coords) = mean(coords)
    # Now perform a second filter.  Look at the center-of-gravity for the matching triplets.
    # Build a list of pairs for which the distance between the cog for triplet1(i) and triplet2(i) are similar.
    res = Tuple{Int,Int,Float64,Float64}[] # index into tripletX(), index into tripletX(), difference in separation, angle between
    clast = min(length(correspondences), 100)
    for ii in 1:clast, jj in ii+1:clast
        dc1, dc2 = cog(triplet1(ii)) - cog(triplet1(jj)), cog(triplet2(ii)) - cog(triplet2(jj))
        ndc1, ndc2 = norm(dc1), norm(dc2)
        # Use determinant to diffentiate 0 to π from π to 2π
        dets = dc1[1] * dc2[2] > dc1[2] * dc2[1] ? -one(T) : one(T)
        # if the norms are exactly equal to zero then they likely represent the same triplets but starting with a different seed particle
        if abs(ndc1 - ndc2) < 6 * tol && ((ndc1 > 1.0e-12) || (ndc2 > 1.0e-12))
            push!(res, (ii, jj, abs(ndc1 - ndc2), mod(dets * acos(clamp(dc1 ⋅ dc2 / (ndc1 * ndc2), -1.0, 1.0)), 2π)))
        end
    end
    # res now contains the indices into ps1 & ps2 for triplets that match in shape and separation
    sort!(res, lt = (a, b) -> a[3] < b[3]) # sort by difference in length
    # Handles cyclic boundaries on angle measurement
    function between_angles(ϕ, ϕmin, ϕmax)
        @assert ϕmin <= ϕmax
        ϕ, ϕmin, ϕmax = mod(ϕ, 2π), mod(ϕmin, 2π), mod(ϕmax, 2π)
        return (ϕmin != ϕmax) && (ϕmin < ϕmax ? (ϕ > ϕmin) && (ϕ < ϕmax) : (ϕ > ϕmin) || (ϕ < ϕmax))
    end
    # Construct a histogram to determine the most probable angle...
    bins = 0.0:π/180.0:2π
    h = fit(Histogram, map(r -> r[4], res), bins)
    θi = findmax(h.weights)[2] # most probably rotation angle
    # Find the indices of corresponding pairs of groups close to this angle
    idxs = filter(i -> between_angles(res[i][4], bins[θi] - π / 180, bins[θi+1] + π / 180), eachindex(res))
    @assert length(idxs) > 0 "θi = $(findmax(h.weights))"
    good = res[idxs]
    # Compute the best estimate of the rotation angle
    θ = mean(r -> r[4], good)
    goodidxs = unique(mapreduce(g -> [g[1], g[2]], append!, good)) # g[1] & g[2] are indices of corresponding particle groups
    # Use these indices to compute the best estimates of the center-of-gravity for ps1 and ps2
    cog1 = cog(map(i -> cog(triplet1(i)), goodidxs))
    cog2 = cog(map(i -> cog(triplet2(i)), goodidxs))
    # Convert this into two affine transformations.
    return (
        LinearMap(Angle2d(zero(T))) ∘ inv(Translation(cog1)),
        LinearMap(Angle2d(T(θ))) ∘ inv(Translation(cog2))
    )
end

"""
    correspondences(ps1::AbstractVector{<:StaticVector{2, T}}, ps2::AbstractVector{<:StaticVector{2, T}}; tol=0.01, invert=false) where { T<: AbstractFloat }

Returns a pair of indices `(idx1, idx2)` into `ps1` and `ps2` respectively such that `ps1[idx1]` are likely the same particle as `ps2[idx2]`.

`invert=true` inverts the meaning so returns particles which don't match well with another.
"""
function correspondences(ps1::AbstractVector{<:StaticVector{2,T}}, ps2::AbstractVector{<:StaticVector{2,T}}; tol = 0.01, invert = false) where {T<:AbstractFloat}
    (idx1, _) = nn(KDTree(ps2), ps1) # idx1 is an index into ps2
    (idx2, dist2) = nn(KDTree(ps1), ps2) # idx2 is an index into ps1
    keep2 = map(enumerate(idx2)) do (i2, i1)
        # Within tolerance and closest alternative in either direction
        (dist2[i2] < tol) && (i2==idx1[i1])
    end
    c1 = idx2[keep2]
    c2 = idx1[c1]
    return invert ? ( Not(c1), Not(c2) ) : ( c1, c2 )
end

function orthogonal_procrustes_alignment(ps1::AbstractVector{<:StaticVector{2,T}}, ps2::AbstractVector{<:StaticVector{2,T}}; tol = 0.01) where {T<:AbstractFloat}
    c1, c2 = correspondences(ps1, ps2; tol = tol, invert = false)
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


"""
    refined_alignment(ps1::Vector{<:StaticVector{2,T}}, ps2::Vector{<:StaticVector{2,T}}) where { T <: AbstractFloat }

Uses a non-linear algorithm to refine a rough alignment by adjusting the rotation and offset to produce the smallest
mean-square distance between particles in `ps1` and `ps2`.  Return the `AffineMap` that best transforms `ps2` to match
`ps1`.

The inputs `ps1` and `ps2` must correspond index-by-index with each other. Meaning they represent the exact same particles
in the same order.  This is where the `correspondences(...)` function comes in handy.

Example:
    
    julia> ct1, ct2 = rough_align(ps1, ps2)                  # Get the alignment transforms
    julia  rs1, rs2 = ct1.(ps1),ct2.(ps2)                    # Rough align the data sets
    julia> ci1, ci2 = correspondences(rs1, rs2)              # Identify corresponding rough-aligned particles by index
    julia> cs1, cs2 = rs1[ci1], rs2[ci2]                     # Extract these particles
    julia> n(ds1, ds2) = sum(x->x^2, norm.(ds1 .- ds2))       # Metric function
    julia> n(rs1[ci1], rs2[ci2])                             # Pre-Error measure
    julia> ct3 = refined_alignment(cs1, cs2)                 # Get the refined alignment transform
    julia> n(rs1, ct3.(rs2))                                 # Post-Error measure
    julia> n(rs1, (ct3∘ct2).(cs2));
    julia> aps1, aps2 = ct1.(ps1), (ct3∘ct2).(ps2)           # Transform all points
"""
function refined_alignment(ps1::AbstractVector{<:StaticVector{2,T}}, ps2::AbstractVector{<:StaticVector{2,T}}) where {T<:AbstractFloat}
    @assert length(ps1) == length(ps2)
    function flatten(ps)
        res = Array{T}(undef, 2 * length(ps))
        foreach(i -> res[2i-1:2i] .= ps[i], eachindex(ps))
        return res
    end
    # Compute the function
    function f(res, ps, param)
        A, B, C, D = cos(param[1]), sin(param[1]), param[2], param[3]
        for i in 1:2:length(ps)
            x, y = ps[i], ps[i+1]
            res[i] = A * (x + C) - B * (y + D)
            res[i+1] = B * (x + C) + A * (y + D)
        end
    end
    # Compute the Jacobian
    function jac(res, ps, param)
        A, B, C, D = cos(param[1]), sin(param[1]), param[2], param[3]
        for i in 1:2:length(ps)
            x, y = ps[i], ps[i+1]
            res[i, :] .= (-B * (x + C) - A * (y + D), A, -B)
            res[i+1, :] .= (-B * (y + D) + (x + C) * A, B, A)
        end
        return res
    end
    fit = curve_fit(f, jac, flatten(ps2), flatten(ps1), [0.0, 0.0, 0.0]; inplace = true)
    return LinearMap(Angle2d(fit.param[1])) ∘ Translation(fit.param[2:3])
end


"""
    align(ps1::Vector{<:StaticVector{2,T}}, ps2::Vector{<:StaticVector{2,T}}; tol=0.001, finealign=true)  where { T <: AbstractFloat }

Perform a `rough_align(...)` followed by a "refined_alignment(...)" on `ps1` and `ps2`.  Return the `AffineMap`
transformations `(ct1, ct2) such that `ct1.(ps1)` and `ct2.(ps2)` are registered.  If `finealign=true` then a 
second-stage non-linear optimization of corresponding points is performed.

Example:

    julia> (ct1, ct2) = align(ps1, ps2)
    julia> ct1.(ps1) # Translates ps1 towards the origin (no rotation)
    julia> ct2.(ps2) # Translates ps2 towards the origin and rotates it to overlay `ct1.(ps1)`
    julia> (inv(ct1)∘ct2).(ps2) # Transforms ps2 to overlay ps1
    julia> (inv(ct2)∘ct1).(ps1) # Transforms ps1 to overlay ps2
"""
function align(ps1::AbstractVector{<:StaticVector{2,T}}, ps2::AbstractVector{<:StaticVector{2,T}}; tol = 0.001, finealign = true) where {T<:AbstractFloat}
    ct1, ct2 = rough_align(ps1, ps2; tol = tol)
    if finealign
        rs1, rs2 = ct1.(ps1), ct2.(ps2)
        ci1, ci2 = correspondences(rs1, rs2)
        cs1, cs2 = rs1[ci1], rs2[ci2]
        ct3 = refined_alignment(cs1, cs2)
        ct2 = ct3 ∘ ct2
    end
    return (ct1, ct2)
end

"""
    identify(pss::AbstractVector{<:AbstractVector{<:StaticVector{2,T}}}; tol=0.001, ctol=0.01 )::DataFrame where { T <: AbstractFloat}

Takes a vector of X-Y coordinate vectors.  Aligns the X-Y coordinate vectors to the first one.  Then it labels each unique particle and
identifies the index (if any) in which it appears in each data set.  Thus we are able to track particle identities through the data sets.
Each row in the returned `DataFrame` represents a "unique particle".  The `:PS?` columns represent the index of the particle in that
X-Y coordinate vectors. The `:COUNT` column represents the number of data sets in which the particle was measured.  The `:FIRST` column
indicates the first data set in which the particle was observed.
"""
function identify(pss::AbstractVector{<:AbstractVector{<:StaticVector{2,T}}}; tol=0.001, ctol=0.01 )::DataFrame where { T <: AbstractFloat}
    # Align all pss relative to pss[1]
    als = map(pss) do psj
        if psj==pss[1]
            pss[1]
        else
            cti, ctj = align(pss[1], psj, tol=tol, finealign=true)
            (inv(cti)∘ctj).(psj)
        end
    end
    df = DataFrame( (Symbol("PS$i")=> Union{Missing,Int}[] for i in eachindex(als))...)
    for i in eachindex(als)
        alsi = als[i]
        # Determine the indices of the previously unseen particles in als[i]
        unseen = fill(true, length(alsi))
        foreach(i-> (!ismissing(i)) && (unseen[i]=false), df[:,i])
        idxs=collect(eachindex(alsi)[unseen])
        # Get the coordinates of these particle
        uni=collect(alsi[idxs])
        # Find correspondences for these particles with later data sets
        cijs = map(i+1:length(als)) do j
            ci, cj = correspondences(uni, als[j], tol=ctol)
            tmp = zeros(Int, length(uni))
            tmp[ci] .= cj
            tmp
        end
        ztom(x) = x==zero(typeof(x)) ? missing : x
        for (ii, ind) in enumerate(idxs)
            # uni[ii]=alsi[ind]
            push!(df, [ ( missing for jj in 1:i-1)..., ind, ( ztom(cij[ii]) for cij in cijs )...])
        end
    end
    f = CategoricalArray(map(eachrow(df)) do r
        findfirst(c->!ismissing(c), Tuple(r))
    end, ordered=true)
    c = CategoricalArray(map(eachrow(df)) do r
        count(c->!ismissing(c), Tuple(r))
    end, ordered=true)
    x = map(eachrow(df)) do r
        mean(skipmissing(map(i->ismissing(r[i]) ? missing : als[i][r[i]][1], 1:length(r))))
    end
    y = map(eachrow(df)) do r
        mean(skipmissing(map(i->ismissing(r[i]) ? missing : als[i][r[i]][2], 1:length(r))))
    end
    insertcols!(df, :APPEARS => f, :COUNT => c, :XABS => x, :YABS => y )
    sort!(df, [ :APPEARS, order(:COUNT,rev=true) ], alg=MergeSort) # Maintains order
    df
end