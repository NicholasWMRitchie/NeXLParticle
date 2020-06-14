using LinearAlgebra: dot
using OnlineStatsBase

function rca(img::AbstractArray, start::CartesianIndex, thresh::Function, nchords = 16)
    function measurechord(ci, slope)
        dx, sx, c = abs(slope[2]), slope[2] > 0 ? 1 : -1, ci.I[2]
        dy, sy, r = -abs(slope[1]), slope[1] > 0 ? 1 : -1, ci.I[1]
        err = dx + dy
        last = (r, c)
        while c >= 1 && c <= size(img, 2) && r >= 1 && r <= size(img, 1) && thresh(img[r, c])
            last = (r, c)
            e2 = 2err
            if e2 >= dy
                err += dy
                c += sx
            end
            if e2 <= dx
                err += dx
                r += sy
            end
        end
        return CartesianIndex(last...)
    end
    @assert thresh(img[start]) "The initial point does not meet the threshold in drawchord(...)"
    bisectors = (((1, 0), (-1, 0)), ((0, 1), (0, -1)), ((1, 1), (-1, -1)), ((-1, 1), (1, -1)))
    rat = map(i -> rationalize(Int32, tan(2π * i / (2 * nchords)), tol = 1.0e-8), 1:nchords÷2) # 9.3 μs
    chords = (
        ((-numerator(r), -denominator(r)) for r in rat)..., #
        ((-denominator(r), numerator(r)) for r in rat)..., #
        ((numerator(r), denominator(r)) for r in rat)..., #
        ((denominator(r), -numerator(r)) for r in rat)..., #
    )
    ci = start # Use successive approximation to find the center
    for (pd, nd) in bisectors
        ci = CartesianIndex((measurechord(ci, pd).I .+ measurechord(ci, nd).I) .÷ 2)
    end
    return CartesianIndex[measurechord(ci, s) for s in chords]
end

"""
    rca(img::AbstractArray, search::Function, measure::Function, nchords = 16)
    rca(blob::Vector{Blob})

Search an image for pixels that meet the `search` threshold.  Then using the the rotation chord algorithm
find the approximate center of the particle and draw 2 `nchords` at equally spaced angles from the center
until the `measure` threshold is no longer met.  The resulting set of perimeter points is stored
as a `Vector{CartesianIndex}`.  It is possible that many regios will meet the `search` threshold so the
actual result is a `Vector{Vector{CartesianIndex}}`.

Note: The measure `threshold` should be more permissive than the `search` threshold.

Example

    rca(img, p->p>0.9, p->p>0.7) #
    bs = blob(img, p->p>0.7)
    rca(bs[1])
"""
function rca(img::AbstractArray, search::Function, measure::Function, nchords = 16)
    function bnds(vci::Vector{CartesianIndex})
        l, r, t, b = size(img, 2), 1, size(img, 1), 1
        for ci in vci
            l = ci.I[2] < l ? ci.I[2] : l
            r = ci.I[2] > r ? ci.I[2] : r
            t = ci.I[1] < t ? ci.I[1] : t
            b = ci.I[1] > b ? ci.I[1] : b
        end
        return CartesianIndices((t:b, l:r))
    end
    res, bounds = Vector{Vector{CartesianIndex}}(), Vector{CartesianIndices}()
    skip = false  # Fast skip over columns that have already been checked...
    for ci in CartesianIndices(img)
        if search(img[ci])
            if !(skip || any(map(cis -> ci in cis, bounds)))
                @assert measure(img[ci]) "The 'measure' threshold must be designed such that a pixel that meets the 'search' threshold must also meet the 'measure' threshold."
                r = rca(img, ci, measure, nchords)
                push!(res, r)
                push!(bounds, bnds(r))
                skip = true
            end
        else
            skip = false
        end
    end
    return res
end

function rca(bs::Vector{Blob})
    function offrca(b)  # RCA offset to the original image coordinates
        offci(bbi, r) = CartesianIndex[ CartesianIndex( map( i->ci.I[i] + bbi[i].start - 1, eachindex(ci.I))...) for ci in r ]
        return [ offci(b.bounds.indices, r) for r in rca(b.mask, p->p, p->p) ]
    end
    return mapreduce(b->offrca(b), append!, bs)
end

function colorizedimage(chords::Vector{Vector{CartesianIndex}}, img::AbstractArray)
    ends = RGB(1.0,0.0,0.0)
    function drawchords(res::Array, chords::Vector{CartesianIndex}, cl::Color)
        function drawchord(res, ci1, ci2)
            dx, sx, xa = abs(ci2.I[2] - ci1.I[2]), ci1.I[2] < ci2.I[2] ? 1 : -1, ci1.I[2]
            dy, sy, ya = -abs(ci2.I[1] - ci1.I[1]), ci1.I[1] < ci2.I[1] ? 1 : -1, ci1.I[1]
            err = dx + dy
            while !((xa == ci2.I[2]) && (ya == ci2.I[1]))
                e2 = 2err
                if e2 >= dy
                    err += dy
                    xa += sx
                end
                if e2 <= dx
                    err += dx
                    ya += sy
                end
                res[ya, xa] = cl
            end
            res[ci1]=ends
            res[ci2]=ends
        end
        lo2 = length(chords) ÷ 2
        for i = 1:lo2
            drawchord(res, chords[i], chords[i+lo2])
        end
        return res
    end
    colors =
        convert.(
            RGB,
            distinguishable_colors(
                length(chords) + 2,
                Color[RGB(253 / 255, 255 / 255, 255 / 255), RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
                transform = deuteranopic,
            )[3:end],
        )
    res = RGB.(img)
    foreach(z -> drawchords(res, z[1], z[2]), zip(chords, colors))
    return res
end

"""
    metrics(chords::Vector{CartesianIndex})::Dict{Symbol,Float64}

Take a `Vector{CartesianIndex}` from `rca(...)` and compute morphology metrics.
"""
function metrics(chords::Vector{CartesianIndex}, scale=1.0)::Dict{Symbol,Float64}
    lo2, lo4 = length(chords) ÷ 2, length(chords) ÷ 4
    len(a) = sqrt(dot(a, a))
    leni(i) = len(chords[i].I .- chords[i+lo2].I)
    triarea(a, b) = 0.5 * abs(a[1] * b[2] - a[2] * b[1]) # area of triangle with edges a and b
    s, extrema = Moments(), Extrema()
    dmax, imax = 0.0, -1
    for i in 1:lo2
        d = leni(i)
        if d>dmax
            dmax=d
            imax=i
        end
        fit!(extrema, d)
        fit!(s, d)
    end
    perim = len(chords[1].I .- chords[end].I)
    center = (chords[1].I[1], chords[lo4+1].I[2])
    area = triarea(chords[1].I .- center, chords[end].I .- center)
    yext, xext = Extrema(), Extrema()
    fit!(yext, chords[1].I[1])
    fit!(xext, chords[1].I[2])
    for i in 2:length(chords)
        perim += len(chords[i].I .- chords[i-1].I)
        area += triarea(chords[i].I .- center, chords[i-1].I .- center)
        fit!(yext, chords[i].I[1])
        fit!(xext, chords[i].I[2])
    end
    dperp = leni((imax -1 + lo4) % lo2 + 1)
    return Dict(
        :DMIN => scale * minimum(extrema),
        :DMAX => scale * maximum(extrema),
        :DAVG => scale * Statistics.mean(s),
        :DSIG => scale * Statistics.var(s),
        :DPERP => dperp,
        :ASPECT => dmax / dperp,
        :PERIM => scale * perim,
        :XCENT => center[2],
        :YCENT => center[1],
        :ORIENT2 => 2π * imax / length(chords) - π,
        :XFERET => scale * (maximum(xext) - minimum(xext)),
        :YFERET => scale * (maximum(yext) - minimum(yext)),
        :AREA => scale^2 * area,
    )
end

function area(chords::Vector{CartesianIndex}, scale=1.0)
    triarea(a, b) = 0.5 * abs(a[1] * b[2] - a[2] * b[1]) # area of triangle with edges a and b
    center = (chords[1].I[1], chords[length(chords)÷4+1].I[2])
    area = triarea(chords[1].I .- center, chords[end].I .- center)
    for i in 2:length(chords)
        area += triarea(chords[i].I .- center, chords[i-1].I .- center)
    end
    return area
end
