# Given two particle data sets collected on the same sample but at different
# positions and rotations, find the offset and rotation which best aligns
# datasets.
using LinearAlgebra
using LoopVectorization

"""
    calculate_rings(data, origin, radius, nrings, nθ, sx, sy, dup)

Given a `data` set with `sx`, `sy` columns, compute a matrix which represents
placing each particle within `radius` of `origin` into one of `nrings` equal-area
rings divided into `nθ` angular slices.  The column `:s` is used to weight the
individual particles. `dup` doubles the height of the result matrix and duplicates
the top half below to optimize the shift calculation.
"""
function calculate_rings(data, origin, radius, nrings, nθ, sx, sy, dup)
    res = zeros(Float32, (dup ? 2nθ : nθ, nrings))
    for r in eachrow(data)
        rx, ry = r[sx] - origin[1], r[sy] - origin[2]
        # equal area rings ri ∈ ( 1:nrings )
        ri = floor(Int, nrings * ((ry^2 + rx^2) / radius^2)) + 1
        # @assert ri >= 1 
        if ri <= nrings
            # θi ∈ 1:nθ for θ ∈ 0:2π
            θi = floor(Int, nθ * modf((atan(ry, rx) + 2π) / (2π))[1]) + 1
            # @assert θi >= 1 && θi <= nθ
            res[θi, ri] += r[:s]
        end
    end
    if dup
        res[nθ+1:2nθ, :] = res[1:nθ, :]
    end
    return res
end

function lvdot(a, b)
    s = zero(eltype(a))
    @avx for i ∈ eachindex(a, b)
        s += a[i] * b[i]
    end
    s
end

"""
    calculate_score(data, mask, origin, radius, nrings, nθ)

Compare the `data` against a `mask` calculated using `calculate_rings(...)` to score
the similarity by finding the rotation which gives the maximum score.  The returned
value is a tuple containing the index of the maximum angle range and the associated
score.
"""
function calculate_score(data, mask, origin, radius, nrings, nθ, sx, sy)
    crd = calculate_rings(data, origin, radius, nrings, nθ, sx, sy, false)
    return findmax( [ sum( lvdot(view(crd,:,r), view(mask, off:off+nθ-1, r)) for r in 1:nrings ) for off in 1:nθ ] )
end

"""
`AlignIntermediary` carries the data calculated by `rough_align(...)` associated with 
the scoring of the alignment of two particle data sets.
"""
struct AlignIntermediary
    scores::Array{Float64} # Match scores at the optimal angle
    angles::Array{Int} # The index of the optimal angle
    xsteps::StepRangeLen # Grid steps in the X-dimension
    ysteps::StepRangeLen # Grid steps in the Y-dimension
    nθ::Int # Number of bins into which the range [0, 2π) is divided
end

function Base.show(io::IO, ai::AlignIntermediary) 
    mt = findmax(ai.scores)
    print(io, "AlignIntermediary[sc=$(mt[1]) @ $(coordinate(ai, mt[2])), θ = $(rad2deg(angle(ai, mt[2])))]")
end

"""
    rough_align(data, datap, nrings=8, nθ=16)::AlignIntermediary

Constructs an AlignIntermediary representing the results of comparing all the `data` against `datap`
to attempt to find the best ( offset, rotation ) to match align `data` and `datap`.
"""
function rough_align(data, datap, sx=:x, sy=:y, nrings=8, nθ=32)::AlignIntermediary
    # Find boundaries of datap
    xex, yex, x2ex, y2ex = extrema(data[:,sx]), extrema(data[:,sy]), extrema(datap[:,sx]), extrema(datap[:,sy])
    radius = 0.5 * min(max(x2ex[2] - x2ex[1], y2ex[2] - y2ex[1]), max(xex[2] - xex[1], yex[2] - yex[1]))
    mask = calculate_rings(datap, 0.5 .* (x2ex[2] + x2ex[1], y2ex[2] + y2ex[1]), radius, nrings, nθ, sx, sy, true)
    # Sorting the data means one pass over the data set
    dfx = sort(data, sx)
    xmin, xmax, step = 1, 1, radius / (3.0*nrings)
    # Construct a grid over the data area
    xst, yst = xex[1]:step:xex[2], yex[1]:step:yex[2]
    scores = zeros(Float64, (length(yst), length(xst)))
    angles = zeros(Int, (length(yst), length(xst)))
    for (xi, x) in enumerate(xst)
        xmin = something(
            findnext(i -> dfx[i, sx] >= x - radius, 1:nrow(dfx), xmin),
            nrow(dfx) + 1
        )
        xmax = something( 
                findnext(i -> dfx[i, sx] > x + radius, 1:nrow(dfx), xmax), 
                nrow(dfx) + 1
            ) - 1
        if xmax >= xmin
            # Sorting the data means one pass over the data set
            #@assert issorted(dfx[xmin:xmax, :], sx)
            #@assert xmin == 1 || dfx[xmin - 1, sx] < x - radius "$xmin : $(dfx[xmin - 1, sx]) !< $(x - radius)"
            #@assert xmin == nrow(dfx) + 1 || dfx[xmin, sx] >= x - radius
            #@assert xmax == nrow(dfx) || dfx[xmax, sx] < x + radius
            #@assert xmax == nrow(dfx) || dfx[xmax + 1, sx] > x + radius "$xmax : $(dfx[xmax + 1, sx]) !< $(x + radius)"
            dfy, ymin, ymax = sort(dfx[xmin:xmax, :], sy), 1, 1
            for (yi, y) in enumerate(yst)
                ymin = something( 
                        findnext(i -> dfy[i, sy] >= y - radius, 1:nrow(dfy), ymin),
                        nrow(dfy) + 1
                    )
                ymax = something( 
                        findnext(i -> dfy[i, sy] > y + radius, 1:nrow(dfy), ymax), 
                        nrow(dfy) + 1
                    ) - 1
                if ymax >= ymin
                    #@assert issorted(dfy[ymin:ymax, :], sy)
                    #@assert ymin == 1 || dfy[ymin - 1, sy] < y - radius "$ymin : $(dfy[ymin - 1, sy]) !< $(y - radius)"
                    #@assert dfy[ymin, sy] >= y - radius
                    #@assert dfy[ymax, sy] < y + radius
                    #@assert ymax == nrow(dfy) || dfy[ymax + 1, sy] > y + radius "$ymax : $(dfy[ymax + 1, sy]) !< $(y + radius)"
                    (scores[yi, xi], angles[yi, xi]) = calculate_score(dfy[ymin:ymax, :], mask, (x, y), radius, nrings, nθ, sx, sy)
                end
            end
        end
    end
    return AlignIntermediary(scores, angles, xst, yst, nθ)
end

"""
    rough_align(data, datap, nrings=8, nθ=16)::AlignIntermediate

Constructs an AlignIntermediary representing the results of comparing all the `data` against `datap`
to attempt to find the best ( offset, rotation ) to match align `data` and `datap`.
"""
function rough_align_slow(data, datap, sx=:x, sy=:y, nrings=8, nθ=16)::AlignIntermediary
    # Find boundaries of datap
    x2ex, y2ex = extrema(datap[:,sx]), extrema(datap[:,sy])
    radius = 0.5 * max(x2ex[2] - x2ex[1], y2ex[2] - y2ex[1])
    mask = calculate_rings(datap, 0.5 .* (x2ex[2] + x2ex[1], y2ex[2] + y2ex[1]), radius, nrings, nθ, sx, sy, true)
    # Sorting the data means one pass over the data set
    step = 0.5 * radius / nrings
    xex, yex = extrema(data[:,sx]), extrema(data[:,sy])
    # Construct a grid over the data area
    xst, yst = xex[1]:step:xex[2], yex[1]:step:yex[2]
    scores = zeros(Float64, (length(yst), length(xst)))
    angles = zeros(Int, (length(yst), length(xst)))
    for (xi, x) in enumerate(xst), (yi, y) in enumerate(yst)
        (scores[yi, xi], angles[yi, xi]) = calculate_score(data, mask, (x, y), radius, nrings, nθ, sx, sy)
    end
    return AlignIntermediary(scores, angles, xst, yst, nθ)
end

Base.CartesianIndices(alint::AlignIntermediary) = CartesianIndices(alint.scores)

"""
    offset(alint::AlignIntermediary, ma::CartesianIndex)

Returns the offset at the specified CartesianIndex within `alint`. 
"""
offset(alint::AlignIntermediary, ma::CartesianIndex) = (alint.xsteps[ma[2]], alint.ysteps[ma[1]])

"""
    angle(alint::AlignIntermediary, ma::CartesianIndex)

Returns the center at the specified CartesianIndex within `alint`. 
"""
Base.angle(alint::AlignIntermediary, ma::CartesianIndex) = 2π*modf((alint.angles[ma]-1.0) / alint.nθ + 3.0)[1]

"""
    score(alint::AlignIntermediary, ma::CartesianIndex)

Returns the score at index `ma` 
"""
score(alint::AlignIntermediary, ma::CartesianIndex) = alint.scores[ma]


"""
    findbest(ai::AlignIntermediary, f=0.8)

Returns an array of indices into `ai` for which the score is greater than `f` times the maximum score.
"""
function findbest(ai::AlignIntermediary, f=0.8)
    tmp = sort(reshape(collect(CartesianIndices(ai.scores)), prod(size(ai.scores))), by=i -> ai.scores[i], rev=true)
    i = something(findfirst(i -> ai.scores[tmp[i]] < f * ai.scores[tmp[1]], eachindex(tmp)), length(tmp)+1)
    return tmp[1:i - 1]
end


function align(datap, ai::AlignIntermediary, idx=1, xs=:x, ys=:y)
    res = copy(datap)
    bst = findbest(ai,0.8)[idx]
    θ, off = -angle(ai, bst), offset(ai, bst)
    rot = [ cos(θ) -sin(θ); sin(θ) cos(θ) ]
    xps, yps = Float64[], Float64[]
    for row in eachrow(datap)
        (xp, yp) = (rot * [ row[xs], row[ys] ]) .+ off
        push!(xps, xp)
        push!(yps, yp)
    end
    res[:, xs] = xps
    res[:, ys] = yps
    return res
end
