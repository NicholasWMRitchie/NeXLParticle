using FileIO
using Images
using LinearAlgebra
using StatsBase


"""
A Blob is a mask consisting of blocks of adjacent pixels meeting a threshold.
Only those pixels immediately above, below, left or right are considered to
be adjacent (not diagonals).  Blobs offer easy access to a mask consisting
of pixels in the blob and the perimeter either as a list of steps from
one perimeter point to the next or a list of CartesianIndex with the coordinates
of the perimeter points.
"""
struct Blob
    bounds::CartesianIndices
    mask::BitArray  # Mask of those pixels in the blob
    pstart::CartesianIndex # Start of the perimeter
    psteps::Vector{Tuple{Int,Int}} # Steps around perimeter

    function Blob(bounds::CartesianIndices, mask::BitArray)
        function perimeter(m)  # Computes the perimeter steps
            stps = Tuple{Int,Int}[]
            if prod(size(m))>1
                msk(r,c) = checkbounds(Bool, m, r, c) && m[r,c]
                mod8(x) = (x + 7) % 8 + 1 # maintains 1...8
                next(cur, drs) = mod8(drs[findfirst(r->msk((cur .+ steps[mod8(r)])...), drs)])
                steps = ( (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1),  )
                pre = ( findfirst(r->m[r,1],1:size(m,1)), 1 )
                @assert msk(pre...)
                ff = findfirst(r->msk((pre .+ steps[mod8(r)])...), 1:5)
                if isnothing(ff)
                    FileIO.save(File(format"PNG", "dump.png"),mask)
                    @assert !isnothing(ff) "bounds=>$bounds, m[pre]=>$(m[pre...]) m[b]=>$([msk((pre .+ s)...) for s in steps])"
                end
                prevdir = next(pre, 1:5)
                start=pre .+ steps[prevdir]
                @assert msk(start...)
                curr = start
                while true
                    nextdir = next(curr, prevdir-2:prevdir+4)
                    push!(stps, steps[nextdir])
                    curr = curr .+ steps[nextdir]
                    @assert msk(curr...)
                    prevdir = nextdir
                    if curr == start && steps[next(curr, prevdir-2:prevdir+4)] == stps[1]
                        break # end at the start going the same direction
                    end
                end
                return ( CartesianIndex(start), stps)
            else
                return ( CartesianIndex(1,1), stps)
            end
        end
        @assert ndims(mask)==2
        @assert size(bounds)==size(mask)
        ci, stps = perimeter(mask)
        # cip = CartesianIndex([ci.I[i]+bounds.indices[i].start-1 for i in eachindex(ci.I)]...)
        return new(bounds, mask, ci, stps)
    end
end


"""
    blob(img::AbstractArray, thresh::Function)::Vector{Blob}

Create a Vector of `Blob`s containing discontiguous regions meeting the threshold function.  The
result is sorted by `Blob` area.
"""
function blob(img::AbstractArray, thresh::Function)::Vector{Blob}
    function extractmask(res, i, minR, maxR, minC, maxC)
        mask = BitArray(undef, maxR - minR + 1, maxC - minC + 1)
        foreach(ci->mask[ci] = (res[(ci.I .+ (minR-1, minC-1))...] == i), CartesianIndices(mask))
        return mask
    end
    res = zeros(UInt16, size(img))
    alias = Vector{Set{eltype(res)}}()
    prev, next = zero(eltype(res)), zero(eltype(res))
    for ci in CartesianIndices(img)
        ci.I[1] == 1 && (prev = 0)==0  # Reset prev each column
        if thresh(img[ci])
            above = ci[2] > 1 ? res[(ci.I .- (0,1))...] : zero(eltype(res))
            if above ≠ 0
                if (prev ≠ 0) && (above ≠ prev)
                    # Same blob => different indexes
                    ia = findfirst(al -> above in al, alias)
                    ip = findfirst(al -> prev in al, alias)
                    @assert !(isnothing(ia) || isnothing(ip)) "$prev or $above not in $alias"
                    if ia ≠ ip # merge regions
                        union!(alias[ia], alias[ip])
                        deleteat!(alias, ip)
                    end
                end
                prev = above
            elseif prev == 0
                prev = (next += one(next))
                push!(alias, Set{UInt16}(prev)) # new region
            end
            @assert prev ≠ 0
            res[ci] = prev
        else
            prev = zero(eltype(res))
        end
    end
    # Second pass combine the adjacent indices
    newidx, nblobs = zeros(eltype(res), next), length(alias)
    for i in 1:nblobs, id in alias[i]
        newidx[id] = i
    end
    rects = Dict{eltype(res),NTuple{4,Int}}(i => (100000, -1, 100000, -1) for i in 1:nblobs)
    for ci in CartesianIndices(img)
        if (bidx = res[ci]) ≠ zero(eltype(res)) # belongs to a blob
            ni=(res[ci] = newidx[bidx])
            rect = rects[ni]
            rects[ni] = (min(ci[1], rect[1]), max(ci[1], rect[2]), min(ci[2], rect[3]), max(ci[2], rect[4]))
        end
    end
    blobs = [Blob(CartesianIndices((rect[1]:rect[2], rect[3]:rect[4])), extractmask(res, i, rect...)) for (i, rect) in rects]
    sort!(blobs, lt=(b1,b2) -> area(b1) > area(b2))
    return blobs
end

"""
    Base.CartesianIndices(b::Blob)

Bounds of the blob in the original image's coordinate system
"""
Base.CartesianIndices(b::Blob) = b.bounds

"""
    Base.getindex(b::Blob, ci::CartesianIndex)

Whether a pixel in the original image's coordinate system is in the blob.
"""
Base.getindex(b::Blob, ci::CartesianIndex) = #
    (ci in b.bounds) && b.mask[map(i->ci.I[i] - b.bounds.indices[i].start + 1, eachindex(ci.I))...]

"""
    perimeter(b::Blob)::Vector{CartesianIndex}

Returns a vector of `CartesianIndex` corresponding to the points around the
perimeter of the blob in the original image's coordinate system.
"""
function perimeter(b::Blob)::Vector{CartesianIndex}
    pts, acc = CartesianIndex[ b.pstart ], [ b.pstart.I... ]
    foreach(stp->push!(pts, CartesianIndex((acc .+= stp)...)), b.psteps[1:end-1])
    return pts
end


"""
    perimeterlength(b::Blob)

Compute the length of the blob perimeter.  Diagonals are √2 and straights are 1.
"""
perimeterlength(b::Blob) =
    mapreduce(st->sqrt(dot(st,st)), +, b.psteps)

"""
    ecd(b::Blob, filled = true)

Computes the equivalent circular diameter (by default computes the ecd(..) based on the area including the area
of any interior holes.)  A = πr² => ecd = 2r = 2√(A/π)
"""
ecd(b::Blob, filled=true) = 2.0 * sqrt((filled ? filledarea(b) : area(b)) / π)

"""
    curvature(b::Blob, n::Int)

Compute an array that measures the angular difference between the pixel n before
and n after then current one for each point on the perimeter of the blob.
Negative indicates convex and positive indicates concave with large positive
values indicating a sharp concave angle.
"""
function curvature(b::Blob, n::Int)
    modn(i) = (i+length(b.psteps)-1) % length(b.psteps) + 1
    stepsum(itr) = mapreduce(j->b.psteps[modn(j)], (x,y)->.+(x,y), itr, init=(0,0))
    angles = Float64[]
    for i in eachindex(b.psteps)
        sm, sp = -1 .* stepsum(i-1:-1:i-n), stepsum(i:i+n-1)
        den² = dot(sm,sm)*dot(sp,sp)
        ac = den² > 0 ? dot(sm,sp)/sqrt(den²) : 1.0
        @assert (ac<1.00001) && (ac>-1.00001) "ac=$ac"
        c =  sign(sm[1]*sp[2]-sm[2]*sp[1])/acos(min(1.0,max(-1.0, ac)))
        push!(angles, c)
    end
    return angles
end

"""
    splitblob(b::Blob, p1::CartesianIndex, p2::CartesianIndex)

Split a Blob by drawing a line from p1 to p2 (assumed to be on the perimeter
or outside b) and reblobining.
"""
function splitblob(b::Blob, p1::CartesianIndex, p2::CartesianIndex)
    function drawline(b, x0::Int, y0::Int, x1::Int, y1::Int)
        dx, sx, xa = abs(x1 - x0), x0<x1 ? 1 : -1, x0
        dy, sy, ya = -abs(y1 - y0), y0<y1 ? 1 : -1, y0
        err = dx + dy
        b[xa,ya]=false
        while !((xa==x1) && (ya==y1))
            e2 = 2err
            if e2 >= dy
                err += dy
                xa += sx
                b[xa,ya]=false
            end
            if e2 <= dx
                err += dx
                ya += sy
                b[xa,ya]=false
            end
        end
    end
    mask = copy(b.mask)
    # Draw a line to divide the particles for reblobbing
    drawline(mask, p1.I..., p2.I...)
    res=blob(mask, p->p)
    # Fix up the hoz and vert
    return map(b2->Blob(__offset(b2.bounds, b.bounds), b2.mask),res)
end

"""
    separate(b::Blob, concavity=0.5)::Vector{Blob}

Break b into many Blob(s) by looking for concave regions on the perimeter
and joining them with a short(ish) line.
"""
function separate(b::Blob, concavity=0.42, withinterior=true)::Vector{Blob}
    modn(bb,i) = (i+length(bb.psteps)-1) % length(bb.psteps) + 1
    # In the concatity vector, find peaks of positive regions
    function findmaxes(c::Vector{Float64}, minwidth::Int, concav::Float64)::Vector{Int}
        first, st, maxi, maxes = c[1] > 0.0 ? -1 : 0, 0, 0, Int[]
        for i in eachindex(c)
            if first < 0 # First keeps track of (possible) initial region above zero
                first = c[i] > c[-first] ? -i : (c[i] < 0.0 ? -first : first)
            end
            if st == 0 # In region below zero concavity
                if c[i] > 0.0
                    st, maxi = i, i
                end
            else # In region of positive concavity
                maxi = c[i] > c[maxi] ? i : maxi
                if c[i] < 0.0
                    if i - st >= minwidth && c[maxi] >= concav
                        push!(maxes, maxi)
                    end
                    st, maxi = 0, 0 # start new region of below zero concavity
                end
            end
        end
        if first > 0 # Deal with initial region
            lenfirst = (st ≠ 0 ? length(c) - st : 0) + findfirst(i -> c[i] < 0.0, eachindex(c))
            maxi = st ≠ 0 ? (c[maxi] > c[first] ? maxi : first) : first
            if lenfirst >= minwidth && c[maxi] >= concav
                push!(maxes, maxi)
            end
        end
        return maxes
    end
    besti = -1
    if length(b.psteps) > 20
        p, c = perimeter(b), curvature(b, 4)
        maxes = findmaxes(c, 4, concavity) # maximum
        if withinterior
            p, c = perimeter(b), curvature(b, 4)
            maxes = findmaxes(c, 4, concavity) # maximum
            for ir in interiorregions(b)
                if area(ir)>20
                    append!(p,perimeter(ir))
                    c = 0.0 .- curvature(ir, 4)
                    append!(maxes, findmaxes(c, 4, concavity))
                end
            end
        end
        # Find pairs of concavities to use as splitters
        bestj, bestlen = -1, 100000^2
        for (i, ii) in enumerate(maxes), ij in maxes[i+1:end]
            len = dot(p[ii].I .- p[ij].I, p[ii].I .- p[ij].I)
            # pick the shortest splitter
            if len < bestlen
                besti, bestj, bestlen = ii, ij, len
            end
        end
    end
    # There is a split point so split it and recursively separate the splits
    return besti==-1 ? [ b ] : mapreduce(separate, append!, splitblob(b, p[besti], p[bestj]),init=Blob[])
end


"""
    scorer(b::Blob)

A default function to score a blob as a candidate particle.  Smaller scores are more particle like.
"""
scorer(b::Blob, minarea=100) = # perimeter/π == ecd for a circle
     (area(b) < minarea ? 100.0 : minarea/area(b)) + perimeterlength(b) / (π*ecd(b, false))


"""
    multiseparate(img::Array, threshes, score; concavity=0.42)

Uses multiple thresholds to attempt to find the best separation of the distinct blobs in the image.  The best blob b
is defined as the one that produces particles that produce large `score(b)`.  So a blob will be split if splitting
the blob will produce multiple blobs of lower scores. The default function 'scorer(b::Blob)' looks for more circular
blobs.
"""
function multiseparate(img::Array, threshes; score=scorer, concavity=0.42, minarea=10, diag=true)
    best = Blob[]
    basefn,cx = "C:\\Users\\nritchie\\Desktop\\EGOS Only\\scoring\\multisep", 0
    for th in threshes
        starters = blob(img, p->p>=th)
        if length(starters)>0
            blobs = filter(b->area(b) > minarea, mapreduce(b->separate(b, concavity), append!, starters))
            newbest = Blob[]
            # compare the new blobs to the best previous ones
            for bb in best
                becomes = filter(b->crosscor(b, bb) > 0.1, blobs)
                deleteat!(blobs, map(b1->findfirst(b2->b1===b2, blobs), becomes))
                if length(becomes)==1  # Take the larger one...
                    push!(newbest, area(becomes[1]) > area(bb) ? becomes[1] : bb)
                # Should current best be split?
                elseif length(becomes)>1
                    split = mean(score.(becomes)) < score(bb)
                    open("$basefn[details].txt","a") do io
                        write(io, "$(basename(basefn)): before[$(cx+=1), $(score(bb))] => becomes[$(score.(becomes)))]\n")
                        write(io, "pl = $(perimeterlength(bb)), ecd=$(ecd(bb,false))\n")
                    end
                    FileIO.save(File(format"PNG","$basefn[$cx, before, $split].png"), NeXLParticle.colorizedimage([bb], img))
                    FileIO.save(File(format"PNG","$basefn[$cx, after, $split].png"), NeXLParticle.colorizedimage(becomes, img))
                    if split
                        append!(newbest, becomes) # split it up...
                    end
                else
                    push!(newbest, bb)  # and still champion...
                end
                # delete 'becomes' from 'blobs'
            end
            best = append!(blobs, newbest)
        end
    end
    return best
end



"""
   area(b::Blob)

Area of the Blob in pixel count.
"""
area(b::Blob) = count(b.mask)

"""
    maskedimage(b::Blob, img::Matrix, mark=missing, markvalue=0.5)

Extract the image data in `img` associate the the Blob `b`.
"""
function maskedimage(b::Blob, img::Matrix, mark=missing, markvalue=0.5)
    trimmed = img[b.bounds]
    res = map!(i -> b.mask[i] ? trimmed[i] : 0, zeros(eltype(trimmed), size(trimmed)), eachindex(trimmed))
    if !ismissing(mark)
        res[mark] = markvalue
    end
    return res
end

"""
    colorizedimage(bs::Vector{Blob}, img::AbstractArray)
    colorizedimage(chords::Vector{Vector{CartesianIndex}}, img::AbstractArray)

Create a colorized version of img and draw the blob or chords on it.
"""
function colorizedimage(bs::Vector{Blob}, img::AbstractArray)
    off(ci, bs) = [ci.I[i]+bs.bounds.indices[i].start-1 for i in eachindex(ci.I)]
    colors = convert.(RGB, distinguishable_colors(
        length(bs)+2,
        Color[RGB(253 / 255, 255 / 255, 255 / 255), RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
        transform = deuteranopic,
    )[3:end])
    res = RGB.(img)
    for (i, blob) in enumerate(bs)
        col = colors[i]
        foreach(ci->res[ci] = 0.5*col+0.5*img[ci], filter(c->blob[c], CartesianIndices(blob))) # draw interior
        foreach(ci->res[off(ci,blob)...] = col, perimeter(blob)) # draw perimeter
        res[off(blob.pstart,blob)...] = RGB(1.0,0.0,0.0) # draw start of perimeter...
    end
    return res
end


"""
    intersect(b1::Blob, b2::Blob)

A CartesianIndices with the region in common between b1 and b2.
"""
Base.intersect(b1::Blob, b2::Blob) =
    CartesianIndices( tuple(map(i->intersect(b1.bounds.indices[i], b2.bounds.indices[i]), eachindex(b1.bounds.indices))...))

"""
    crosscorr(b1::Blob, b2::Blob)::Float64

Measures the extent to which `b1` and `b2` represent the same region on the
image. Normalized by area of b2.
"""
StatsBase.crosscor(b1::Blob, b2::Blob)::Float64 =
    count(ci->b1[ci] && b2[ci], intersect(b1,b2))/area(b2)

"""
    soille_watershed(img::Matrix, mask::BitArray{2}, connectity4::Bool = true)

Implements the watershed algorithm described in Soille, Pierre, and Luc M. Vincent.
"Determining watersheds in digital pictures via flooding simulations." Lausanne, D.L.
International Society for Optics and Photonics, 1990
"""
function soille_watershed(img::Matrix, mask::BitArray{2}, connectity4::Bool = true)
    @assert size(img)==size(mask)
    onimg(ci::CartesianIndex)::Bool =
        (ci[1] >= 1) && (ci[1] <= height(img)) && (ci[2] >= 1) && (ci[2] <= width(img))
    function neighbors4(ci::CartesianIndex)
        res = (
            CartesianIndex(ci[1] - 1, ci[2]),
            CartesianIndex(ci[1], ci[2] - 1),
            CartesianIndex(ci[1], ci[2] + 1),
            CartesianIndex(ci[1] + 1, ci[2]),
        )
        return (ci[1] > 1) && (ci[2] > 1) && (ci[1] < height(img)) && (ci[2] < width(img)) ? res : filter(onimg, res)
    end
    function neighbors8(ci::CartesianIndex)
        res = (
            CartesianIndex(ci[1] - 1, ci[2] - 1),
            CartesianIndex(ci[1] - 1, ci[2]),
            CartesianIndex(ci[1] - 1, ci[2] + 1),
            CartesianIndex(ci[1], ci[2] - 1),
            CartesianIndex(ci[1], ci[2] + 1),
            CartesianIndex(ci[1] + 1, ci[2] - 1),
            CartesianIndex(ci[1] + 1, ci[2]),
            CartesianIndex(ci[1] + 1, ci[2] + 1),
        )
        return (ci[1] > 1) && (ci[1] < height(img)) && (ci[2] > 1) && (ci[2] < width(img)) ? res : filter(onimg, res)
    end
    neighbors = connectity4 ? neighbors4 : neighbors8
    INQUEUE, MASK, INIT, WSHED = -3, -2, -1, 0
    tabLabels = fill!(zeros(Int16, size(img)), INIT)
    currentLabel, flag = 0, false
    pixelList = Tuple{CartesianIndex,typeof(img[1, 1])}[]
    for r = 1:height(img), c = 1:width(img)
        if mask[r, c]
            push!(pixelList, (CartesianIndex(r, c), img[r, c]))
        end
    end
    # Ascending order by intensity
    sort!(pixelList, lt = (p1, p2) -> p1[2] < p2[2])
    fifo = Deque{CartesianIndex}()
    currentIndex, heightIndex1, heightIndex2 = 1, 1, 1
    while currentIndex <= length(pixelList)
        h = pixelList[currentIndex][2]
        for pixelIndex = heightIndex1:length(pixelList)
            pij1, v1 = pixelList[pixelIndex]
            if v1 ≠ h
                heightIndex1 = pixelIndex
                break
            end
            tabLabels[pij1] = MASK
            for cuv1 in neighbors(pij1)
                # initialize queue with neighbors at level h of current basins or watersheds
                if (tabLabels[cuv1] >= WSHED) && mask[cuv1]
                    push!(fifo, pij1)
                    tabLabels[pij1] = INQUEUE
                    break
                end
            end
        end
        while !isempty(fifo)
            pij2 = popfirst!(fifo)
            @assert pij2 isa CartesianIndex
            for cuv2 in neighbors(pij2)
                # labeling current point by inspecting neighbors
                if mask[cuv2]
                    if (tabLabels[cuv2] > 0)  # i.e. the pixel belongs to an already labeled basin
                        if (tabLabels[pij2] == INQUEUE) || ((tabLabels[pij2] == WSHED) && flag)
                            tabLabels[pij2] = tabLabels[cuv2]
                        elseif (tabLabels[pij2] > 0) && (tabLabels[pij2] != tabLabels[cuv2])
                            tabLabels[pij2] = WSHED
                            flag = false
                        end
                    elseif (tabLabels[cuv2] == WSHED)
                        if tabLabels[pij2] == INQUEUE
                            tabLabels[pij2] = WSHED
                            flag = true
                        end
                    elseif tabLabels[cuv2] == MASK
                        tabLabels[cuv2] = INQUEUE
                        push!(fifo, cuv2)
                    end
                end
            end
        end
        # check for new minima at level h
        pixelIndex = heightIndex2
        while pixelIndex <= length(pixelList)
            pij3, v3 = pixelList[pixelIndex]
            if v3 ≠ h
                # this pixel is at level h+1
                heightIndex2 = pixelIndex
                break
            end
            if tabLabels[pij3] == MASK  # the pixel is inside a new minimum
                currentLabel += 1
                push!(fifo, pij3)
                tabLabels[pij3] = currentLabel
                while isempty(fifo) == false
                    p2 = popfirst!(fifo)
                    for cuv3 in neighbors(p2) # inspect neighbors of p2
                        if (tabLabels[cuv3] == MASK) && mask[cuv3]
                            push!(fifo, cuv3)
                            tabLabels[cuv3] = currentLabel
                        end
                    end
                end
            end
            currentIndex += 1
            pixelIndex += 1
        end
    end
    fp = zeros(UInt8, size(img))
    for c = 1:width(img), r = 1:height(img)
        fp[r, c] = tabLabels[r, c] == INIT ? 0 : tabLabels[r, c]
    end
    return fp
end

function interiorregions(b::Blob)
    onedge(b1, b2) = any(map(i->b1[i].start==1 ||  b1[i].stop==b2[i].stop-b2[i].start+1, eachindex(b2)))
    tmp = filter(ib->!onedge(ib.bounds.indices, b.bounds.indices), blob(b.mask,p->!p))
    return Blob[ Blob(__offset(tb.bounds, b.bounds), tb.mask) for tb in tmp ]
end

filledarea(b::Blob) = area(b) + sum(area.(interiorregions(b)))

function __offset(ci::CartesianIndices, base::CartesianIndices)
    rs = map(i->ci.indices[i].start + base.indices[i].start -1:ci.indices[i].stop + base.indices[i].start -1, eachindex(ci.indices))
    return CartesianIndices(tuple(rs...))
end
