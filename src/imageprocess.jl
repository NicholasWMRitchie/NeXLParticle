using Images
using FileIO
using LinearAlgebra


"""
A Blob is a mask consisting of blocks of adjacent pixels meeting a threshold.
Only those pixels immediately above, below, left or right are considered to
be adjacent (not diagonals).  Blobs offer easy access to a mask consisting
of pixels in the blob and the perimeter either as a list of steps from
one perimeter point to the next or a list of CartesianIndex with the coordinates
of the perimeter points.
"""
struct Blob
    vertical::UnitRange{Int} # vertical extent of the blob in the original image
    horizontal::UnitRange{Int} # horizontal extent of the blob in the original image
    mask::BitArray  # Mask of those pixels in the blob
    pstart::CartesianIndex # Start of the perimeter
    psteps::Vector{Tuple{Int,Int}} # Steps around perimeter

    function Blob(vert::UnitRange{Int}, horz::UnitRange{Int}, mask::BitArray)
        function perimeter(m)  # Computes the perimeter steps
            stps = Tuple{Int,Int}[]
            if prod(size(m))>1
                msk(r,c) = checkbounds(Bool, m, r, c) && m[r,c]
                mod8(x) = (x + 7) % 8 + 1 # maintains 1...8
                steps = ( (0, 1), (-1, 1), (-1, 0), (-1, -1), #
                          (0, -1), (1, -1), (1, 0), (1, 1) )
                pre, dirs = ( findfirst(r->m[r,1],1:size(m,1)), 1 ), 7:10
                @assert msk(pre...)
                ff = findfirst(r->msk((pre .+ steps[mod8(r)])...), dirs)
                @assert !isnothing(ff)
                prevdir = mod8(dirs[ff])
                start=pre .+ steps[prevdir]
                @assert msk(start...)
                curr = start
                while true
                    dirs = prevdir-2:prevdir+4
                    nextdir = mod8(dirs[findfirst(r->msk((curr .+ steps[mod8(r)])...), dirs)])
                    push!(stps, steps[nextdir])
                    curr = curr .+ steps[nextdir]
                    @assert msk(curr...)
                    if curr == start
                        break
                    end
                    prevdir = nextdir
                end
                return ( CartesianIndex(start), stps)
            else
                return ( CartesianIndex(1,1), stps)
            end
        end
        @assert ndims(mask)==2
        @assert length(horz)==size(mask,2)
        @assert length(vert)==size(mask,1)
        return new(vert, horz, mask, perimeter(mask)...)
    end
end


"""
    blob(img::AbstractArray, thresh::Function)::Vector{Blob}

Create a Vector of Blob containing discontiguous regions meeting the
threshold function.
"""
function blob(img::AbstractArray, thresh::Function)::Vector{Blob}
    function extractmask(res, i, minR, maxR, minC, maxC)
        mask = BitArray(undef, maxR - minR + 1, maxC - minC + 1)
        for c in minC:maxC, r in minR:maxR
            mask[r-minR+1, c-minC+1] = (res[r, c] == i)
        end
        return mask
    end
    imgH, imgW = size(img)
    alias = Vector{Set{Integer}}()
    res = zeros(UInt16, size(img))
    prev, next = img[1, 1], 0
    for r = 1:imgH
        for c = 1:imgW
            if thresh(img[r, c])
                above = r>1 ? res[r-1, c] : 0
                if prev == 0
                    if above == 0
                        prev = (next += 1) # new region
                        push!(alias, Set{Int}(prev))
                    else
                        prev = above
                    end
                elseif (above ≠ 0) && (above ≠ prev)
                    # Same blob => different indexes
                    ia = findfirst(al -> above in al, alias)
                    ip = findfirst(al -> prev in al, alias)
                    @assert !isnothing(ia) "$above not in $alias"
                    @assert !isnothing(ip) "$prev not in $alias"
                    if ia ≠ ip # merge regions
                        union!(alias[ia], alias[ip])
                        deleteat!(alias, ip)
                    end
                    prev = above
                end
                @assert prev != 0
                res[r, c] = prev
            else
                prev = 0
            end
        end
    end
    # Second pass combine the adjacent indices
    newidx, nblobs = zeros(UInt8, next), length(alias)
    for i = 1:nblobs
        for id in alias[i]
            newidx[id] = i
        end
    end
    rects = Dict{Int,NTuple{4,Int}}(i => (100000, -1, 100000, -1) for i = 1:nblobs)
    for r = 1:imgH, c = 1:imgW
        bidx = res[r, c]
        if bidx ≠ 0
            newi = newidx[bidx]
            res[r, c] = newi
            rect = rects[newi]
            rects[newi] = (min(r, rect[1]), max(r, rect[2]), min(c, rect[3]), max(c, rect[4]))
        end
    end
    res = [Blob(rect[1]:rect[2], rect[3]:rect[4], extractmask(res, i, rect...)) for (i, rect) in rects]
    sort!(res, lt=(b1,b2) -> area(b1) > area(b2))
    return res
end

"""
    perimeter(b::Blob)::Vector{CartesianIndex}

Returns a vector of `CartesianIndex` corresponding to the points around the
perimeter of the blob.
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
    ecd(b::Blob)

Computes the equivalent circular diameter.
"""
ecd(b::Blob) = 2.0 * sqrt(area(b) / π)

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
        ac = dot(sm,sp)/(sqrt(dot(sm,sm))*sqrt(dot(sp,sp)))
        @assert (ac<1.000001) && (ac>-1.00001)
        c = (sm[1]*sp[2]-sm[2]*sp[1] < 0.0 ? -1.0 : 1.0) / #
            acos(min(1.0,max(-1.0, ac)))
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
    offset(ui, off) = ui.start+off-1:ui.stop+off-1
    mask = copy(b.mask)
    # Draw a line to divide the particles for reblobbing
    drawline(mask, p1.I..., p2.I...)
    res=blob(mask, p->p)
    # Fix up the hoz and vert
    return map(b2->Blob(offset(b2.vertical, b.vertical.start), #
                       offset(b2.horizontal, b.horizontal.start), b2.mask),res)
end

"""
    separate(b::Blob, concavity=0.5)::Vector{Blob}

Break b into many Blob(s) by looking for concave regions on the perimeter
and joining them with a short(ish) line.
"""
function separate(b::Blob, concavity=0.5)::Vector{Blob}
    modn(i) = (i+length(b.psteps)-1) % length(b.psteps) + 1
    stepsum(itr) = mapreduce(j->b.psteps[modn(j)], (x,y)->.+(x,y), itr, init=(0,0))
    # point where the line from p0->p1 intersects the line from p2->p3
    function intersection(p0, p1, p2, p3)
        function check(t)
            s = (t*d01[1] - d02[1])/d23[1]
            i1, i2 = (p2 .- s .* d23), (p0 .- t .* d01)
            all(map(x->isapprox(x,0.0,atol=1.0e-8), i1 .- i2))
        end
        d23, d01, d02 = p2 .- p3, p0 .- p1, p0 .- p2
        den = d23[2]*d01[1]-d23[1]*d01[2]
        if isapprox(den, 0.0, atol=1.0e-8) # parallel
            # check distance between parallel lines (a*x+b*y+c=0 form)
            aa, bb = d01[2], d01[1]
            c1, c2 = -aa*p0[1] - bb*p0[2], -aa*p2[1] - bb*p2[2]
            sep2 = (c1-c2)^2/(aa^2+bb^2)
            if sep2 <= dot(p0 .- p2, p0 .- p2) && #
                (sqrt(sep2) < 0.1*max(length(b.horizontal),length(b.vertical)))
                # parallel, offset by less than 10%
                ii = p2
            else # parallel but not close
                ii = (-1.0, -1.0)
            end
        else
            t = (d23[2]*d02[1] - d23[1]*d02[2]) / den
            @assert check(t)
            ii = p0 .- t .* d01
        end
        return CartesianIndex(map(x->round(Int, x),ii))
    end
    besti = -1
    if length(b.psteps) > 20
        n = max(3, 4*length(b.psteps)÷100)
        p, c = perimeter(b), curvature(b, max(3, 4*length(b.psteps)÷100))
        beg = -1;
        if c[end]>concavity
            for i in length(c)-1:-1:1
                if c[i]<thresh
                    beg = i+1
                    break
                end
            end
        end
        # Find the perimeter points with the maximum concavity
        maxes = Int[] # Concave region's max
        for i in 1:(beg==-1 ? length(c) : beg-1)
            if c[i]>concavity
                beg = (beg==-1 ? i : beg)
            elseif beg ≠ -1
                mc = findmax(c[beg:i-1])
                push!(maxes, beg+mc[2]-1)
                beg=-1
            end
        end
        # Find pairs of concavities to use as splitters
        bestj, bestlen = -1, 100000^2
        for i in eachindex(maxes)
            ii = maxes[i]
            ci, mi = 1.0 .* p[ii].I, (p[modn(ii-n)].I .+ p[modn(ii+n)].I) ./ 2.0
            for j in i+1:length(maxes)
                ij = maxes[j]
                cj, mj = 1.0 .* p[ij].I, (p[modn(ij-n)].I .+ p[modn(ij+n)].I) ./ 2.0
                pi=intersection(ci, mi, cj, mj)
                # Is the intersection point inside the blob??
                if checkbounds(Bool, b.mask, pi) && b.mask[pi]
                    len = dot(p[ii].I .- p[ij].I, p[ii].I .- p[ij].I)
                    # pick the shortest splitter
                    if len < bestlen
                        besti, bestj, bestlen = ii, ij, len
                    end
                end
            end
        end
    end
    # There is a split point so split it and recursively separate the splits
    return besti==-1 ? [ b ] : #
        mapreduce(separate, append!, splitblob(b, p[besti], p[bestj]),init=Blob[])
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
    trimmed = img[b.vertical, b.horizontal]
    res = map!(i -> b.mask[i] ? trimmed[i] : 0, zeros(eltype(trimmed), size(trimmed)), eachindex(trimmed))
    if !ismissing(mark)
        res[mark] = markvalue
    end
    return res
end

"""
    crosscorr(b1::Blob, b2::Blob)::Float64

Measures the extent to which `b1` and `b2` represent the same region on the
image.
"""
function crosscorr(b1::Blob, b2::Blob)::Float64
    vo(b, r, c)::Bool = b.mask[r-b.vertical.start+1, c-b.horizontal.start+1]
    vi, hi = intersect(b1.vertical, b2.vertical), intersect(b1.horizontal, b2.horizontal)
    all(l->length(l)>0, (vi,hi)) ? #
        sum(map(c->count(map(r->all(map(b->vo(b, r, c),(b1, b2))), vi)), hi))/max(area.((b1,b2))...) : 0.0
end

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
