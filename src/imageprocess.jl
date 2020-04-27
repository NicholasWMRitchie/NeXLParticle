using Images
using FileIO

struct Blob
    vertical::UnitRange{Int}
    horizontal::UnitRange{Int}
    mask::Matrix{Bool}
end

area(b::Blob) = count(b.mask)

function image(b::Blob, img::Matrix)
    trimmed = img[b.vertical, b.horizontal]
    map!(i -> b.mask[i] ? trimmed[i] : 0, zeros(eltype(trimmed), size(trimmed)), eachindex(trimmed))
end

function blob(img::Matrix, thresh)::Vector{Blob}
    function extractmask(res, i, minR, maxR, minC, maxC)
        mask = zeros(Bool, (maxR - minR + 1, maxC - minC + 1))
        for c = minC:maxC, r = minR:maxR
            mask[r-minR+1, c-minC+1] = (res[r, c] == i)
        end
        return mask
    end
    imgH, imgW = size(img)
    alias = Vector{Set{Integer}}()
    res = zeros(UInt8, size(img))
    prev, next = img[1, 1], 0
    for r = 1:imgH
        for c = 1:imgW
            if thresh(img[r, c])
                above = r>1 ? res[r-1, c] : 0
                if prev == 0
                    if above == 0
                        prev = (next += 1)
                        push!(alias, Set{Int}(prev))
                    else
                        prev = above
                    end
                elseif (above ≠ 0) && (above ≠ prev)
                    # Same region => different indices
                    ia = findfirst(al -> above in al, alias)
                    ip = findfirst(al -> prev in al, alias)
                    @assert !isnothing(ia) "$above not in $alias"
                    @assert !isnothing(ip) "$prev not in $alias"
                    if ia ≠ ip
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
        #@show currentIndex
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

"""
    insidemask(img)

Determine which pixels are inside a watershed wall.  Returns a BitArray with the pixels inside a wall.
"""
function insidemask(img)
  function pixelinside(r0, c0, img)
    cx1, cx2=0, 0
    for c1 in 1:c0
      if img[r0,c1] == 0 # skip to first wall
        prev=0
        for c2 in c1+1:c0
          if img[r0,c2] ≠ prev
            cx1+=1
            prev=img[r0,c2]
          end
        end
        break
      end
    end
    for c1 in width(img):-1:c0
      if img[r0,c1] == 0 # skip to first wall
        prev=0
        for c2 in c1-1:-1:c0
          if img[r0,c2] ≠ prev
            cx2+=1
            prev=img[r0,c2]
          end
        end
        break
      end
    end
    return (cx1%2==1) && (cx2%2==1)
  end
  pi1=BitArray(undef, size(img))
  for c in 1:width(img), r in 1:height(img)
    pi1[r,c]=pixelinside(r,c,img)
  end
  return pi1
end
