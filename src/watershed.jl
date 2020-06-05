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
