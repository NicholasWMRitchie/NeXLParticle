# These are methods to implement particle segregation - separating particles that touch for counting etc.


struct Watershed
    sz::Tuple
    steps:: LinRange
    bins::Vector{Vector{CartesianIndex}}

    function Watershed(data::Array, steps::LinRange)
        indexof(x,lr::LinRange)=min(lr.len,max(1,trunc(Int,(lr.lendiv*(x-lr.start))/(lr.stop-lr.start))+1))
        bins = [ CartesianIndex[] for i in 1:length(steps) ]
        foreach(ci->push!(bins[indexof(data[ci], steps)], ci), CartesianIndices(data))
        return new(size(data), steps, bins)
    end
end

bin(ws::Watershed, idx::Int) = ( ws.steps[idx-1], ws.steps[idx] )
NeXLSpectrum.depth(ws::Watershed) = length(ws.steps)

function masks(ws::Watershed)::BitArray
    prev = fill(true, ws.sz)
    res = BitArray(undef, (ws.sz..., length(ws.steps)))
    for (i, vci) in enumerate(ws.bins)
        foreach(ci->prev[ci]=false,vci)
        res[:,:,i]=prev
    end
    return res
end

function areas(ws::Watershed)::Vector{Int}
    res = [size(vci,1) for vci in reverse(ws.bins)]
    sum = 0
    for i in eachindex(res)
        res[i] = (sum += res[i])
    end
    return res
end
