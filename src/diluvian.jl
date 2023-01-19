using Colors
import OnlineStats
using DataFrames

# Islands are collects of adjacent CartesianIndex with their population count.  The first CartesianIndex is maintained
# as the one with the largest population. Internal use only.
Island = Vector{Tuple{CartesianIndex,Int}}

# This function converts the data values into bin index
defbin(x) = floor(Integer, min(1.0, max(0.0, x)) * 25.0)

# Are two CartesianIndex(s) adjacent?
adj(ci1::CartesianIndex, ci2::CartesianIndex) = sum(abs(ci1.I[i] - ci2.I[i]) for i in eachindex(ci1.I)) <= 1
# Is this CartesianIndex adjacent to the island cis2?
adj(ci::CartesianIndex, cis2::Island) = any(adj(ci, ci2[1]) for ci2 in cis2)
# Are two islands of CartesianIndex(s) adjacent by any CI within an island being adjacent to a CI from the other island?
adj(cis1::Island, cis2::Island) = any(adj(ci[1], cis2) for ci in cis1)

# The height of an island is the count in the first (highest) bin
height(isl::Island) = isl[1][2]

# Are two islands of heights `h1` and `h2` connected by the index `ci` of height `hi`?
defconnects(ci::CartesianIndex, hi::Int, h1::Int, h2::Int, σ::Float64 = 3.0) = hi > min(h1, h2) - σ * sqrt(min(h1, h2))


"""
Clusters N-dimensional data using an algorithm based on a N-dimensional histogram.  The data is summarized in an
N-dimensional histogram.  The bins are sorted according to population count.  Starting with the most populous bin,
adjacent bins are joined into islands according to the following rules.  If a bin is adjacent to one island then the
bin is added to the island.  If the bin is adjacent to two islands and it's population exceeds a threshold of the
highest bin in each island, then the islands are joined together into a larger island.  If the bin height is less than
the threshold, the bin is added to the less island with the smaller highest bin.  So constructed, each island represents
a cluster.  The clusters are assigned integer indexes starting at 1 and cluster assignments for each pixel are written
to an `Array`.

This algorithm is best for data like compositional data which is naturally normalized to a unity sum and for which
each layer in the data is on an equivalent scale.  The algorithm can be extended to work with other data types or the
data can be preprocessed to bring it onto an equivalent scale.
"""
struct DiluvianCluster
  labels::Vector
  # An Vector of planes of data Arrays (one for each label)
  data::Vector{Array}
  # An Array containing the cluster assignment
  clusters::Array{Int16}
  # Cached cluster count
  count::Int

  function DiluvianCluster(df::AbstractDataFrame, bin::Function = defbin, connects::Function = defconnects)
    nms = names(df)
    return DiluvianCluster(nms, [ df[:, name] for name in nms ], bin=bin, connects=connects)
  end

  function DiluvianCluster(
    labels::Vector,
    data::Vector{<:AbstractArray};
    bin::Function = defbin,
    connects::Function = defconnects,
  )
    @assert length(labels) == length(data) "length(labels) must equal length(data)"
    @assert all(size(data[1]) == size(datum) for datum in data[2:end]) "The arrays in data must all be the same size."
    islands = flood(buildhistogram(data, bin), connects)
    clusters = buildclusters(islands, data, bin)
    return new(labels, data, clusters, length(islands))
  end

  function DiluvianCluster(labeleddata::Tuple{Any,AbstractArray}...; bin::Function = defbin, connects::Function = defconnects)
    lbls = [lbl for (lbl, _) in labeleddata]
    data = [arr for (_, arr) in labeleddata]
    return DiluvianCluster(lbls, data, bin = bin, connects = connects)
  end
end

Base.show(io::IO, dc::DiluvianCluster) = print(io, "DiluvianCluster[$(dc.count) clusters]")

# Bin the data in a N-dimensional histogram
function buildhistogram(data::Vector{<:AbstractArray}, bin::Function)
  binindex = (ci, data, bin) -> CartesianIndex((bin(datum[ci]) for datum in data)...)
  hist = Dict{CartesianIndex,Integer}()
  for ci in CartesianIndices(data[1])
    cc = binindex(ci, data, bin)
    hist[cc] = get(hist, cc, 0) + 1
  end
  return hist
end

# Construct islands by joining adjacent bins unless there is a valley between peaks
function flood(hist::Dict{CartesianIndex,Integer}, connects::Function)::Vector{Island}
  dist(ci, isl) = sum((ci.I .- isl[1][1].I) .^ 2)
  # Sort by bins by height (highest first)
  shist = sort([(ci, n) for (ci, n) in hist], lt = (a, b) -> a[2] > b[2])
  islands = Vector{Island}()  # Maintained in height(island) order (highest first)
  for (ci, n) in shist
    added = false
    for (i1, isl1) in enumerate(islands)
      if adj(ci, isl1)
        rm, h1 = [], height(isl1)
        # Add ci to the closest adjacent peak
        closest, cdist = isl1, dist(ci, isl1)
        # Is ci adjacent to any other islands?
        for i2 in i1+1:length(islands)
          isl2 = islands[i2]
          if adj(ci, isl2)
            if connects(ci, n, h1, height(isl2))
              append!(isl1, isl2)
              push!(rm, i2)
              # The bin ci connects isl1 to isl2 so always add it to isl1
              closest, cdist = isl1, 0
            else
              if dist(ci, isl2) < cdist
                # ci is closer to the peak of isl2 than closest
                closest, cdist = isl2, dist(ci, isl2)
              end
            end
          end
        end
        added = true
        push!(closest, (ci, n))
        deleteat!(islands, rm)
        break
      end
    end
    if !added # create a new Island which is shorter than all previous Islands
      push!(islands, [(ci, n)])
    end
  end
  binsum = isl -> sum(n for (_, n) in isl)
  return sort!(islands, lt = (a, b) -> binsum(a) < binsum(b), rev = true)
end

# Build an array summarizing cluster membership
function buildclusters(islands::Vector{Island}, data::Vector{<:Array}, bin::Function)::Array{UInt16}
  binindex = (ci, data, bin) -> CartesianIndex((bin(datum[ci]) for datum in data)...)
  idxforci = Dict{CartesianIndex,UInt16}(ci => i for (i, island) in enumerate(islands) for (ci, _) in island)
  return map(ci -> idxforci[binindex(ci, data, bin)], CartesianIndices(data[1]))
end

Base.length(dc::DiluvianCluster) = dc.count
Base.eachindex(dc::DiluvianCluster) = Base.OneTo(dc.count)

"""
    asimage(dc::DiluvianCluster, colors::Vector{Colorant}, other = colorant"black")::Array{Colorant}
    asimage(dc::DiluvianCluster)

Get a colorized image when the first `length(colors)` clusters are colored according to `colors` and the
remainder are `other`. The second version uses `distinguishable_colors(…)` to generate the enough colors
automatically.
"""
asimage(dc::DiluvianCluster, colors::Vector{<:Colorant}, other = colorant"black")::Array{<:Colorant} =
  map(v -> get(colors, v, other), dc.clusters)

defaultpalette(dc::DiluvianCluster) = distinguishable_colors(length(dc), colorant"yellow")

asimage(dc::DiluvianCluster) = asimage(dc, defaultpalette(dc))

"""
    asa(::Type{DataFrame}, dc::DiluvianCluster, cluster::Int)::DataFrame

Create a DataFrame with the data from each pixel in `cluster`.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, dc::DiluvianCluster, cluster::Int)::DataFrame
  cis = collect(filter(ci -> dc.clusters[ci] == cluster, CartesianIndices(dc.clusters)))
  df = DataFrame(Cluster = [cluster for _ in cis], Index = cis)
  for (i, lbl) in enumerate(dc.labels)
    datum = dc.data[i]
    insertcols!(df, "$lbl" => [datum[ci] for ci in cis])
  end
  return df
end

"""
    asmask(dc::DiluvianCluster, cluster::Int)

Get a Bool mask for the specified cluster
"""
asmask(dc::DiluvianCluster, cluster::Int)::BitArray =
  BitArray(dc.clusters[ci] == cluster for ci in CartesianIndices(dc.clusters))


"""
  clusters(dc::DiluvianCluster)::Array{UInt16}

Returns an Array containing the cluster assignment as an integer index.
"""
clusters(dc::DiluvianCluster) = dc.clusters

"""
    membercounts(dc::DiluvianCluster)

How many data points are there in each cluster?
"""
membercounts(dc::DiluvianCluster) = [count(v -> v == cl, dc.clusters) for cl = 1:length(dc)]
Base.count(dc::DiluvianCluster, i::Integer) = count(v -> v == i, dc.clusters)

function clusterstats(dc::DiluvianCluster, cluster::Int)
  stats = Dict(lbl => Series(Mean(), Variance(), Extrema()) for lbl in dc.labels)
  cis = filter(ci -> dc.clusters[ci] == cluster, CartesianIndices(dc.clusters))
  for (i, lbl) in enumerate(dc.labels)
    fit!(stats[lbl], (dc.data[i][ci] for ci in cis))
  end
  return stats
end

function clusterstats(dc::DiluvianCluster, clusters::Vector{Int})
  cls = Set(clusters)
  stats = Dict(lbl => Series(Mean(), Variance(), Extrema()) for lbl in dc.labels)
  cis = filter(ci -> dc.clusters[ci] in cls, CartesianIndices(dc.clusters))
  for (i, lbl) in enumerate(dc.labels)
    fit!(stats[lbl], (dc.data[i][ci] for ci in cis))
  end
  return stats
end

function summarizeclusters(dc::DiluvianCluster; statistics = false, sorted = true)
  dcidx = Base.OneTo(length(dc))
  stats = [Series(Mean(), Variance(), Extrema()) for cl in dcidx, lbl in dc.labels]
  for (ii, cl) in enumerate(dc.clusters), lbli in eachindex(dc.labels)
    fit!(stats[cl, lbli], dc.data[lbli][ii])
  end
  res = DataFrame(Cluster = collect(dcidx), Count = membercounts(dc))
  sl = [ (i, dc.labels[i]) for i in eachindex(dc.labels)]
  if sorted
    sort!(sl , lt=(a,b)->a[2]<b[2])
  end
  for i in eachindex(sl)
    (lbli, lbl) = sl[i]
    insertcols!(res, "$lbl" => [OnlineStats.value(stats[cl, lbli].stats[1]) for cl in dcidx])
    if statistics
      insertcols!(res, "σ[$lbl]" => [sqrt(OnlineStats.value(stats[cl, lbli].stats[2])) for cl in dcidx])
      insertcols!(res, "min[$lbl]" => [minimum(stats[cl, lbli].stats[3]) for cl in dcidx])
      insertcols!(res, "max[$lbl]" => [maximum(stats[cl, lbli].stats[3]) for cl in dcidx])
    end
  end
  return res
end

"""
    multiternary(dc::DiluvianCluster, cluster::Int; maxitems=1000, norm=1.0, palette = nothing, variance = false)

Plot the data from the specified cluster from most significant dimensions as ternary diagram.
`maxitems` limits the total number of points plotted when there are a very large number of data points.
`norm` is useful if the data is normalized to something other than unity.
`palette` defaults to the same default palette as toimage(…) or can be specified as `Colorant[]`
`variance` determines whether maximum mean value or maximum variance is used to select elements to plot
"""
function multiternary(dc::DiluvianCluster, clusters::Vector{Int}; maxitems = 1000, norm = 1.0, palette = nothing, variance = false)
  sts = [(OnlineStats.value(s.stats[variance ? 2 : 1]), lbl) for (lbl, s) in clusterstats(dc, clusters)]
  sort!(sts, lt = (a, b) -> a[1] > b[1])
  cols = [Symbol(st[2]) for st in sts]
  df = mapreduce(cluster -> asa(DataFrame, dc, cluster), vcat, clusters)
  insertcols!(df, :Randomizer => randperm(nrow(df)))
  sort!(df, :Randomizer)
  df = nrow(df) > maxitems ? df[1:maxitems, :] : df
  if isnothing(palette)
    palette = distinguishable_colors(length(dc), colorant"yellow")[clusters]
  end
  return multiternary(df, cols, :Cluster, norm = norm, withcount = false, palette = palette)
end
