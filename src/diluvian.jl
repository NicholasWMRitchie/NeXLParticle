
using Colors
using OnlineStats
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

# Are two islands connected by the index `ci`?  Islands are connected if `ci` is adjacent to both and `n` is at least
# thresh times the population of the height of both `isl1` and `isl2`.
defconnects(ci::CartesianIndex, n::Int, isl1::Island, isl2::Island, thresh = 0.5) =
  adj(ci, isl1) && adj(ci, isl2) && (n > thresh * height(isl1)) && (n > thresh * height(isl2))


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

  function DiluvianCluster(
    labels::Vector,
    data::Vector{<:Array};
    bin::Function = defbin,
    connects::Function = defconnects,
  )
    @assert length(labels) == length(data) "length(labels) must equal length(data)"
    @assert all(size(data[1])==size(datum) for datum in data[2:end]) "The arrays in data must all be the same size."
    islands = flood(buildhistogram(data, bin), connects)
    clusters = buildclusters(islands, data, bin)
    return new(labels, data, clusters, length(islands))
  end

  function DiluvianCluster(labeleddata::Tuple{Any,Array}...; bin::Function=defbin, connects::Function = defconnects)
    lbls = [ lbl for (lbl, arr) in labeleddata ]
    data = [ arr for (lbl, arr) in labeleddata ]
    return DiluvianCluster(lbls, data, bin=bin, connects=connects)
  end
end

Base.show(io::IO, dc::DiluvianCluster) = print(io, "DiluvianCluster[$(dc.count) clusters]")

# Bin the data in a N-dimensional histogram
function buildhistogram(data::Vector{<:Array}, bin::Function)
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
  # Sort by height
  shist = sort([(ci, n) for (ci, n) in hist], lt = (a, b) -> a[2] > b[2])
  islands = Vector{Island}()
  for (ci, n) in shist
    added = false
    for i1 in eachindex(islands)
      if added
        break
      end
      for i2 in i1+1:length(islands)
        if added
          break
        end
        if connects(ci, n, islands[i1], islands[i2])
          # add to taller and join islands
          added = true
          addto, delfrom = height(islands[i1]) > height(islands[i2]) ? (i1, i2) : (i2, i1)
          append!(islands[addto], islands[delfrom])
          push!(islands[addto], (ci, n))
          deleteat!(islands, delfrom)
        else
          a1, a2 = adj(ci, islands[i1]), adj(ci, islands[i2])
          added = a1 || a2
          if a1 && a2
            # Add to less high island
            push!(islands[height(islands[i1]) < height(islands[i2]) ? i1 : i2], (ci, n))
          elseif a1 # Add to island 1
            push!(islands[i1], (ci, n))
          elseif a2 # Add to island 2
            push!(islands[i2], (ci, n))
          end
        end
      end
    end
    if !added # Create a e=new island
      push!(islands, [(ci, n)])
    end
  end
  pixelcount = (island) -> sum( n for (_, n) in island)
  return sort!(islands, lt = (a, b) -> pixelcount(a) > pixelcount(b))
end

# Build an array summarizing cluster membership
function buildclusters(islands::Vector{Island}, data::Vector{<:Array}, bin::Function)::Array{UInt16}
  binindex = (ci, data, bin) -> CartesianIndex((bin(datum[ci]) for datum in data)...)
  idxforci = Dict{CartesianIndex,UInt16}()
  for i in eachindex(islands), ci in islands[i]
    idxforci[ci[1]] = i
  end
  return map(ci -> idxforci[binindex(ci, data, bin)], CartesianIndices(data[1]))
end

Base.length(dc::DiluvianCluster) = dc.count

"""
    asimage(dc::DiluvianCluster, colors::Vector{Colorant}, other = colorant"black")::Array{Colorant}
    asimage(dc::DiluvianCluster)

Get a colorized image when the first `length(colors)` clusters are colored according to `colors` and the
remainder are `other`. The second version uses `distinguishable_colors(…)` to generate the enough colors
automatically.
"""
asimage(dc::DiluvianCluster, colors::Vector{<:Colorant}, other = colorant"black")::Array{<:Colorant} =
  map(v -> get(colors, v, other), dc.clusters)

asimage(dc::DiluvianCluster) = asimage(dc, distinguishable_colors(length(dc), colorant"yellow"))

"""
    asa(::Type{DataFrame}, dc::DiluvianCluster, cluster::Int)::DataFrame

Create a DataFrame with the data from each pixel in `cluster`.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, dc::DiluvianCluster, cluster::Int)::DataFrame
  cis = collect(filter(ci->dc.clusters[ci]==cluster, CartesianIndices(dc.clusters)))
  df = DataFrame(Cluster = [ cluster for _ in cis ], Index=cis)
  for (i, lbl) in enumerate(dc.labels)
    datum = dc.data[i]
    insertcols!(df, "$lbl" => [ datum[ci] for ci in cis ])
  end
  return df
end

"""
    asmask(dc::DiluvianCluster, cluster::Int)

Get a Bool mask for the specified cluster
"""
asmask(dc::DiluvianCluster, cluster::Int) = map(v -> v == cluster, dc.clusters)


"""
  clusters(dc::DiluvianCluster)::Array{UInt16}

Returns an Array containing the cluster assignment as an integer index.
"""
clusters(dc::DiluvianCluster) = dc.clusters

"""
    membercounts(dc::DiluvianCluster)

How many data points are there in each cluster?
"""
membercounts(dc::DiluvianCluster) =
  [ count(v -> v == cl, dc.clusters) for cl in 1:length(dc) ]


function clusterstats(dc::DiluvianCluster, cluster::Int)
  stats = Dict(lbl => Series(Mean(), Variance(), Extrema()) for lbl in dc.labels)
  cis = filter(ci->dc.clusters[ci]==cluster, CartesianIndices(dc.clusters))
  for (i, lbl) in enumerate(dc.labels)
    fit!(stats[lbl], (dc.data[i][ci] for ci in cis))
  end
  return stats
end

function clusterstats(dc::DiluvianCluster, clusters::Vector{Int})
  cls = Set(clusters)
  stats = Dict(lbl => Series(Mean(), Variance(), Extrema()) for lbl in dc.labels)
  cis = filter(ci->dc.clusters[ci] in cls, CartesianIndices(dc.clusters))
  for (i, lbl) in enumerate(dc.labels)
    fit!(stats[lbl], (dc.data[i][ci] for ci in cis))
  end
  return stats
end

function summarizeclusters(dc::DiluvianCluster; statistics=false)
  len=length(dc)
  stats = [ Series(Mean(), Variance(), Extrema()) for cl in 1:length(dc), lbl in dc.labels ]
  for (ii, cl) in enumerate(dc.clusters), lbli in eachindex(dc.labels)
    fit!(stats[cl,lbli],dc.data[lbli][ii])
  end
  res = DataFrame(Cluster=collect(1:len), Count=membercounts(dc))
  for lbli in eachindex(dc.labels)
    lbl = dc.labels[lbli]
    insertcols!(res,
      "$lbl" => [ OnlineStats.value(stats[cl, lbli].stats[1]) for cl in 1:len ])
    if statistics
      insertcols!(res,
        "σ[$lbl]" => [ sqrt(OnlineStats.value(stats[cl, lbli].stats[2])) for cl in 1:len ])
      insertcols!(res,
        "min[$lbl]" => [ minimum(stats[cl, lbli].stats[3]) for cl in 1:len ])
      insertcols!(res,
        "max[$lbl]" => [ maximum(stats[cl, lbli].stats[3]) for cl in 1:len ])
    end
  end
  return res
end

"""
    multiternary(dc::DiluvianCluster, cluster::Int)

Plot the data from the specified cluster from most significant dimensions as ternary diagram.
"""
function multiternary(dc::DiluvianCluster, clusters::Vector{Int}; maxitems=1000, norm=1.0)
  sts = [ ( OnlineStats.value(s.stats[1]), lbl) for (lbl, s) in clusterstats(dc, clusters) ]
  sort!(sts, lt=(a,b)->a[1]>b[1])
  cols = [ Symbol(st[2]) for st in sts ]
  df = mapreduce(cluster->asa(DataFrame,dc,cluster), vcat, clusters)
  insertcols!(df,:Randomizer=>randperm(nrow(df)))
  sort!(df, :Randomizer)
  df = nrow(df)>maxitems ? df[1:maxitems, :] : df
  return multiternary(df, cols, :Cluster, norm=norm)
end
