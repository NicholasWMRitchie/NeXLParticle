module NeXLParticle

using DataFrames
using FileIO
using Reexport
using Requires
using StatsBase
using AxisArrays
using ThreadsX

@reexport using NeXLSpectrum


include("zeppelin.jl")
export Zeppelin # struct holding a hdz/pxz file pair
export classes # An ordered list of particle classes
export header # A alphabetically sorted list of header items
export data # Provides access to the DataFrame containing the particle data
export iszeppelin # Is this IOStream or filename a Zeppelin file?
export spectrum # read a spectrum associate with a Zeppelin file (or missing)
export eachparticle # The range of partice row indices
export rowsmax, rowsmin, rowsclass
export MORPH_COLS, COMP_COLS, CLASS_COLS
export beamenergy, probecurrent, magdata
export maxparticle
export writeZep

export ParticleClassifier

include("signature.jl")
export Signature
export signature # Compures the particle signature
export quantify # Quantifies the particle spectrum data and constructs a fresh Zeppelin with the results
export NSigmaCulling
export NoCulling

include("classrule.jl")
export classify
export OrderedRuleSet
export GSRRules

include("baseRules.jl")
export BaseRules
include("nullRules.jl")
export NullRules

include("multiternary.jl")
export multiternary
export TernPalette
export TernColorblind

include("blob.jl")
export Blob
export blob
export separate # Split agglomerated features in a Blob in many Blob(s)
export multiseparate # Multithreshold algorithm to separate agglomerated blobs
export curvature
export perimeter
export crosscorr
export area
export perimeterlength
export colorizedimage
export maskedimage
export perimeterlength
export ecd
export scorer
export interiorregions

include("rca.jl")
export rca
export metrics

include("watershed.jl")
export soille_watershed

include("particlesegregation.jl")
export Watershed

include("diluvian.jl")
export DiluvianCluster # Constructs the clusters from the data
export asmask # As a binary mask
export membercounts # Number of items in each cluster
export summarizeclusters # Statistics summarizing the clusters as a DataFrame
export clusterstats # Statistics summarizing one or more clusters
export clusters # The raw cluster assignments
export asimage
export defaultpalette

include("align.jl")
export rough_align
export AlignIntermediary
export offset, score, findbest, align

include("translate.jl")
export XYPosition
export translate

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
end


end
