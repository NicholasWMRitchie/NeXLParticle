module NeXLParticle

using DataFrames
using FileIO
using Reexport
using Requires

@reexport using PeriodicTable
@reexport using NeXLCore
@reexport using NeXLSpectrum


include("zeppelin.jl")
export Zeppelin # struct holding a hdz/pxz file pair
export classes # An ordered list of particle classes
export header # A alphabetically sorted list of header items
export data # Provides access to the DataFrame containing the particle data
export iszeppelin # Is this IOStream or filename a Zeppelin file?
export spectrum # read a spectrum associate with a Zeppelin file (or missing)
export eachparticle # The range of partice row indices
export rowsMax, rowsClass
export MORPH_COLS, COMP_COLS, CLASS_COLS
export allelms
export beamenergy, probecurrent
export maxparticle

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
export BaseRules
export NullRules

function __init__()
    @require Gadfly = "c91e804a-d5a3-530f-b6f0-dfbca275c004" include("gadflysupport.jl")
    @info "Loading Gadfly support into NeXLParticle."
end


end
