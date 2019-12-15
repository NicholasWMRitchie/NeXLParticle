module NeXLParticle

using CSV
using DataFrames
using PeriodicTable
using DataStructures
using Reexport

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

include("signature.jl")
export signature # Compures the particle signature

end
