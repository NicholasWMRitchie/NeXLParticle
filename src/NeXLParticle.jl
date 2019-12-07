module NeXLParticle

using CSV
using DataFrames
using PeriodicTable
using DataStructures

include("zepellin.jl")
export Zepellin # struct holding a hdz/pxz file pair
export classes # An ordered list of particle classes
export elements # An ordered list of elements in quant table
export header # A alphabetically sorted list of header items
export data # Provides access to the DataFrame containing the particle data

include("aspextiff.jl")
# Exposes FileIO.load and FileIO.save for the format"ASPEX TIFF"

end
