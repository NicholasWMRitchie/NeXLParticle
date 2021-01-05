

using Test
using NeXLParticle
using StaticArrays
using CoordinateTransformations

@testset "Translate" begin
    oldPts = [ SA[ 10.0 * rand(Float64,2)...] for _ in 1:20 ]
    θ, sc, off = π/3, SA[ 1.1 0.0; 0.0 0.9], SA[10.0,-10.0]
    tr = AffineMap( SA[ cos(θ) -sin(θ) ; sin(θ) cos(θ) ] * sc, off)
    newPts = (x->SA[ (tr(x) .+ 0.1*rand(2))...]).(oldPts)

    trp = translate(oldPts, newPts)
    diff = inv(trp).(newPts) .- oldPts
    for d in diff
        @test abs(d[1])<0.1
        @test abs(d[2])<0.1
    end
    diff = trp.(oldPts) .- newPts
    for d in diff
        @test abs(d[1])<0.1
        @test abs(d[2])<0.1
    end
end