using Test
using NeXLParticle
using StaticArrays
using CoordinateTransformations
using Random

@testset "Translate" begin
    r = MersenneTwister(0xBADF00D)
    oldPts = [ SA[ 10.0 * rand(r, Float64,2)...] for _ in 1:20 ]
    θ, sc, off = π/3, SA[ 1.1 0.0; 0.0 0.9], SA[10.0,-10.0]
    tr = AffineMap( SA[ cos(θ) -sin(θ) ; sin(θ) cos(θ) ] * sc, off)
    newPts = (x->SA[ (tr(x) .+ 0.1*rand(r, 2))...]).(oldPts)
    # Estimate the transform that takes `oldPts` to `newPts`
    # tr and tr′ should be very similar
    tr′ = translate(oldPts, newPts)
    diff = tr′.(oldPts) .- newPts
    @test all(d->abs(d[1])<0.1, diff)
    @test all(d->abs(d[2])<0.1, diff)
    # Check that inverting the transform brings `newPts` back to `oldPts`
    itr′ = inv(tr′)
    diff = itr′.(newPts) .- oldPts
    @test all(d->abs(d[1])<0.1, diff)
    @test all(d->abs(d[2])<0.1, diff)

    zep = Zeppelin(joinpath(datadep"ZepTestArtifact","test.hdz"))
    zeptr = translate(zep, tr)
    # Check that `tr` and `tr′` both produce a similar forward transformation
    zeptr′ = translate(zep, tr′)
    @test all(rr->abs(rr[1].XABS-rr[2].XABS)<0.1, zip(eachrow(zeptr.data), eachrow(zeptr′.data)))
    @test all(rr->abs(rr[1].YABS-rr[2].YABS)<0.1, zip(eachrow(zeptr.data), eachrow(zeptr′.data)))

    # Check that `inv(tr′)` brings us back to close to `zep`
    izeptr′ = translate(zeptr, inv(tr′))
    @test all(rr->abs(rr[1].XABS-rr[2].XABS)<0.1, zip(eachrow(zep.data), eachrow(izeptr′.data)))
    @test all(rr->abs(rr[1].YABS-rr[2].YABS)<0.1, zip(eachrow(zep.data), eachrow(izeptr′.data)))
end