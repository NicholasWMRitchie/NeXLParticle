using NeXLParticle
using Test

@testset "Blob" begin
    spec=loadspectrum(ASPEXTIFF, joinpath(@__DIR__, "images", "test1.tif"), withImgs=true)
    bse = spec[:Image1]

    bs = blob(bse, i->i>0.7)
    @test area.(bs) == [ 2882, 43, 1, 1, 1, 1]
    @test perimeter(bs[1])[24]==CartesianIndex(66,12)
    @test length(perimeter(bs[1]))==240
    @test size(CartesianIndices(bs[1]))==(69, 66)
    @test isapprox(perimeterlength(bs[1]),292.19,atol=0.01)

    @test isapprox(ecd(bs[1],true), 60.81, atol=0.01)
    @test isapprox(ecd(bs[1]), 60.81, atol=0.01)
    @test isapprox(ecd(bs[1],false), 60.58, atol=0.01)

    @test length(curvature(bs[1], 8))==length(perimeter(bs[1]))
    bs2=separate(bs[1])
    @test length(bs2)==10
    @test area.(bs2) == [ 1415, 836, 183, 227, 57, 2, 3, 3, 1, 1 ]
    @test_nowarn colorizedimage(bs, bse)
    @test_nowarn colorizedimage(bs2, bse)

    ms=multiseparate(bse, 0.6:0.05:1.0)

    ir=interiorregions(bs[1])
    @test area.(ir)== [15,2,1,1,1,1,1,1]
    @test NeXLParticle.filledarea(bs[1]) > area(bs[1])
    @test NeXLParticle.filledarea(bs[1]) == area(bs[1]) + sum(area.(ir))
    @test isapprox(sum(maskedimage(bs[1], bse), init=1.0), 2798.8, atol=0.1)
end