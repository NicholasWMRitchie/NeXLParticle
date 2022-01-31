using Test
using NeXLParticle
using Random
using StaticArrays
using CoordinateTransformations, Rotations


@testset "ICP32" begin
    r, npts, ty = MersenneTwister(0xBADFEED), 100, Float32
    xy = SVector.(zip(rand(r, ty, npts), rand(r, ty, npts)))
    @test eltype(eltype(xy)) == ty
    # Muddle the points up
    xyp = AffineMap(Angle2d{ty}(deg2rad(1.2f0)), [0.1f0, -0.1f0]).(xy) .+ 0.02f0 * SVector.(collect(zip(0.5f0 .- rand(r, ty, npts), 0.5f0 .- rand(r, ty, npts))))
    xyps = shuffle(r, xyp)
    @test eltype(eltype(xyps)) == ty
    @test isapprox(NeXLParticle.icperror(xy, xyps), 5.79987f0, atol=0.001f0)
    @test isapprox(NeXLParticle.icperror(xy[1:2:100], xyps),7.61191f0,atol=0.001f0)
    @test isapprox(NeXLParticle.icperror(xy, xyps[1:2:100]), 2.85196f0, atol=0.001f0)
    res=NeXLParticle.icp(xy, xyps)
    @test isapprox(NeXLParticle.icperror(xy, res), 0.13382f0, atol=0.001f0)
end

@testset "ICP64_1" begin
    r, npts, ty = MersenneTwister(0xBADF00D), 100, Float64
    xy = SVector.(zip(rand(r, ty, npts), rand(r, ty, npts)))
    @test eltype(eltype(xy)) == ty
    # Muddle the points up
    xyp = AffineMap(Angle2d{ty}(deg2rad(1.2f0)), [0.1f0, -0.1f0]).(xy) .+ 0.02f0 * SVector.(collect(zip(0.5f0 .- rand(r, ty, npts), 0.5f0 .- rand(r, ty, npts))))
    xyps = shuffle(r, xyp)[1:npts÷2] # Pick only fifty points to match
    @test eltype(eltype(xyps)) == ty
    @test isapprox(NeXLParticle.icperror(xy, xyps), 2.5320, atol=0.0001)
    @test isapprox(NeXLParticle.icperror(xy[1:2:100], xyps), 3.032293, atol=0.0001)
    res=icp(xy, xyps)
    @test isapprox(NeXLParticle.icperror(xy, res), 0.075434, atol=0.0001)
end

@testset "ICP64_2" begin
    r, npts, ty = MersenneTwister(0xBADADDA), 100, Float64
    xy = SVector.(zip(rand(r, ty, npts), rand(r, ty, npts)))
    @test eltype(eltype(xy)) == ty
    # Muddle the points up
    xyp = AffineMap(Angle2d{ty}(deg2rad(1.2f0)), [0.1f0, -0.1f0]).(xy) .+ 0.02f0 * SVector.(collect(zip(0.5f0 .- rand(r, ty, npts), 0.5f0 .- rand(r, ty, npts))))
    xyps = shuffle(r, xyp)[1:npts÷2] # Pick only fifty points to match
    @test eltype(eltype(xyps)) == ty
    @test isapprox(NeXLParticle.icperror(xy, xyps), 3.458189, atol=0.0001)
    @test isapprox(NeXLParticle.icperror(xy[1:2:100], xyps), 4.235188, atol=0.0001)
    res=icp(xy, xyps)
    @test isapprox(NeXLParticle.icperror(xy, res), 0.172820, atol=0.0001)
end