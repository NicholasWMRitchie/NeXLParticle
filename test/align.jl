using Test
using NeXLParticle
using StaticArrays
using GeometryBasics
using Random
using LinearAlgebra
using Rotations
#using Gadfly

@testset "rough_align" begin
    @testset "rough_align (BADFOOD)" begin
        rnd = MersenneTwister(0xBADF00D)
        npts = 1000
        bounds = Rect2((-5.0, -1.5), (10.0,3.0))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[-10.0,5.0], deg2rad(12.0), 0.8, 0.003)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[-2.0,-1.0], deg2rad(213.0), 0.78, 0.003)

        (ct1, ct2) = NeXLParticle.rough_align(ps1, ps2, tol=0.001)
        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [9.463509135855892, -5.084006487177425], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.933518161752927 -0.35853011265085666; 0.35853011265085666 -0.933518161752927], atol=1.0e-8)
        @test isapprox(ct2.translation, [-2.7619203122160743, -0.2997257033588756], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(12.0-213.0)), ct2.linear, atol = 0.002)

        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = correspondences(ts1, ts2; tol=0.01)
        @test norm(ts1[c1[1]]-ts2[c2[1]]) < 0.01
        @test norm(ts1[c1[2]]-ts2[c2[2]]) < 0.01
        @test norm(ts1[c1[3]]-ts2[c2[3]]) < 0.01
        @test norm(ts1[c1[end]]-ts2[c2[end]]) < 0.01

        # plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
    end

    @testset "rough_align (0xABCDEF)" begin
        rnd = MersenneTwister(0xABCDEF)
        npts = 2000
        bounds = Rect2((-10.0, -10.0), (20.0,20.0))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[2.0,4.0], deg2rad(230.0), 0.99, 0.003)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[-2.0,-6.0], deg2rad(23.0), 0.99, 0.003)

        (ct1, ct2) = NeXLParticle.rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-2.439381678523621, -4.83348719078295], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.8909816833778803 0.4540392492782079; -0.4540392492782079 -0.8909816833778803], atol=1.0e-8)
        @test isapprox(ct2.translation, [0.5030208959078659, -7.088192398044117], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(230.0-23.0)), ct2.linear, atol = 0.001)


        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = NeXLParticle.correspondences(ts1, ts2; tol=0.01)
        @test norm(ts1[c1[1]]-ts2[c2[1]]) < 0.01
        @test norm(ts1[c1[2]]-ts2[c2[2]]) < 0.01
        @test norm(ts1[c1[3]]-ts2[c2[3]]) < 0.01
        @test norm(ts1[c1[end]]-ts2[c2[end]]) < 0.01

        # plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
    end

    @testset "rough_align (0xF156F00D)" begin
        rnd = MersenneTwister(0xF156F00D)
        npts = 100
        bounds = Rect2((-10.0, -10.0), (20.0,20.0))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[2.0,4.0], deg2rad(230.0), 0.89, 0.003)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[-2.0,-6.0], deg2rad(23.0), 0.94, 0.003)

        (ct1, ct2) = NeXLParticle.rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-3.600666188303255, -3.823165241126716], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.8910316678635093 0.4539411491199855; -0.4539411491199855 -0.8910316678635093], atol=1.0e-8)
        @test isapprox(ct2.translation, [-0.6588973878728955, -6.076703669619125], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(230.0-23.0)), ct2.linear, atol = 0.001)

        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = NeXLParticle.correspondences(ts1, ts2; tol=0.01)
        @test norm(ts1[c1[1]]-ts2[c2[1]]) < 0.01
        @test norm(ts1[c1[2]]-ts2[c2[2]]) < 0.01
        @test norm(ts1[c1[3]]-ts2[c2[3]]) < 0.01
        @test norm(ts1[c1[end]]-ts2[c2[end]]) < 0.01

        # plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
    end

    @testset "rough_align (0x7357)" begin
        rnd = MersenneTwister(0x7357)
        npts = 10000
        bounds = Rect2((-5.0, -5.0), (10.0,10.0))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[2.0,4.0], deg2rad(113.0), 0.7, 0.001)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[-2.0,-6.0], deg2rad(213.0), 0.65, 0.001)

        (ct1, ct2) = NeXLParticle.rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-1.9797624240598497, -5.2545418722493595], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.17353087223345756 0.9848284299216262; -0.9848284299216262 -0.17353087223345756], atol=1.0e-8)
        @test isapprox(ct2.translation, [5.582172374918754, -4.264693091414748], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(113.0-213.0)), ct2.linear, atol = 0.001)
        

        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = NeXLParticle.correspondences(ts1, ts2; tol=0.01)
        @test length(c1) == 4668 # approx 10000*0.7*0.67
        @test norm(ts1[c1[1]]-ts2[c2[1]]) < 0.01
        @test norm(ts1[c1[2]]-ts2[c2[2]]) < 0.01
        @test norm(ts1[c1[3]]-ts2[c2[3]]) < 0.01
        @test norm(ts1[c1[end]]-ts2[c2[end]]) < 0.01

        #plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
    end
end