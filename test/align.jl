using Test
using NeXLParticle
using StaticArrays
using GeometryBasics
using Random
using LinearAlgebra
using Rotation
#using Gadfly

@testset "rough_align" begin
    @testset "rough_align (BADFOOD)" begin
        rnd = MersenneTwister(0xBADF00D)
        npts = 1000
        bounds = Rect2((-5.0, -1.5), (10.0,3.0))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[-10.0,5.0], deg2rad(12.0), 0.8, 0.003)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[-2.0,-1.0], deg2rad(213.0), 0.78, 0.003)

        (ct1, ct2) = rough_align(ps1, ps2, tol=0.001)
        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [9.334994200384267, -5.198016372692697], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.9339031207175221 -0.35752616843256846; 0.35752616843256846 -0.9339031207175221], atol=1.0e-8)
        @test isapprox(ct2.translation, [-2.8904928250171156, -0.4162956186599809], atol=1.0e-8 )
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

        (ct1, ct2) = rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-1.038503062049208, -5.350538616445855], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.8909411399146158 0.45411880076434297; -0.45411880076434297 -0.8909411399146158], atol=1.0e-8)
        @test isapprox(ct2.translation, [1.904393007377276, -7.6051699023220465], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(230.0-23.0)), ct2.linear, atol = 0.001)


        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = correspondences(ts1, ts2; tol=0.01)
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

        (ct1, ct2) = rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-2.664472483448743, -4.768618139930582], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.8910676616790995 0.45387049067960117; -0.45387049067960117 -0.8910676616790995], atol=1.0e-8)
        @test isapprox(ct2.translation, [0.2768702232176849, -7.022161641453611], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(230.0-23.0)), ct2.linear, atol = 0.001)

        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = correspondences(ts1, ts2; tol=0.01)
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

        (ct1, ct2) = rough_align(ps1, ps2, tol=0.001)

        @test isapprox(ct1.linear, LinearAlgebra.I, atol=1.0e-8)
        @test isapprox(ct1.translation, [-3.570780539551327, -6.321500211234103], atol=1.0e-8)
        @test isapprox(ct2.linear, [-0.1735628929293961 0.9848227872048753; -0.9848227872048753 -0.1735628929293961], atol=1.0e-8)
        @test isapprox(ct2.translation, [3.990875748179417, -5.3317038874444185], atol=1.0e-8 )
        @test isapprox(RotMatrix{2}(deg2rad(113.0-213.0)), ct2.linear, atol = 0.001)
        

        ts1, ts2 = ct1.(ps1), ct2.(ps2)
        c1, c2 = correspondences(ts1, ts2; tol=0.01)
        @test length(c1) == 4667 # approx 10000*0.7*0.67
        @test norm(ts1[c1[1]]-ts2[c2[1]]) < 0.01
        @test norm(ts1[c1[2]]-ts2[c2[2]]) < 0.01
        @test norm(ts1[c1[3]]-ts2[c2[3]]) < 0.01
        @test norm(ts1[c1[end]]-ts2[c2[end]]) < 0.01

        #plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
    end
end