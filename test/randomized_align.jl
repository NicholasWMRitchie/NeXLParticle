using NeXLParticle
using StaticArrays
using GeometryBasics
using Random
using LinearAlgebra
using Rotations
using Gadfly

r0 = MersenneTwister(0xEA7F00D)
for npts in (100, 1000, 10000, 100000)
    println("Working on $npts points...")
    for i in 1:10
        seed = rand(r0, UInt)
        rnd = MersenneTwister(seed)
        snp = sqrt(npts)
        bounds = Rect2((-0.5*snp, -0.5*snp), (snp,snp))
        xy = NeXLParticle.generate_ground_truth(rnd, bounds, npts, x->norm(x)<0.5*snp)
        θ1, θ2 = 2π*rand(rnd, 2)
        f1, f2 = 0.5 .+ 0.5 * rand(rnd, 2)
        ps1 = NeXLParticle.measure_particles(rnd, xy, SA[100.0*rand(rnd),100.0*rand(rnd)], θ1, f1, 0.001)
        ps2 = NeXLParticle.measure_particles(rnd, xy, SA[100.0*rand(rnd),100.0*rand(rnd)], θ2, f2, 0.001)

        (ct1, ct2) = @time rough_align(ps1, ps2, tol=0.001)

        r = RotMatrix{2}(θ1-θ2)
        if !isapprox(r, ct2.linear, atol = 0.01)
            
            7println("* * * * rough_align(...) failure * * * *")
            println(" Failure:  $seed")
            println(" NPoints: $npts")
            println("Fraction:f1 = $f1, f2 = $f2")
            println("   Angle: θ1 = $θ1, θ2 = $θ2")
            println("   Truth: $r")
            println("Estimate: $ct2.linear")
            ts1, ts2 = ct1.(ps1), ct2.(ps2)
            p = plot( ( layer(x=map(p->p[1], t), y=map(p->p[2], t), Geom.point, Theme(default_color=c)) for (t,c) in zip((ts1, ts2), ("red", "blue"))) ...)
            p |> SVG(joinpath(homedir(),"Desktop","RA[$seed].svg"), 20inch, 20inch)
        end
        GC.gc()
    end
end
#=  
The algorithm seems to fail when:
  1) the number of datapoints is small and the measurement fraction is low
  2) density of particles is too high
  3) variability in measured positions is high

  Offset and angles don't matter.

Timings:
Working on 100 points... (Average: 4.3 ms)
  0.002455 seconds (34.60 k allocations: 1.318 MiB)
  0.002719 seconds (56.89 k allocations: 2.134 MiB)
  0.003065 seconds (63.53 k allocations: 2.380 MiB)
  0.001486 seconds (29.57 k allocations: 1.121 MiB)
  0.006415 seconds (134.60 k allocations: 5.160 MiB)
  0.004236 seconds (67.86 k allocations: 2.573 MiB)
  0.004809 seconds (104.54 k allocations: 3.893 MiB)
  0.004412 seconds (51.79 k allocations: 1.950 MiB)
  0.006890 seconds (139.50 k allocations: 5.320 MiB)
  0.006069 seconds (113.35 k allocations: 4.214 MiB)
Working on 1,000 points...  (Average: 9.99 ms)
  0.010888 seconds (185.80 k allocations: 7.901 MiB)
  0.013183 seconds (183.53 k allocations: 7.804 MiB)
  0.011232 seconds (182.69 k allocations: 7.709 MiB)
  0.008922 seconds (180.93 k allocations: 7.608 MiB)
  0.010046 seconds (182.18 k allocations: 7.679 MiB)
  0.006969 seconds (134.49 k allocations: 5.283 MiB)
  0.010058 seconds (185.00 k allocations: 7.854 MiB)
  0.008368 seconds (158.51 k allocations: 6.278 MiB)
  0.010471 seconds (185.53 k allocations: 7.887 MiB)
  0.009765 seconds (181.95 k allocations: 7.672 MiB)
Working on 10,000 points...  (Average: 29.5 ms)
  0.037703 seconds (360.80 k allocations: 18.659 MiB)
  0.027629 seconds (305.02 k allocations: 15.226 MiB)
  0.027257 seconds (291.41 k allocations: 14.402 MiB)
  0.027362 seconds (310.73 k allocations: 15.577 MiB)
  0.025490 seconds (300.53 k allocations: 14.935 MiB)
  0.026743 seconds (297.59 k allocations: 14.743 MiB)
  0.034434 seconds (353.36 k allocations: 18.208 MiB)
  0.034588 seconds (354.96 k allocations: 18.272 MiB)
  0.031228 seconds (320.44 k allocations: 16.160 MiB)
  0.022824 seconds (269.03 k allocations: 13.118 MiB)
Working on 100,000 points...  (Average:  285.3 ms)
  0.270818 seconds (1.52 M allocations: 89.571 MiB, 6.31% gc time)
  0.349587 seconds (1.98 M allocations: 118.098 MiB, 5.10% gc time)
  0.231989 seconds (1.49 M allocations: 87.577 MiB)
  0.335018 seconds (1.89 M allocations: 112.289 MiB, 5.34% gc time)
  0.256628 seconds (1.50 M allocations: 88.605 MiB, 7.25% gc time)
  0.283030 seconds (1.64 M allocations: 96.868 MiB, 6.65% gc time)
  0.298009 seconds (1.64 M allocations: 96.943 MiB, 7.23% gc time)
  0.348799 seconds (2.01 M allocations: 120.019 MiB, 5.38% gc time)
  0.221129 seconds (1.37 M allocations: 80.252 MiB)
  0.257807 seconds (1.51 M allocations: 89.376 MiB, 6.76% gc time)
=#