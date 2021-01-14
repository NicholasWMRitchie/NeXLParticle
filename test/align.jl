using NeXLParticle
using DataFrames

const align_example = false
const align_example2 = true
const hist_test = false

function subsample(df, center, θ, dims)
    df2 = DataFrame(x=Float64[], y=Float64[], s=Int[])
    rot = [ cos(θ) -sin(θ); sin(θ) cos(θ) ]
    for r in eachrow(df)
        pp = rot * [ r[:x] - center[1], r[:y] - center[2] ]
        if abs(pp[1]) < dims[1] && abs(pp[2]) < dims[2]
            push!(df2, [pp..., r[:s]])
        end
    end
    df2
end

if align_example
    using Cairo
    using Fontconfig 
    using Gadfly
    println("Starting")
    const align_slow = false
    n, nθ, nrings = 10000, 180, 8
    # Construct a randomized "particle position" dataset
    df = DataFrame(x=100.0 .* (0.5 .- rand(n)), y=100.0 .* (0.5 .- rand(n)), s=[1 for _ in 1:n ])
    # Construct a localized subset that is offset and rotate 
    center, θ = 50.0 .* (0.5 .- rand(2)), 2π * rand()
    df2 = subsample(df, center, θ, (5.0, 5.0))

    ai = rough_align(df, df2, :x, :y, nrings, nθ)
    ai = @time rough_align(df, df2, :x, :y, nrings, nθ)

    print("Fast: ")
    println([ score(ai, bb) for bb in findbest(ai) ])
    println("$center vs. $([ offset(ai, bb) for bb in findbest(ai) ])")
    println("$(rad2deg(θ)) vs. $([ rad2deg(angle(ai, bb)) for bb in findbest(ai) ])")

    spy(ai.scores) |> PNG(joindir(homedir(),"Desktop","example.png"), 10inch, 8inch)
    if align_slow
        ai2 = @time rough_align_slow(df, df2, :x, :y, nrings, nθ)
        print("Slow: ")
        println([ score(ai2, bb) for bb in findbest(ai2) ])
        println("$center vs. $([ coordinate(ai2, bb) for bb in findbest(ai2) ])")
        println("$(rad2deg(θ)) vs. $([ rad2deg(angle(ai2, bb)) for bb in findbest(ai2) ])")
        spy(ai.scores) |> PNG(joinpath(homedir(),"Desktop","example_slow.png"), 10inch, 8inch)
    end
end

if align_example2
    using Cairo
    using Fontconfig 
    using Gadfly
    println("Starting")
    const align_slow = false
    n, nθ, nrings = 10000, 180, 8
    # Construct a randomized "particle position" dataset
    df = DataFrame(x=100.0 .* (0.5 .- rand(n)), y=100.0 .* (0.5 .- rand(n)), s=[1 for _ in 1:n ])
    # Construct a localized subset that is offset and rotate 
    center, θ = 50.0 .* (0.5 .- rand(2)), 2π * rand()
    df2 = subsample(df, center, θ, (35.0, 35.0))

    ai = rough_align(df, df2, :x, :y, nrings, nθ)
    ai = @time rough_align(df, df2, :x, :y, nrings, nθ)

    print("Fast: ")
    println([ score(ai, bb) for bb in findbest(ai) ])
    println("$center vs. $([ offset(ai, bb) for bb in findbest(ai) ])")
    println("$(rad2deg(θ)) vs. $([ rad2deg(angle(ai, bb)) for bb in findbest(ai) ])")

    spy(ai.scores) |> PNG(joinpath(homedir(),"Desktop","example.png"), 10inch, 8inch)
    plot(ai, df, df2) |> SVG(joinpath(homedir(), "Desktop", "overlaid.svg"))
    if align_slow
        ai2 = @time rough_align_slow(df, df2, :x, :y, nrings, nθ)
        print("Slow: ")
        println([ score(ai2, bb) for bb in findbest(ai2) ])
        println("$center vs. $([ coordinate(ai2, bb) for bb in findbest(ai2) ])")
        println("$(rad2deg(θ)) vs. $([ rad2deg(angle(ai2, bb)) for bb in findbest(ai2) ])")
        spy(ai.scores) |> PNG(joinpath(homedir(),"Desktop","example_slow.png"), 10inch, 8inch)
    end
end

if hist_test
    using Gadfly, CSV
    n = 10000
    nθ, nrings = 32, 8
    res = DataFrame(Iteration=Int[], Index=Int[], Score=Float64[], Ratio=Float64[], xOff=Float64[], yOff=Float64[], θOff= Float64[])
    pd = DataFrame(xOff=Float64[], yOff=Float64[], θOff=Float64[], score=Float64[])
    for i in 1:100
        # Construct a randomized "particle position" dataset
        df = DataFrame(x=100.0 .* (0.5 .- rand(n)), y=100.0 .* (0.5 .- rand(n)), s=[1 for _ in 1:n ])
        # Construct a localized subset that is offset and rotate 
        center, θ = 50.0 .* (0.5 .- rand(2)), 2π * rand()
        df2 = subsample(df, center, θ, (5.0, 5.0))

        ai = @time rough_align(df, df2, :x, :y, nrings, nθ)
        b = findbest(ai)
        for j in eachindex(b)
            off = offset(ai, b[j])
            push!(res, ( i, j, score(ai,b[j]), score(ai,b[j])/score(ai,b[1]), off[1]-center[1], off[2]-center[2], nθ*(modf((θ - angle(ai,b[j]) + 5π)/(2π))[1]-0.5)))
        end
        off = offset(ai, b[1])
        push!(pd, (off[1]-center[1], off[2]-center[2], nθ*(modf((θ - angle(ai,b[1]) + 5π)/(2π))[1]-0.5), score(ai,b[1])))
    end
    plot(
        layer(x=filter(x->abs(x)<1.0, pd[:,:xOff]), Geom.histogram(bincount=10), Theme(default_color="blue")),
        layer(x=filter(y->abs(y)<1.0, pd[:,:yOff]), Geom.histogram(bincount=10), Theme(default_color="red"))
    ) |> SVG(joinpath(homedir(),"Desktop","xy.svg"))
    plot(x=filter(x->abs(x)<1.0, pd[:,:θOff]), Geom.histogram(bincount=10), Theme(default_color="green")) |> SVG(joinpath(homedir(),"Desktop","off.svg"))
    plot(x=pd[:,:score], Geom.histogram(bincount=10), Theme(default_color="green")) |> SVG(joinpath(homedir(),"Desktop","score.svg"))
    CSV.write(joinpath(homedir(),"Desktop","align_hist.csv"), res)
end

