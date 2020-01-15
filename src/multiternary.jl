using Compose
using Colors
using DataFrames

const TernBack = RGB(253 / 255, 253 / 255, 241 / 255)

const TernPalette = distinguishable_colors(
    66,
    Color[TernBack, RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)],
)[3:end]
const TernColorblind = distinguishable_colors(
    66,
    Color[TernBack, RGB(0, 0, 0), colorant"DodgerBlue4"],
    transform = deuteranopic,
)[3:end]

"""
    multiternary(data::DataFrame, cols::AbstractArray{Symbol}, clscol::Union{Symbol, Missing}=missing)

Uses Compose.jl to draw a simple multi-ternary diagram from the `data`, a `DataFrame`.
The columns are drawn in the order in `cols` and the vertices labeled with the names
of the columns.
"""
function multiternary(
        data::DataFrame,
        cols::AbstractArray{Symbol},
        clscol::Union{Symbol, Missing}=missing;
        palette = TernPalette)
    # Draw an individual ternary diagram for cola and colb in data
    function ternary(ctx, data, cola, colb, v1, L, θ)
        # Compute the intersection of the line (P1:P2) with (P3:P4)
        function xy(P1,P2,P3,P4)
            s=(P1[2]*(P2[1] - P3[1]) + P2[2]*P3[1] - P2[1]*P3[2] + P1[1]*(P3[2] - P2[2]))/
            (P2[2]*(P3[1] - P4[1]) + P1[2]*(P4[1] - P3[1]) + (P1[1] - P2[1])*(P3[2] - P4[2]))
            return ( P3[1] + s*(P4[1] - P3[1]), P3[2] + s*(P4[2] - P3[2]) )
        end
        # From a, b on [0,1] compute P1, P2, P3 and P4
        function cc(a, b, va, vb, vc)
            return xy(
                vb .+ a .* ( va .- vb),
                vc .+ a .* ( va .- vc),
                va .+ b .* ( vb .- va),
                vc .+ b .* ( vb .- vc) )
        end
        v3 = ( v1[1]+L*cos(θ-π/6), v1[2]-L*sin(θ-π/6) )
        v2 = ( v1[1]+L*cos(θ+π/6), v1[2]-L*sin(θ+π/6) )
        # 0.1 increment lines
        vx = [  [ v1 .+ f.*(v2 .- v1), v1 .+ f.*(v3 .- v1) ] for f in 0.1:0.1:0.91 ]
        append!(vx, [ [ v2 .+ f.*(v3 .- v2), v2 .+ f.*(v1 .- v2) ] for f in 0.1:0.1:0.91 ])
        append!(vx, [ [ v3 .+ f.*(v1 .- v3), v3 .+ f.*(v2 .- v3) ] for f in 0.1:0.1:0.91 ])
        # Data points
        datapts, range = [], ismissing(clscol) ? (1:1) : levels(data[:,clscol])
        for (ilvl, lvl) in enumerate(range)
            f(row) = ismissing(clscol) || (row[clscol] == lvl)
            pts = [ cc(r[cola], r[colb], v3, v2, v1) for r in filter(f, eachrow(data)) ]
            pa, pb = [ pts[i][1] for i in eachindex(pts)], [ pts[i][2] for i in eachindex(pts)]
            push!(datapts, ( ctx, circle(pa, pb, [L*0.01]), stroke(palette[ilvl]), clip([v1,v2,v3]) ))
            push!(datapts, (ctx, text(0.2, 0.5+ilvl*0.05, lvl), stroke(palette[ilvl])))
        end
        return ( ctx,
            datapts...,
            ( ctx, polygon([v1, v2, v3]), stroke("black"), fill("transparent") ) ,
            ( ctx, line(vx), stroke("dark gray") ) ,
            ( ctx, polygon([v1, v2, v3]), stroke("transparent"), fill(TernBack) ) ,
        )
    end
    # Label the ends of the axes
    function label(ctx, v0, labels, L)
        xl, yl, ll = Float64[], Float64[], String[]
        for i in eachindex(labels)
            θ = -π+π*i/3
            vv = v(v0, θ, L)
            push!(xl, vv[1])
            push!(yl, vv[2])
            push!(ll, labels[i])
        end
        return (ctx, text(xl, yl, ll, [hcenter], [vcenter]), stroke("black") )
    end
    v(v0, θ, off=0.02)=(v0[1] + off*cos(θ), v0[2]-off*sin(θ))
    center = (0.5, 0.5)
    ops = []
    for i in 2:length(cols)
        sa, sb = cols[1], cols[2]
        θ = -π/2+π/3*(i-2)
        append!(ops, ternary(context(), data, sa, sb, v(center,θ), 0.42, θ ) )
    end
    compose(context(), label(context(), center, string.(cols), 0.48), ops...)
end
