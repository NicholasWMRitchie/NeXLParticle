using Compose
using Colors
using DataFrames

const TernBack = RGB(253 / 255, 253 / 255, 241 / 255)

const TernPalette = distinguishable_colors(66, Color[TernBack, RGB(0, 0, 0), RGB(0 / 255, 168 / 255, 45 / 255)])[3:end]
const TernColorblind =
    distinguishable_colors(66, Color[TernBack, RGB(0, 0, 0), colorant"DodgerBlue4"], transform = deuteranopic)[3:end]

"""
    multiternary(data::DataFrame, cols::AbstractArray{Symbol}, clscol::Union{Symbol, Missing}=missing; palette=TernPalette, norm=1.0)

Uses Compose.jl to draw a simple multi-ternary diagram from the `data`, a `DataFrame`.
The columns are drawn in the order in `cols` and the vertices labeled with the names
of the columns.
"""
function multiternary(
    data::DataFrame,
    cols::AbstractArray{Symbol},
    clscol::Union{Symbol,Missing} = missing;
    title = missing,
    palette = TernPalette,
    deffont = "Verdana",
    fontsz = 12Compose.pt,
    norm = 1.0,
)
    # Draw an individual ternary diagram for cola and colb in data
    function ternary(ctx, data, cola, colb, v1, L, θ)
        # Compute the intersection of the line (P1:P2) with (P3:P4)
        function xy(P1, P2, P3, P4)
            s =
                (P1[2] * (P2[1] - P3[1]) + P2[2] * P3[1] - P2[1] * P3[2] + P1[1] * (P3[2] - P2[2])) /
                (P2[2] * (P3[1] - P4[1]) + P1[2] * (P4[1] - P3[1]) + (P1[1] - P2[1]) * (P3[2] - P4[2]))
            return (P3[1] + s * (P4[1] - P3[1]), P3[2] + s * (P4[2] - P3[2]))
        end
        # From a, b on [0,1] compute P1, P2, P3 and P4
        function cc(a, b, va, vb, vc)
            return xy(vb .+ a .* (va .- vb), vc .+ a .* (va .- vc), va .+ b .* (vb .- va), vc .+ b .* (vb .- vc))
        end
        v3 = (v1[1] + L * cos(θ - π / 6), v1[2] - L * sin(θ - π / 6))
        v2 = (v1[1] + L * cos(θ + π / 6), v1[2] - L * sin(θ + π / 6))
        # 0.1 increment lines
        vx = [[v1 .+ f .* (v2 .- v1), v1 .+ f .* (v3 .- v1)] for f = 0.1:0.1:0.91]
        append!(vx, [[v2 .+ f .* (v3 .- v2), v2 .+ f .* (v1 .- v2)] for f = 0.1:0.1:0.91])
        append!(vx, [[v3 .+ f .* (v1 .- v3), v3 .+ f .* (v2 .- v3)] for f = 0.1:0.1:0.91])
        # Data points
        datapts, range = [], ismissing(clscol) ? (1:1) : levels(data[:, clscol])
        ilvl = 1
        for lvl in range
            f(row) = ismissing(clscol) || (row[clscol] == lvl)
            pts = [cc(r[cola] / norm, r[colb] / norm, v3, v2, v1) for r in filter(f, eachrow(data))]
            if length(pts) > 0
                pa, pb, col = [pt[1] for pt in pts], [pt[2] for pt in pts], palette[ilvl]
                push!(
                    datapts,
                    (
                        ctx,
                        Compose.circle(pa, pb, [L * 0.01]),
                        Compose.stroke(col),
                        Compose.fill("transparent"),
                        Compose.clip([v3, v2, v1]),
                    ),
                )
                if !ismissing(clscol)
                    push!(
                        datapts,
                        (
                            ctx,
                            Compose.text(0.02, 0.54 + ilvl * 0.03, "$lvl - $(length(pts))"),
                            Compose.stroke(col),
                            Compose.font(deffont),
                            Compose.fontsize(fontsz),
                        ),
                    )
                end
                ilvl = (ilvl % length(palette)) + 1
            end
        end
        return (
            ctx,
            datapts...,
            (ctx, Compose.polygon([v1, v2, v3]), Compose.stroke("black"), Compose.fill("transparent")),
            (ctx, Compose.line(vx), Compose.stroke("darkgray")),
            (ctx, Compose.polygon([v1, v2, v3]), Compose.stroke("transparent"), Compose.fill(TernBack)),
        )
    end
    # Label the ends of the axes
    function label(ctx, v0, labels, L)
        xl, yl, ll = Float64[], Float64[], String[]
        for i in eachindex(labels)
            θ = -π + π * i / 3
            vv = v(v0, θ, L)
            push!(xl, vv[1])
            push!(yl, vv[2])
            push!(ll, titlecase(labels[i]))
        end
        return (
            ctx,
            Compose.text(xl, yl, ll, [hcenter], [vcenter]),
            Compose.stroke("black"),
            Compose.font(deffont),
            Compose.fontsize(fontsz),
        )
    end
    v(v0, θ, off = 0.02) = (v0[1] + off * cos(θ), v0[2] - off * sin(θ))
    center = (0.5, 0.5)
    ops = []
    for i = 2:length(cols)
        sa, sb = cols[i-1], cols[i]
        θ = -π / 2 + π / 3 * (i - 2)
        append!(ops, ternary(context(), data, sa, sb, v(center, θ), 0.42, θ))
    end
    if !ismissing(title)
        append!(
            ops,
            (
                Compose.text(0.5, 0.03, title, hcenter, vtop),
                Compose.stroke("black"),
                Compose.font(deffont),
                Compose.fontsize(1.3*fontsz),
            ),
        )
    end
    compose(context(), label(context(), center, string.(cols), 0.48), ops...)
end
