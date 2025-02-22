using Compose
using Colors
using DataFrames
using Random: shuffle!

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
    cols::AbstractArray,
    clscol::Any = missing;
    title = missing,
    palette = TernPalette,
    font = "Verdana",
    fontsize = 5.0,
    norm = 1.0,
    withcount = true,
)
    # Draw an individual ternary diagram for cola and colb in data
    function ternary(data, cola, colb, v1, L, θ)
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
        datapts, labels, range = [], [], ismissing(clscol) ? (1:1) : levels(data[:, clscol])
        ilvl = 1
        for lvl in range
            f(row) = ismissing(clscol) || (row[clscol] == lvl)
            pts = [cc(r[cola] / norm, r[colb] / norm, v3, v2, v1) for r in filter(f, eachrow(data))]
            if length(pts) > 0
                pa, pb, col = [pt[1] for pt in pts], [pt[2] for pt in pts], palette[ilvl]
                push!(
                    datapts,
                    (
                        context(),
                        Compose.circle(pa, pb, [L * 0.003]),
                        Compose.stroke(col),
                        Compose.fill("transparent"),
                        Compose.clip([v3, v2, v1]),
                    ),
                )
                if !ismissing(clscol)
                    push!(
                        labels,
                        (
                            context(),
                            Compose.text(0.02, 0.54 + ilvl * 0.03, withcount ? "$lvl - $(length(pts))" : "$lvl"),
                            Compose.fill(col),
                            Compose.font(font),
                            Compose.fontsize(0.8 * fontsize),
                        ),
                    )
                end
                ilvl = (ilvl % length(palette)) + 1
            end
        end
        shuffle!(datapts) # randomize datapoint order
        if length(labels) > 15
            last = labels[end]
            labels = labels[1:12]
            push!(
                labels,
                (last[1], Compose.text(0.02, 0.54 + 13 * 0.03, "…"), Compose.fill(colorant"black"), last[4], last[5]),
            )
            push!(
                labels,
                (last[1], Compose.text(0.02, 0.54 + 14 * 0.03, last[2].primitives[1].value), last[3], last[4], last[5]),
            )
        end
        return (
            context(),
            (context(), datapts...),
            (context(), Compose.polygon([v1, v2, v3]), Compose.stroke("black"), Compose.fill("transparent")),
            (context(), Compose.line(vx), Compose.stroke("darkgray")),
            (context(), Compose.polygon([v1, v2, v3]), Compose.stroke("transparent"), Compose.fill(TernBack)),
            (context(), labels...),
        )
    end
    # Label the ends of the axes
    function label(v0, labels, L)
        xl, yl, ll = Float64[], Float64[], String[]
        function ss(l)
            if startswith(l,"Element(")
                try
                    parse(Element,l[9:end-1]).symbol
                catch
                    titlecase(l)
                end
            else
                titlecase(l)
            end
        end
        for i in eachindex(labels)
            θ = -π + π * i / 3
            vv = v(v0, θ, L)
            push!(xl, vv[1])
            push!(yl, vv[2])
            push!(ll, ss(labels[i]))
        end
        return (
            context(),
            Compose.text(xl, yl, ll, [hcenter], [vcenter]),
            Compose.fill("black"),
            Compose.font(font),
            Compose.fontsize(fontsize),
        )
    end
    v(v0, θ, off = 0.02) = (v0[1] + off * cos(θ), v0[2] - off * sin(θ))
    center, ops, ncols = (0.5, 0.5), [], min(length(cols), 6)
    for i = 2:ncols
        sa, sb = cols[i-1], cols[i]
        θ = -π / 2 + π / 3 * (i - 2)
        append!(ops, ternary(data, sa, sb, v(center, θ), 0.42, θ))
    end
    if !ismissing(title)
        append!(
            ops,
            (
                context(),
                Compose.text(0.5, 0.03, title, hcenter, vtop),
                Compose.fill("black"),
                Compose.font(font),
                Compose.fontsize(1.3 * fontsize),
            ),
        )
    end
    Compose.compose(context(0.05, 0.0, 0.9, 0.9), (context(), ops...), (context(), label(center, string.(cols[1:ncols]), 0.48)))
end

"""
    multiternary(
        zep::Zeppelin;
        rows = missing,
        omit = [ n"C", n"O" ],
        palette = TernPalette,
        fontsize = 12pt,
        font = "Verdana",
    )

Plot the elemental data in `zep` to a multi-ternary diagram.
"""
function multiternary(
        zep::Zeppelin;
        rows = missing,
        omit = [ n"C", n"O" ],
        palette = TernPalette,
        fontsize = 12pt,
        font = "Verdana",
)
    # Determine which elements to plot...
    zd = ismissing(rows) ? zep.data : zep.data[rows, :]
    df=DataFrame(Elm=Symbol[], Mean=Float64[])
    foreach(elm->push!(df, [ convert(Symbol,elm), mean(zd[:,convert(Symbol, elm)])]),  filter(elm->!(elm in omit), NeXLCore.elms(zep)))
    sort!(df, :Mean, rev=true)
    elms = df[:,:Elm][1:min(size(df,1),6)]
    multiternary(zd, elms, :CLASS, title=zep.header["DESCRIPTION"], palette=palette, norm=100.0, fontsize=fontsize, font=font)
end