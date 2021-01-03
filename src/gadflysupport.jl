using .Gadfly
using Colors

function Gadfly.plot(
    zep::Zeppelin,
    rows;
    autoklms = true,
    klms = [],
    xmin = 0.0,
    xmax = missing,
    norm = :None,
    yscale = 1.05,
    ytransform = identity,
    style = NeXLSpectrum.NeXLSpectrumStyle,
    palette = NeXLCore.NeXLPalette,
)
    plot(
        Spectrum[filter(s->!ismissing(s), map(r -> spectrum(zep, r, false),rows))...],
        autoklms = autoklms,
        klms = klms,
        xmin = xmin,
        xmax = xmax,
        norm = norm,
        yscale = yscale,
        ytransform = ytransform,
        style = style,
        palette = palette,
    )
end


function Gadfly.plot(ai::AlignIntermediary, df, df2, idx=1, xs=:x, ys=:y)
    df2p = align(df2, ai, idx, xs, ys)
    xex, yex, x2ex, y2ex = extrema(df[:,xs]), extrema(df[:,ys]), extrema(df2p[:,xs]), extrema(df2p[:,ys])
    xmin, xmax = max(xex[1], x2ex[1]), min(xex[2], x2ex[2])
    ymin, ymax = max(yex[1], y2ex[1]), min(yex[2], y2ex[2])
    extra = 0.1*max(xmax-xmin, ymax-ymin)
    plot(
        layer(x=xs, y=ys, color=[ colorant"red" ], df2p, Geom.point, Theme(alphas=[0.4], highlight_width=0pt, point_size=2pt)),
        layer(x=xs, y=ys, color=[ colorant"blue" ], df, Geom.point, Theme(alphas=[1.0], highlight_width=0pt, point_size=2pt)),
        Guide.xlabel("X"), Guide.ylabel("Y"), 
        Coord.cartesian(xmin=xmin-extra, xmax=xmax+extra, ymin=ymin-extra, ymax=ymax+extra))
end