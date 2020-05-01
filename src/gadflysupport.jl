using .Gadfly
using Random

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
