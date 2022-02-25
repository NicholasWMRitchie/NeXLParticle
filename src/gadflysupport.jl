using .Gadfly
using Colors
using StatsBase

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

function Gadfly.plot(zeps::AbstractArray{Zeppelin}; offset=(0.0,0.0), point_size=1.5pt, coords=missing)
    names = map(z->get(z.header,"DESCRIPTION",splitpath(z.headerfile)[end-1]), zeps)
    colors = map(i->NeXLPalette[1 + (i-1) % length(NeXLPalette)], eachindex(zeps))
    clsCx=collect(merge( [ StatsBase.countmap(repr.(z[:,:CLASS])) for z in zeps]...))
    sort!(clsCx, lt=(a,b)->a[2]<b[2], rev=true)
    shlen = length(Theme().point_shapes)
    if length(clsCx) >= shlen
        clsCx = clsCx[1:(shlen-1)]
        push!(clsCx,"Other"=>0)
    end
    shIdx=Dict(clsCx[i][1]=>i for i in eachindex(clsCx))
    xmin = minimum(z->minimum(z.data[:,:XABS]), zeps)
    xmax = maximum(z->maximum(z.data[:,:XABS]), zeps)
    ymin = minimum(z->minimum(z.data[:,:YABS]), zeps)
    ymax = maximum(z->maximum(z.data[:,:YABS]), zeps)
    dmax = ceil(max(xmax-xmin, ymin-ymax)/2)
    xc, yc = round((xmax+xmin)/2), round((ymin+ymax)/2)
    coords = Coord.cartesian(xmin=xc-dmax,xmax=xc+dmax,ymin=yc-dmax,ymax=yc+dmax)
    plot(
        Theme(point_size=1.5pt),
        map(enumerate(zeps)) do (i,z)
            layer(
                x=map(r->r.XABS+(i-1)*offset[1], eachrow(z.data)),
                y=map(r->r.YABS+(i-1)*offset[2], eachrow(z.data)), 
                shape=map(r->get(shIdx, repr(r.CLASS), length(clsCx)), eachrow(z.data)),
                Geom.point, 
                Theme(default_color=colors[i], point_size=point_size)
            )
        end...,
        Guide.xlabel("X (mm)"), Guide.ylabel("Y (mm)"), 
        Guide.manual_color_key("Particle Data Set", names, colors),
        Guide.shapekey(title="Particle Class", labels=map(cc->cc[1],clsCx)),
        coords
    )
end

function Gadfly.plot(::Type{Histogram}, zeps::AbstractArray, sym::Symbol; bincount=40, limits=missing)
    names = map(z->get(z.header,"DESCRIPTION",splitpath(z.headerfile)[end-1]), zeps)
    df = vcat(
        map(enumerate(zeps)) do (i, z) 
            DataFrame(Dataset=fill(names[i],nrow(z.data)), v=z.data[:,sym])
        end...
    )
    if ismissing(limits)
        limits = ( min=floor(minimum(z->minimum(z.data[:,sym]), zeps)), max=ceil(maximum(z->maximum(z.data[:,sym]), zeps)) )
    end
    h = Geom.histogram(bincount=bincount, position=:identity, limits=limits)
    plot(df, x=:v, color=:Dataset, alpha=[0.5], h,
        Guide.xlabel(repr(sym)[2:end]), Guide.ylabel("Count")        
    )
end