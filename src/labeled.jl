using Compose
using FileIO
using ImageIO

struct LabeledView
    label::String
    image::Array
end

function labeledimage(label::String, image::Array)
    io = IOBuffer(maxsize=10*1024*1024)
    save(Stream(format"PNG",io), image)
    pix = max(size(image,1),size(image,2))
    scaleX, scaleY = size(image,1)/pix, size(image,2)/pix
    return compose(context(),
            (context(), bitmap("image/png", take!(io), 0.5*(1.0 - 0.8*scaleX), 0.5*(1.0-0.8*scaleY), 0.8*scaleX, 0.8*scaleY)),
            (context(), text(0.5, 0.96, label, hcenter, vbottom), stroke("transparent"), fill("black"), fontsize(0.08h))
    )
end

function labeledimages(labels::AbstractVector{<:AbstractString}, images::AbstractVector{<:AbstractArray}; ncols=3, halign = hleft)
    @assert length(labels)==length(images)
    nrows = (length(images)+ncols-1) ÷ ncols
    sc = 1.0/max(ncols, nrows)
    tmp = []
    for i in eachindex(labels)
        r, c = (i-1) ÷ ncols, (i-1) % ncols
        push!(tmp, (context(c*sc, r*sc, sc, sc), labeledimage(labels[i], images[i])))
    end
    return compose(context(), tmp...)
end
