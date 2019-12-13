
using FileIO

struct Zepellin
    headerfile::String
    header::Dict{String,String}
    columns::Vector{Vector{String}}
    elms::Array{Element}
    data::DataFrame

    function Zepellin(headerfilename::String)
        new(headerfilename, loadZep(headerfilename)...)
    end
    function Zepellin( #
        headerfile::String,
        header::Dict{String,String},
        columns::Vector{Vector{String}},
        data::DataFrame,
    )
        new(
            headerfile,
            header,
            columns,
            filter(col -> col in names(data), map(z -> Symbol(elements[z].symbol), 1:94)),
            data,
        )
    end
end

Base.copy(z::Zepellin) = Zepellin(z.headerfile, copy(z.header), copy(z.columns), copy(z.data))

function loadZep(headerfilename::String)
    remapcolumnnames = Dict{String,String}(
        "PART#" => "NUMBER",
        "PARTNUM" => "NUMBER",
        "FIELD#" => "FIELD",
        "FIELDNUM" => "FIELD",
        "MAGFIELD#" => "MAGFIELD",
        "MAGFIELDNUM" => "MAGFIELD",
        "X_ABS" => "XABS",
        "Y_ABS" => "YABS",
        "X_DAC" => "XDAC",
        "Y_DAC" => "YDAC",
        "XCENT" => "XDAC",
        "YCENT" => "YDAC",
        "X_FERET" => "XFERET",
        "Y_FERET" => "YFERET",
        "DVG" => "DAVG",
        "DAVE" => "DAVG",
        "PERIM" => "PERIMETER",
        "ORIENT" => "ORIENTATION",
        "MAG_INDEX" => "MAGINDEX",
        "FIRST_ELEM" => "FIRSTELM",
        "SECOND_ELEM" => "SECONDELM",
        "THIRD_ELEM" => "THIRDELM",
        "FOURTH_ELEM" => "FOURTHELM",
        "ATOMICNUMBER1" => "FIRSTELM",
        "ATOMICNUMBER2" => "SECONDELM",
        "ATOMICNUMBER3" => "THIRDELM",
        "ATOMICNUMBER4" => "FOURTHELM",
        "FIRST_CONC" => "COUNTS1",
        "SECOND_CONC" => "COUNTS2",
        "THIRD_CONC" => "COUNTS3",
        "FOURTH_CONC" => "COUNTS4",
        "FIRST_PCT" => "FIRSTPCT",
        "SECOND_PCT" => "SECONDPCT",
        "THIRD_PCT" => "THIRDPCT",
        "FOURTH_PCT" => "FOURTHPCT",
        "PCT1" => "FIRSTPCT",
        "PCT2" => "SECONDPCT",
        "PCT3" => "THIRDPCT",
        "PCT4" => "FOURTHPCT",
        "TYPE(4ET#)" => "TYPE4ET",
        "VOID_AREA" => "VOIDAREA",
        "RMS_VIDEO" => "RMSVIDEO",
        "FIT_QUAL" => "FIT_QUAL",
        "VERIFIED_CLASS" => "VERIFIEDCLASS",
        "EDGE_ROUGHNESS" => "EDGEROUGHNESS",
        "COMP_HASH" => "COMPHASH",
        "PSEM_CLASS" => "CLASS",
    )
    columnnames(cols) = map(cn -> get(remapcolumnnames, cn, cn), map(c -> c[1], cols))
    header, columns, hdr = Dict{String,String}(), [], true
    for line in readlines(headerfilename)
        if hdr
            p = findfirst(c -> c == '=', line)
            if !isnothing(p)
                (k, v) = line[1:p-1], line[p+1:end]
                header[k] = v
            end
            hdr = !isequal(uppercase(k), "PARTICLE_PARAMETERS")
        else
            push!(columns, string.(strip.(split(line, "\t"))))
        end
    end
    pxz = CSV.File(
        replace(headerfilename, r".[h|H][d|D][z|Z]$" => ".pxz"),
        header = columnnames(columns),
        normalizenames = true,
    ) |> DataFrame

    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    sortedkeys = sort(collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(header))), lt = sortclasses)
    clsnames = map(c -> header[c], sortedkeys)
    pxz = hcat(pxz, DataFrame(CLASSNAME = map(cl -> get(clsnames, convert(Int, cl) + 1, "####"), pxz[:, :CLASS])))
    categorical!(pxz, :CLASS, compress = true) # Convert class column to pxz
    categorical!(pxz, :CLASSNAME, compress = true) # Convert class column to pxz
    cols = names(pxz)
    elms = filter(z -> Symbol(z.symbol) in cols, PeriodicTable.elements[1:94])
    return (header, columns, elms, pxz)
end

function Base.show(io::IO, zep::Zepellin)
    print(io, "Zepellin[$(zep.headerfile),$(size(zep.data))]")
end

function classes(zep::Zepellin)
    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    return (c -> get(zep.header, c, "")).(sort(
        collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(zep.header))),
        lt = sortclasses,
    ))
end

elements(zep::Zepellin) = zep.elms

function header(zep::Zepellin)
    keep(k) = !mapreduce(ty -> startswith(k, uppercase(ty)), (a, b) -> a || b, ("CLASS", "ELEM", "MAG"))
    kys = filter(keep, keys(zep.header))
    return DataStructures.SortedDict(filter(kv -> kv[1] in kys, zep.header))
end

data(zep::Zepellin) = zep.data


"""
    zep[123] # where zep is a Zepellin

Returns the Spectrum (with images) associated with the particle at row
"""
Base.getindex(zep::Zepellin, row::Int) = spectrum(zep, row, true)

"""
    spectrum(zep::Zepellin, row::Int, withImgs = true)::Union{Spectrum, missing}

Returns the Spectrum (with images) associated with the particle at row.  If withImgs
is true, the associated image or images are read.
"""
function spectrum(zep::Zepellin, row::Int, withImgs = true)::Union{Spectrum,Missing}
    mag = hasproperty(zep, :MAG) ? zep[row, :MAG] : 0
    file = joinpath(dirname(zep.headerfile), "MAG$(mag)", "00000$(row)"[end-4:end] * ".tif")
    if !isfile(file)
        file = joinpath(dirname(zep.headerfile), "MAG$(mag)", "00000$(row)"[end-3:end] * ".tif")
    end
    at = missing
    if isfile(file)
        try
            at = readAspexTIFF(file, withImgs = withImgs)
            at[:Name] = "P[$(zep[row, :Number]), $(zep[row, :ClassName])]"
            at[:Signature] = Dict(sig[z] => zep.data[row, Symbol(elm.symbol)] for elm in zep.elms)
        catch
            @info "$(file) does not appear to be a valid ASPEX spectrum TIFF."
        end
    end
    return at
end
