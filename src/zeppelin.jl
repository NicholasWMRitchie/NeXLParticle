
using FileIO

Base.convert(::Type{Symbol}, elm::Element) = Symbol(uppercase(elm.symbol))

struct Zeppelin
    headerfile::String
    header::Dict{String,String}
    columns::Vector{Vector{String}}
    elms::Array{Element}
    data::DataFrame

    function Zeppelin(headerfilename::String)
        new(headerfilename, loadZep(headerfilename)...)
    end

    function Zeppelin( #
        headerfile::String,
        header::Dict{String,String},
        columns::Vector{Vector{String}},
        data::DataFrame,
    )
        new(
            headerfile,
            header,
            columns,
            filter(col -> col in names(data), map(elm -> convert(Symbol,elm)), PeriodicTable.elements[1:94]),
            data,
        )
    end
end

Base.copy(z::Zeppelin) = Zeppelin(z.headerfile, copy(z.header), copy(z.columns), copy(z.data))

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
        "LIVE_TIME" => "LIVETIME",
        "FIT_QUAL" => "FITQUAL",
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
    elms = filter(z -> convert(Symbol, z) in cols, PeriodicTable.elements[1:94])
    return (header, columns, elms, pxz)
end

function Base.show(io::IO, zep::Zeppelin)
    print(io, "Zeppelin[$(zep.headerfile),$(size(zep.data))]")
end

function classes(zep::Zeppelin)
    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    return (c -> get(zep.header, c, "")).(sort(
        collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(zep.header))),
        lt = sortclasses,
    ))
end

NeXLCore.elms(zep::Zeppelin) = zep.elms

function header(zep::Zeppelin)
    keep(k) = !mapreduce(ty -> startswith(k, uppercase(ty)), (a, b) -> a || b, ("CLASS", "ELEM", "MAG"))
    kys = filter(keep, keys(zep.header))
    return DataStructures.SortedDict(filter(kv -> kv[1] in kys, zep.header))
end

data(zep::Zeppelin) = zep.data

eachparticle(zep::Zeppelin) = 1:size(zep.data, 1)

"""
    zep[123] # where zep is a Zeppelin

Returns the Spectrum (with images) associated with the particle at row
"""
Base.getindex(zep::Zeppelin, row::Int) = spectrum(zep, row, true)

"""
    spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum, missing}

Returns the Spectrum (with images) associated with the particle at row.  If withImgs
is true, the associated image or images are read.
"""
function spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum,Missing}
    mag = hasproperty(zep.data, :MAG) ? convert(Int,trunc(zep.data[row, :MAG])) : 0
    part = zep.data[row, :NUMBER]
    file = joinpath(dirname(zep.headerfile), "MAG$(mag)", "00000$(part)"[end-4:end] * ".tif")
    if !isfile(file)
        file = joinpath(dirname(zep.headerfile), "MAG$(mag)", "00000$(part)"[end-3:end] * ".tif")
    end
    at = missing
    if isfile(file)
        try
            at = readAspexTIFF(file, withImgs = withImgs)
        catch
            @info "$(file) does not appear to be a valid ASPEX spectrum TIFF."
        end
        try
            at[:Name] = "P[$part, $(zep.data[row, :CLASSNAME])]"
            at[:Signature] = filter(kv->kv[2]>0.0, Dict(elm => zep.data[row, convert(Symbol,elm)] for elm in zep.elms))
        catch
            @info "Error adding properties to ASPEX TIFF file."
        end
    else
        @info "Can't find $file."
    end
    return at
end

function iszeppelin(filename::String)
    res=false
    open(filename) do ios
        seekstart(ios)
        if isequal(uppercase(String(read(ios,11))),"PARAMETERS=")
            readline(ios) # read the rest of the line
            res = isequal(uppercase(readline(ios)),"HEADER_FMT=ZEPP_1")
        end
    end
    return res
end

RJLG_ZEPPELIN=format"RJLG Zeppelin"

load(file::File{RJLG_ZEPPELIN}) = Zeppelin(file.filename)

FileIO.add_format(RJLG_ZEPPELIN, iszeppelin, [ ".hdz" ])

const COMPOSITIONAL_COLUMNS = ( :FIRST, :FIRSTELM, :SECONDELM, :THIRDELM, :FOURTHELM, :FIRSTELM, :SECONDELM, :THIRDELM,
   :FOURTHELM, :COUNTS1, :COUNTS2, :COUNTS3, :COUNTS4, :FIRSTPCT, :SECONDPCT, :THIRDPCT, :FOURTHPCT, :FIRSTPCT,
:SECONDPCT, :THIRDPCT, :FOURTHPCT, :TYPE4ET, :FITQUAL, :COMPHASH )

const CLASS_COLUMNS =  ( :VERIFIEDCLASS, :CLASS, :CLASSNAME )
