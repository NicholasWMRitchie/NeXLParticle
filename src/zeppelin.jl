
using FileIO

Base.convert(::Type{Symbol}, elm::Element) = Symbol(uppercase(elm.symbol))

struct Zeppelin
    headerfile::String
    header::Dict{String,String}
    elms::Array{Element}
    data::DataFrame

    function Zeppelin(headerfilename::String)
        new(headerfilename, loadZep(headerfilename)...)
    end

    function Zeppelin( #
        headerfile::String,
        header::Dict{String,String},
        data::DataFrame,
    )
        new(
            headerfile,
            header,
            filter(elm -> convert(Symbol,elm) in names(data), PeriodicTable.elements[1:94]),
            data,
        )
    end
end

Base.copy(z::Zeppelin) = Zeppelin(z.headerfile, copy(z.header), copy(z.data))

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
    columnnames(cols) = uppercase.(map(cn -> get(remapcolumnnames, cn, cn), map(c -> c[1], cols)))
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
    return (header, elms, pxz)
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

data(zep::Zeppelin; sortCol=:None, rev=false) =
    sortCol≠:None ? sort(zep.data, sortCol, rev=rev) : zep.data

function data(zep::Zeppelin, rows; sortCol=:None, rev=false)
    res = zep.data[intersect(rows, eachparticle(zep)), :]
    return sortCol≠:None ? sort(res, sortCol, rev=rev) : res
end

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
            at[:BeamEnergy] = beamenergy(zep, get(at, :BeamEnergy, 20.0e3))
            at[:ProbeCurrent] = get(at, :ProbeCurrent, probecurrent(zep, 1.0))
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

function beamenergy(zep::Zeppelin, def=missing)
    number(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)",v)[1])
    val = def
    try
        v = get(zep.header, "ACCELERATING_VOLTAGE", missing)
        val = !ismissing(v) ? number(v)*1.0e3 : def
    catch
        # Ignore it...
    end
    return val
end

function probecurrent(zep::Zeppelin, def=missing)
    number(v) = parse(Float64, match(r"([+-]?[0-9]+[.]?[0-9]*)",v)[1])
    val = def
    try
        v = get(zep.header, "PROBE_CURRENT", missing)
        val = !ismissing(v) ? number(v) : def
    catch
        # Ignore it...
    end
    return val
end


"""
    maxparticle(zep::Zeppelin, rows=1:100000000)

Computes the maxparticle spectrum from the specified rows in a Zeppelin dataset.
"""
function maxparticle(zep::Zeppelin, rows=1:100000000)
    rows = intersect(rows,eachparticle(zep))
    maxi(a, b) = collect(max(a[i],b[i]) for i in eachindex(a))
    mp = reduce(maxi, filter(s->!ismissing(s),(r->spectrum(zep, r, false)).(rows)))
    # Now copy it out as a spectrum
    i = findfirst(s->!ismissing(s),(r->spectrum(zep, r, false)).(rows))
    sp = spectrum(zep, i, false)
    props = copy(sp.properties)
    props[:Name] = "MaxParticle[$(basename(zep.headerfile))]"
    return Spectrum(sp.energy, mp, props)
end

"""
    rowsMax(zep::Zeppelin, col::Symbol, n=20)

Returns the row indices associated with the `n` maximum values in the `col` column.

Examples:

    # Plot ten largest by :DAVG
    plot(zep, rowsMax(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    data(zep, rowsMax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsMax(zep::Zeppelin, col::Symbol, n::Int=20)
    d = collect(zip(1:250,zep.data[:,col]))
    s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
    return getindex.(s[intersect(eachparticle(zep),1:n)],1)
end

"""
    rowsMin(zep::Zeppelin, col::Symbol, n=20)

Returns the row indices associated with the `n` minimum values in the `col` column.

Examples:

    # Plot ten smallest by :DAVG
    plot(zep, rowsMin(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    data(zep, rowsMax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsMin(zep::Zeppelin, col::Symbol, n::Int=20)
    d = collect(zip(eachparticle(zep),zep.data[:,col]))
    s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
    return getindex.(s[intersect(eachparticle(zep),1:n)],1)
end


"""
    rowsMax(zep::Zeppelin, col::Symbol, classname::AbstractString, n=20)

Returns the row indices associated with the `n` maximum values in the `col` column that
are also members of the specified `classname`.

Example:

    # Plot ten largest "Calcite" particles by :DAVG
    plot(zep, rowsMax(zep, :DAVG, "Calcite", 10))
"""
function rowsMax(zep::Zeppelin, col::Symbol, classname::AbstractString, n::Int=20)
    cr = filter(r->zep.data[r,:CLASSNAME]==classname, eachparticle(zep))
    d = collect(zip(cr, zep.data[cr,col]))
    s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
    return getindex.(s[intersect(eachindex(cr),1:n)],1)
end

"""
    rowsMin(zep::Zeppelin, col::Symbol, classname::AbstractString, n=20)

Returns the row indices associated with the `n` minimum  values in the `col` column that
are also members of the specified `classname`.

Example:

    # Plot ten smallest "Calcite" particles by :DAVG
    plot(zep, rowsMin(zep, :DAVG, "Calcite", 10))
"""
function rowsMin(zep::Zeppelin, col::Symbol, classname::AbstractString, n::Int=20)
    # Extract the rows for the class
    cr = filter(r->zep.data[r,:CLASSNAME]==classname, eachparticle(zep))
    # Combine (row, datum) for each row
    d = collect(zip(cr, zep.data[cr,col]))
    # Sort by datum
    s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
    # Return up to the first n rows
    return getindex.(s[intersect(eachindex(cr),1:n)],1)
end

"""
    rowsClass(zep::Zeppelin, classname::AbstractString, shuffle=false)

Returns the row indices associated with the specified `classname`.  If `shuffle` is true,
the row indices are shuffled so that a randomized n can be plucked off to plot or other.

Example:

    # Plot a randomized selection of 10 "Calcite" particles
    plot(zep, rowsClass(zep,"Calcite",true)[1:10], xmax=8.0e3)
"""
function rowsClass(zep::Zeppelin, classname::AbstractString, shuffle=false)
    res = filter(i->zep.data[i,:CLASSNAME]==classname, eachparticle(zep))
    return shuffle ? Random.shuffle(res) : res
end

const MORPH_COLS = [ :NUMBER, :XABS, :YABS, :DAVG, :DMIN, :DMAX, :DPERP, :PERIMETER, :AREA ]
const CLASS_COLS = [ :CLASS, :CLASSNAME, :VERIFIEDCLASS, :IMPORTANCE ]
const COMP_COLS = [ :FIRSTELM, :FIRSTPCT, :SECONDELM, :SECONDPCT, :THIRDELM,
    :THIRDPCT, :FOURTHELM, :FOURTHPCT, :COUNTS, :CLASSNAME ]
const ALL_ELMS = convert.(Symbol, elements[1:95])

allelms(zep::Zeppelin) = intersect(ALL_ELMS, names(zep.data))

const ALL_COMPOSITIONAL_COLUMNS = ( :FIRST, :FIRSTELM, :SECONDELM, :THIRDELM, :FOURTHELM, :FIRSTELM, :SECONDELM, :THIRDELM,
   :FOURTHELM, :COUNTS1, :COUNTS2, :COUNTS3, :COUNTS4, :FIRSTPCT, :SECONDPCT, :THIRDPCT, :FOURTHPCT, :TYPE4ET, :COUNTS,
   :FITQUAL, :COMPHASH )
const ALL_CLASS_COLS = ( :CLASS, :CLASSNAME, :VERIFIEDCLASS, :IMPORTANCE )

RJLG_ZEPPELIN=format"RJLG Zeppelin"

load(file::File{RJLG_ZEPPELIN}) = Zeppelin(file.filename)

FileIO.add_format(RJLG_ZEPPELIN, iszeppelin, [ ".hdz" ])
