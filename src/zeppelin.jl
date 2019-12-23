using CSV
using DataStructures
using Random

Base.convert(::Type{Symbol}, elm::Element) = Symbol(uppercase(elm.symbol))

struct Zeppelin
    headerfile::String
    header::Dict{String,String}
    elms::Array{Element}
    data::DataFrame

    function Zeppelin( hdzfilename::String)
        new( hdzfilename, loadZep( hdzfilename)...)
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

function loadZep( hdzfilename::String)
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
        "FIT_QUAL" => "FITQUAL",
        "VERIFIED_CLASS" => "VERIFIEDCLASS",
        "EDGE_ROUGHNESS" => "EDGEROUGHNESS",
        "COMP_HASH" => "COMPHASH",
        "PSEM_CLASS" => "CLASS",
    )
    columnnames(cols) = uppercase.(map(cn -> get(remapcolumnnames, cn, cn), map(c -> c[1], cols)))
    header, columns, hdr = Dict{String,String}(), [], true
    for line in readlines( hdzfilename)
        if hdr
            p = findfirst(c -> c == '=', line)
            if !isnothing(p)
                (k, v) = line[1:p-1], line[p+1:end]
                header[k] = v
                hdr = !isequal(uppercase(k), "PARTICLE_PARAMETERS")
            end
        else
            push!(columns, string.(strip.(split(line, "\t"))))
        end
    end
    foreach(hi -> if haskey(header, hi) delete!(header, hi) end, ( "PARTICLE_PARAMETERS", "PARAMETERS", "HEADER_FMT" ))
    # Remove and regenerate "ELEMXX" data
    foreach(f->if startswith(f,"ELEM") delete!(header,f) end, keys(header))
    pxz = CSV.File(
        replace( hdzfilename, r".[h|H][d|D][z|Z]$" => ".pxz"),
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


function spectrumfilename(zep::Zeppelin, row::Int)
    mag = hasproperty(zep.data, :MAG) ? convert(Int,trunc(zep.data[row, :MAG])) : 0
    tmp = "$(zep.data[row, :NUMBER])"
    for n in 6:-1:4
        fn = joinpath(dirname(zep.headerfile), "MAG$(mag)", repeat('0',max(0,n-length(tmp)))*tmp*".tif")
        if isfile(fn)
            return fn
        end
    end
    return joinpath(dirname(zep.headerfile), "MAG$(mag)", repeat('0',max(0,5-length(tmp)))*tmp*".tif")
end



"""
    spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum, missing}

Returns the Spectrum (with images) associated with the particle at row.  If withImgs
is true, the associated image or images are read.
"""
function spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum,Missing}
    file, at = spectrumfilename(zep, row), missing
    if isfile(file)
        try
            at = readAspexTIFF(file, withImgs = withImgs)
        catch
            @info "$(file) does not appear to be a valid ASPEX spectrum TIFF."
        end
        try
            at[:BeamEnergy] = beamenergy(zep, get(at, :BeamEnergy, 20.0e3))
            at[:ProbeCurrent] = get(at, :ProbeCurrent, probecurrent(zep, 1.0))
            at[:Signature] = filter(kv->kv[2]>0.0, Dict(elm => zep.data[row, convert(Symbol,elm)] for elm in zep.elms))
            at[:Name] = "P[$(zep.data[row, :NUMBER]), $(zep.data[row, :CLASSNAME])]"
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

function writeZep(zep::Zeppelin,  hdzfilename::String)
        remapcolumnnames = Dict{String,String}(
            :NUMBER => "PART# 1 INT16",
            :FIELD  => "FIELD# 1 INT16",
            :MAGFIELD  => "MAGFIELD# 1 INT16",
            :XABS  => "X_ABS mm FLOAT",
            :YABS  => "Y_ABS mm FLOAT",
            :XDAC  => "X_DAC 1 INT16",
            :YDAC  => "Y_DAC 1 INT16",
            :XFERET  => "X_FERET µm FLOAT",
            :YFERET  => "Y_FERET µm FLOAT",
            :DAVG  => "DAVE µm FLOAT",
            :DMAX  => "DMAX µm FLOAT",
            :DMIN  => "DMIN µm FLOAT",
            :DPERP  => "DPERP µm FLOAT",
            :ASPECT  =>"ASPECT 1 FLOAT",
            :AREA  => "AREA µm² FLOAT",
            :PERIMETER  => "PERIMETER µm FLOAT",
            :ORIENTATION  => "ORIENTATION deg FLOAT",
            :LIVETIME  => "LIVE_TIME s FLOAT",
            :FITQUAL  => "FIT_QUAL 1 FLOAT",
            :MAG  => "MAG 1 INT16",
            :VIDEO  => "VIDEO 1 INT16",
            :IMPORTANCE  => "IMPORTANCE 1 INT16",
            :COUNTS  => "COUNTS 1 FLOAT",
            :MAGINDEX  => "MAG_INDEX 1 INT16",
            :FIRSTELM  => "FIRST_ELEM 1 INT16",
            :SECONDELM  => "SECOND_ELEM 1 INT16",
            :THIRDELM  => "THIRD_ELEM 1 INT16",
            :FOURTHELM  => "FOURTH_ELEM 1 INT16",
            :COUNTS1  => "FIRST_CONC counts FLOAT",
            :COUNTS2  => "SECOND_CONC counts FLOAT",
            :COUNTS3  => "THIRD_CONC counts FLOAT",
            :COUNTS4  => "FOURTH_CONC counts FLOAT",
            :FIRSTPCT  => "FIRST_PCT % FLOAT",
            :SECONDPCT  => "SECOND_PCT % FLOAT",
            :THIRDPCT  => "THIRD_PCT % FLOAT",
            :FOURTHPCT  => "FOURTH_PCT % FLOAT",
            :TYPE4ET  => "TYPE(4ET#) 1 LONG",
            :VOIDAREA  => "VOID_AREA µm² FLOAT",
            :RMSVIDEO  => "RMS_VIDEO 1 INT16",
            :FITQUAL  => "FIT_QUAL 1 FLOAT",
            :VERIFIEDCLASS  => "VERIFIED_CLASS 1 INT16",
            :EDGEROUGHNESS  => "EDGE_ROUGHNESS 1 FLOAT",
            :COMPHASH  => "COMP_HASH 1 LONG",
            :CLASS  => "PSEM_CLASS 1 INT16",
        )
        merge!(remapcolumnnames, Dict( convert(Symbol,elm)  => "$(elm.symbol) %(k) FLOAT" for elm in elms))
        merge!(remapcolumnnames, Dict( Symbol("U[$(uppercase(elm.symbol))]")  => "$(elm.symbol) %(k) FLOAT"  for elm in elms))
        headeritems = copy(zep.header)
        # add back the element tags
        foreach((i,elm) -> headeritems["ELEM%i"] = "$(elm.symbol) $(z(elm)) 1", enumerate(zep.elms))
        header["ELEMENTS"] = "$(length(zep.elms))"
        header["TOTAL_PARTICLES"] = "$(size(zep.data,1))"
        # write out the header
        colnames = filter(n-> n != :CLASSNAME, names(zep.data))
        IOStream( hdzfilename,"w") do ios
            println(ios, "PARAMETERS=$(length(header)+size(zep.data,2)+3)")
            println(ios, "HEADER_FMT=ZEPP_1")
            foreach(hk->println(ios,"$(hk)=$(headeritems[hk])"), sort(keys(headeritems)))
            println(ios, "PARTICLE_PARAMETERS=$(size(zep.data,2))")
            foreach(col->println(ios,remapcolumnnames[col]),colnames)
        end
        pxzfilename = replace( hdzfilename,r".[h|H][d|D][z|Z]",".pxz")
        CSV.write(pxzfilename, zep.data[:, colnames], delim="\t", missingstring="-", writeheader=false)
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
    maxparticle(zep::Zeppelin, maxrows=100000000)

Computes the maxparticle spectrum from up to `maxrows` spectra.  If the number of particles in `zep` is
larger than `maxrows` then a randomized subset is chosen.
"""
function maxparticle(zep::Zeppelin, maxrows=100000000)
    if maxrows<size(zep.data,1)
        rows = Random.shuffle(collect(eachparticle(zep)))[1:maxrows]
    else
        rows = eachparticle(zep)
    end
    maxi(a, b) = collect(max(a[i],b[i]) for i in eachindex(a))
    mp, firstspec = missing, missing
    for r in rows
        try
            spec = spectrum(zep, r, false)
            if !ismissing(spec)
                firstspec = ismissing(firstspec) ? spec : firstspec
                cx = counts(spec, Float64)
                mp = ismissing(mp) ? cx : maxi(mp, cx)
            end
        catch err
            @info "Failed $(err)"
            # Ignore it...
        end
    end
    # Now copy it out as a spectrum
    props = copy(firstspec.properties)
    props[:Name] = "MaxParticle[$(basename(zep.headerfile)[1:end-4])]"
    return Spectrum(firstspec.energy, mp, props)
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
save(f::File{RJLG_ZEPPELIN}, zep) = writeZep(zep, file.filename)

FileIO.add_format(RJLG_ZEPPELIN, iszeppelin, [ ".hdz" ])
