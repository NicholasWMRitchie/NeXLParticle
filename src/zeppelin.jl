using CSV
using DataStructures
using Random
using StatsBase
using DataAPI
using StringEncodings

Base.convert(::Type{Symbol}, elm::Element) = Symbol(uppercase(elm.symbol))

struct Zeppelin
    headerfile::String
    header::Dict{String,String}
    elms::Vector{Element}
    classes::Dict{String, String}
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
    function _massagehdz(header)
        res = copy(header)
        foreach(hi -> if haskey(header, hi) delete!(res, hi) end, ( "PARTICLE_PARAMETERS", "PARAMETERS", "HEADER_FMT" ))
        foreach(f->if startswith(f,"ELEM") delete!(res,f) end, keys(header))
        foreach(f->if startswith(f,"CLASS") delete!(res,f) end, keys(header))
        return res
    end
    # Mostly converts columns to categorical as desired...
    function _massagepxz(header, pxz)::DataFrame
        res = copy(pxz)
        sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
        sortedkeys = sort(collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(header))), lt = sortclasses)
        clsnames = append!(["--None--"], map(c -> header[c], sortedkeys))
        for col in (:CLASS, :VERIFIED_CLASS)
            ic = findfirst(nm->nm==col, names(res))
            if !isnothing(ic)
                cls = categorical(Union{String,Missing}[ clsnames... ], ordered=true)[1:0] # empty but with correct levels
                foreach(name->push!(cls, name), map(cl -> get(clsnames, convert(Int, cl) + 2, "####"), pxz[:, col]))
                select!(res, Not(ic))
                insertcols!(res, ic, col=>cls)
            end
        end
        for col in ( :FIRSTELM, :SECONDELM, :THIRDELM, :FOURTHELM)
            ic = findfirst(nm->nm==col, names(res))
            if !isnothing(ic)
                fee=categorical( Union{Element,Missing}[ elements[1:95]...],false,ordered=true)[1:0]
                foreach(z->push!(fee, z>0 ? elements[z] : missing), pxz[:,col])
                select!(res, Not(ic))
                insertcols!(res, ic, col=>fee)
            end
        end
        for col in ( :XDAC, :YDAC, :TYPE_4ET_, :TYPE4ET )
            ic = findfirst(nm->nm==col, names(res))
            if !isnothing(ic)
                select!(res, Not(ic))
            end
        end
        return res
    end
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
        "TYPE(4ET)#" => "TYPE4ET",
        "TYPE(4ET)" => "TYPE4ET",
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
    open(hdzfilename, enc"WINDOWS-1252","r") do f
        for line in readlines(f)
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
    end
    pxz = CSV.File(
        replace( hdzfilename, r".[h|H][d|D][z|Z]$" => ".pxz"),
        header = columnnames(columns),
        delim='\t',
        normalizenames = true,
        missingstring="-",
        decimal='.'
    ) |> DataFrame
    elms = filter(z -> convert(Symbol, z) in names(pxz), PeriodicTable.elements[1:94])
    classes = Dict{String, String}(key=>header[key] for key in filter(f->startswith(f, "CLASS") && (f ≠ "CLASSES"),keys(header)))
    return (_massagehdz(header), elms, classes, _massagepxz(header, pxz))
end

function Base.show(io::IO, zep::Zeppelin)
    print(io, "Zeppelin[$(zep.headerfile),$(size(zep.data))]")
end

function classes(zep::Zeppelin)
    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    return (c -> zep.classes[c]).(sort(collect(keys(zep.classes)),lt = sortclasses,))
end

NeXLCore.elms(zep::Zeppelin) = zep.elms

header(zep::Zeppelin) = SortedDict(zep.header)

"""
    asa(::Type{DataFrame}, zep::Zeppelin; rows=missing, columns=missing, sortcol=:None, rev=false)

Create a DataFrame containing the particle data.  Optionally specify which rows and columns to include and
whether to sort the resulting table.
"""
function NeXLUncertainties.asa(::Type{DataFrame}, zep::Zeppelin; rows=missing, sortcol=:None, columns=missing, rev=false)
    rs = ismissing(rows) ? eachparticle(zep) : intersect(eachparticle(zep), rows)
    cols = ismissing(columns) ? names(zep.data) : intersect(columns,names(zep.data))
    res = zep.data[rs, cols]
    return sortcol≠:None ? sort(res, sortcol, rev=rev) : res
end

eachparticle(zep::Zeppelin) = 1:size(zep.data, 1)

function DataAPI.describe(zeps::AbstractVector{Zeppelin}; dcol=:DAVG, nelms=2)
    function build(zep)
        df=describe(zep, dcol=dcol, nelms=nelms)
        insertcols!(df, 1, :Dataset=>[zep.header["DESCRIPTION"] for _ in 1:size(df,1)])
        return df
    end
    return vcat(map(z->build(z),zeps))
end

function DataAPI.describe(zep::Zeppelin; dcol=:DAVG, nelms=3)
    function sortedstats(rows)
        tmp = [ ( elm, summarystats(zep[rows,convert(Symbol, elm)])) for elm in zep.elms ]
        return sort!(tmp, lt=(v1,v2)->isless(v1[2].mean,v2[2].mean),rev=true)
    end
    cnelms = min(length(zep.elms), nelms)
    ds = summarystats(zep[:,dcol])
    clss = String["All"]
    szs = Int[size(zep.data,1)]
    minds = Float64[ ds.min ]
    medds = Float64[ ds.median ]
    maxds = Float64[ ds.max ]
    sbcm = StatsBase.countmap(zep[:,:CLASS])
    for (cn, cx) in sbcm
        push!(clss,string(cn))
        push!(szs, cx)
        rows = rowsclass(zep,cn)
        @assert length(rows)==cx
        ds = summarystats(zep[rows,dcol])
        push!(minds, ds.min)
        push!(medds, ds.median)
        push!(maxds, ds.max)
    end
    dfres = DataFrame(Symbol("Class")=>clss, Symbol("Count")=>szs,Symbol("Min[$dcol]")=>minds, Symbol("Median[$dcol]")=>medds, Symbol("Max[$dcol]")=>maxds)
    # Add All elemental row
    ss = sortedstats(eachparticle(zep))
    elmdfs=DataFrame[ ]
    uem, uef = Union{Element,Missing}, Union{Float64,Missing}
    for ne in 1:cnelms
        cols = Symbol.( [ "Elm[$ne]", "Min[$ne]", "Median[$ne]", "Max[$ne]" ])
        sst = ss[ne][2]
        push!(elmdfs, DataFrame(cols[1]=>uem[ ss[ne][1] ], cols[2]=>uef[sst.min], cols[3]=>uef[sst.median], cols[4]=>uef[sst.max]))
    end
    for (cn, cx) in sbcm
        try
            ss = sortedstats(rowsclass(zep, cn))
            for ne in 1:cnelms
                sst = ss[ne][2]
                push!(elmdfs[ne], [ ss[ne][1], sst.min, sst.median, sst.max ])
            end
        catch
            for ne in 1:cnelms
                push!(elmdfs[ne], [ missing, missing, missing, missing ])
            end
        end
    end
    for elmdf in elmdfs
        for col in names(elmdf)
            insertcols!(dfres, size(dfres,2)+1, col=>elmdf[:,col])
        end
    end
    return dfres # sort!(res,(:Count, :Class), rev=true)
end


"""
    zep[123] # where zep is a Zeppelin

Returns the Spectrum (with images) associated with the particle at row
"""
Base.getindex(zep::Zeppelin, row::Int) = spectrum(zep, row, true)

Base.getindex(zep::Zeppelin, rows, cols) = Base.getindex(zep.data, rows, cols)

Base.lastindex(zep::Zeppelin, axis::Integer) = Base.lastindex(zep.data, axis)

# Replace the default in PeriodicTable because it is too verbose...
# Base.show(io::IO, elm::Element) = print(io, elm.symbol)

"""
    spectrumfilename(zep::Zeppelin, row::Int, dir::AbstractString="MAG", ext::AbstractString=".tif")

Returns the name of the spectrum/image file for the particle in the specified row.
"""
function spectrumfilename(zep::Zeppelin, row::Int, dir::AbstractString="MAG", ext::AbstractString=".tif")
    mag = hasproperty(zep.data, :MAG) ? convert(Int,trunc(zep.data[row, :MAG])) : 0
    tmp = "$(zep.data[row, :NUMBER])"
    # First check if a spectrum file exists
    for n in 6:-1:4
        fn = joinpath(dirname(zep.headerfile), "$(dir)$(mag)", repeat('0',max(0,n-length(tmp)))*tmp*ext)
        if isfile(fn)
            return fn
        end
    end
    # Default case
    return joinpath(dirname(zep.headerfile), "$(dir)$(mag)", repeat('0',max(0,5-length(tmp)))*tmp*ext)
end

"""
    ParticleClassifier

A type that implements:

    classify(zep::Zeppelin, sr::ParticleClassifier)::CategoricalArray{String}
"""
abstract type ParticleClassifier end

"""
    spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum, missing}

Returns the Spectrum (with images) associated with the particle at row.  If withImgs
is true, the associated image or images are read.
"""
function spectrum(zep::Zeppelin, row::Int, withImgs = true)::Union{Spectrum,Missing}
    file, at = spectrumfilename(zep, row, "MAG", ".tif"), missing
    if isfile(file)
        try
            at = readAspexTIFF(file, withImgs = withImgs)
        catch err
            showerror(stderr, err)
            @info "$(file) does not appear to be a valid ASPEX spectrum TIFF."
        end
        try
            at[:BeamEnergy] = beamenergy(zep, get(at, :BeamEnergy, 20.0e3))
            at[:ProbeCurrent] = get(at, :ProbeCurrent, probecurrent(zep, 1.0))
            at[:Signature] = filter(kv->kv[2]>0.0, Dict(elm => zep.data[row, convert(Symbol,elm)] for elm in zep.elms))
            at[:Name] = "P[$(zep.data[row, :NUMBER]), $(zep.data[row, :CLASS])]"
        catch
            @info "Error adding properties to ASPEX TIFF file."
        end
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
    remapcolumnnames = Dict{Symbol,String}(
        :NUMBER => "PART#\t1\tINT16",
        :FIELD  => "FIELD#\t1\tINT16",
        :MAGFIELD  => "MAGFIELD#\t1\tINT16",
        :XABS  => "X_ABS\tmm\tFLOAT",
        :YABS  => "Y_ABS\tmm\tFLOAT",
        :XDAC  => "X_DAC\t1\tINT16",
        :YDAC  => "Y_DAC\t1\tINT16",
        :XFERET  => "X_FERET\tµm\tFLOAT",
        :YFERET  => "Y_FERET\tµm\tFLOAT",
        :DAVG  => "DAVE\tµm\tFLOAT",
        :DMAX  => "DMAX\tµm\tFLOAT",
        :DMIN  => "DMIN\tµm\tFLOAT",
        :DPERP  => "DPERP\tµm\tFLOAT",
        :ASPECT  =>"ASPECT\t1\tFLOAT",
        :AREA  => "AREA\tµm²\tFLOAT",
        :PERIMETER  => "PERIMETER\tµm\tFLOAT",
        :ORIENTATION  => "ORIENTATION\tdeg\tFLOAT",
        :LIVETIME  => "LIVE_TIME\ts\tFLOAT",
        :FITQUAL  => "FIT_QUAL\t1\tFLOAT",
        :MAG  => "MAG\t1\tINT16",
        :VIDEO  => "VIDEO\t1\tINT16",
        :IMPORTANCE  => "IMPORTANCE\t1\tINT16",
        :COUNTS  => "COUNTS\t1\tFLOAT",
        :MAGINDEX  => "MAG_INDEX\t1\tINT16",
        :FIRSTELM  => "FIRST_ELEM\t1\tINT16",
        :SECONDELM  => "SECOND_ELEM\t1\tINT16",
        :THIRDELM  => "THIRD_ELEM\t1\tINT16",
        :FOURTHELM  => "FOURTH_ELEM\t1\tINT16",
        :COUNTS1  => "FIRST_CONC\tcounts\tFLOAT",
        :COUNTS2  => "SECOND_CONC\tcounts\tFLOAT",
        :COUNTS3  => "THIRD_CONC\tcounts\tFLOAT",
        :COUNTS4  => "FOURTH_CONC\tcounts\tFLOAT",
        :FIRSTPCT  => "FIRST_PCT\t%\tFLOAT",
        :SECONDPCT  => "SECOND_PCT\t%\tFLOAT",
        :THIRDPCT  => "THIRD_PCT\t%\tFLOAT",
        :FOURTHPCT  => "FOURTH_PCT\t%\tFLOAT",
        :TYPE4ET  => "TYPE(4ET#)\t1\tLONG",
        :VOIDAREA  => "VOID_AREA\tµm²\tFLOAT",
        :RMSVIDEO  => "RMS_VIDEO\t1\tINT16",
        :FITQUAL  => "FIT_QUAL\t1\tFLOAT",
        :VERIFIEDCLASS  => "VERIFIED_CLASS\t1\tINT16",
        :EDGEROUGHNESS  => "EDGE_ROUGHNESS\t1\tFLOAT",
        :COMPHASH  => "COMP_HASH\t1\tLONG",
        :CLASS  => "PSEM_CLASS\t1\tINT16",
        :TYPE_4ET_ => "Type[4ET]\t1\tLONG",
    )
    merge!(remapcolumnnames, Dict( convert(Symbol,elm)  => "$(uppercase(elm.symbol))\t%(k)\tFLOAT" for elm in zep.elms))
    merge!(remapcolumnnames, Dict( Symbol("U_$(uppercase(elm.symbol))_")  => "U[$(uppercase(elm.symbol))]\t%(k)\tFLOAT" for elm in zep.elms))
    headeritems = copy(zep.header)
    # add back the element tags
    for (i,elm) in enumerate(zep.elms)
        headeritems["ELEM$(i-1)"] = "$(elm.symbol) $(z(elm)) 1"
    end
    headeritems["ELEMENTS"] = "$(length(zep.elms))"
    if :CLASS in names(zep.data)
        clsdata=zep.data[:,:CLASS]
        for (i,cls) in enumerate(clsdata.pool.index[1:end])
            # println("CLASS$(i-1) => $cls")
            headeritems["CLASS$(i-1)"]=cls
            headeritems["CLASSES"]="$(i+1)"
        end
    end
    headeritems["TOTAL_PARTICLES"] = "$(size(zep.data,1))"
    # write out the header
    colnames = names(zep.data)
    open(hdzfilename,"w") do ios
        println(ios, "PARAMETERS=$(length(headeritems)+size(zep.data,2)+2)")
        println(ios, "HEADER_FMT=ZEPP_1")
        foreach(hk->println(ios,"$(hk)=$(headeritems[hk])"), sort(collect(keys(headeritems))))
        println(ios, "PARTICLE_PARAMETERS=$(size(zep.data,2))")
        foreach(col->println(ios,remapcolumnnames[col]),colnames)
    end
    pxzfilename = replace( hdzfilename, r".[h|H][d|D][z|Z]"=>".pxz")
    zd = DataFrame(zep.data)
    # Replace categorical data with the index (either into CLASS# in header or z(elm))
    for (ic, col) in enumerate(names(zd))
        coldata=zd[:,col]
        if coldata isa CategoricalArray
            select!(zd,Not(ic))
            et = eltype(coldata)
            if et == CategoricalArrays.CategoricalString{UInt32}
                insertcols!(zd, ic, col=>[ coldata.refs[i]-1  for i in eachindex(coldata) ]) # CLASS
            else
                insertcols!(zd, ic, col=>[ coldata.refs[i]  for i in eachindex(coldata) ]) # Element
            end
        end
    end
    CSV.write(pxzfilename, zd, delim="\t", missingstring="-", writeheader=false)
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

function magdata(zep::Zeppelin, index::Int)
    h = split(zep.header["MAG_FMT"],isspace)
    v = split(zep.header["MAG$index"],isspace)
    return Dict( string(h[i])=>parse(Float64, v[i]) for i in 1:min(length(h),length(v)) )
end

"""
    randomsubset(zep::Zeppelin, maxRows::Int)

Creates an ordered random subselection of `maxRows` without replacement of the rows in the `zep` dataset.  If maxRows
is larger than the number of rows in `zep` then a UnitRange with all rows is returned.
"""
randomsubset(zep::Zeppelin, maxRows::Int) =
    return maxRows<size(zep.data,1) ? #
        sort(Random.shuffle(collect(eachparticle(zep)))[1:maxRows]) :
        eachparticle(zep)

"""
    residual(zep::Zeppelin, row::Int, withImgs=false)::Union{Spectrum,Missing}

Returns the residual spectrum for the particle at `row` or missing if it does not exist. Never returns images.
"""
function residual(zep::Zeppelin, row::Int, withImgs=false)::Union{Spectrum,Missing}
    res = missing
    try
        filename = spectrumfilename(zep, row, "Residual", ".msa")
        res = isfile(filename) ? readEMSA(filename) : missing
    catch
        # Ignore errors, just return missing
    end
    return res
end

"""
    maxparticle(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})

Computes the maxparticle spectrum from the specified particles by `rows`.
"""
function maxparticle(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    mp, firstspec = missing, missing
    for spec in filter(s->!ismissing(s), map(r->spectrum(zep,r,false),rows))
        @assert spec isa Spectrum
        firstspec = ismissing(firstspec) ? spec : firstspec
        cx = NeXLSpectrum.counts(spec)
        mp = ismissing(mp) ? cx : max.(mp, cx)
    end
    # Now copy it out as a spectrum
    props = copy(firstspec.properties)
    props[:Name] = "MaxParticle"
    return Spectrum(firstspec.energy, mp, props)
end


"""
    maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}}=1:1000000)

Computes the max-residual spectrum from the specified particles by `rows`.  The residual must have been
previously computed in quantify(zep, ....)
"""
function maxresidual(zep::Zeppelin, rows::Union{AbstractVector{Int},UnitRange{Int}})
    mp, firstspec = missing, missing
    for spec in filter(s->!ismissing(s), map(r->residual(zep,r,false),rows))
        @assert spec isa Spectrum
        firstspec = ismissing(firstspec) ? spec : firstspec
        cx = NeXLSpectrum.counts(spec)
        mp = ismissing(mp) ? cx : max.(mp, cx)
    end
    # Now copy it out as a spectrum
    props = copy(firstspec.properties)
    props[:Name] = "MaxResidual"
    return Spectrum(firstspec.energy, mp, props)
end

"""
    rowsmax(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)

Returns the row indices associated with the `n` maximum values in the `col` column.

Examples:

    # Plot ten largest by :DAVG
    plot(zep, rowsmax(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    asa(DataFrame, zep, rowsmax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsmax(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)
    if ismissing(classname)
        d = collect(zip(1:250,zep.data[:,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
        return getindex.(s[intersect(eachparticle(zep),1:n)],1)
    else
        cr = filter(r->zep.data[r,:CLASS]==classname, eachparticle(zep))
        d = collect(zip(cr, zep.data[cr,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=true)
        return getindex.(s[intersect(eachindex(cr),1:n)],1)
    end
end
"""
    rowsmin(zep::Zeppelin, col::Symbol; n::Int=10000000, classname::Union{Missing,AbstractString}=missing)

Returns the row indices associated with the `n` minimum values in the `col` column.

Examples:

    # Plot ten smallest by :DAVG
    plot(zep, rowsmin(zep, :DAVG, 10))
    # data from 100 most Iron-rich sorted by DAVG
    asa(DataFrame, zep, rowsmax(zep, :FE, 100), sortCol=:DAVG)
"""
function rowsmin(zep::Zeppelin, col::Symbol; n::Int=20, classname::Union{Missing,AbstractString}=missing)
    if ismissing(classname)
        d = collect(zip(eachparticle(zep),zep.data[:,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
        return getindex.(s[intersect(eachparticle(zep),1:n)],1)
    else
        cr = filter(r->zep.data[r,:CLASS]==classname, eachparticle(zep))
        d = collect(zip(cr, zep.data[cr,col]))
        s = sort(d, lt=(a1,a2)->isless(a1[2],a2[2]),rev=false)
        return getindex.(s[intersect(eachindex(cr),1:n)],1)
    end
end

"""
    rowsclass(zep::Zeppelin, classname::AbstractString; shuffle=false, n=1000000)

Returns the row indices associated with the specified `classname`.  If `shuffle` is true,
the row indices are shuffled so that a randomized n can be plucked off to plot or other.

Example:

    # Plot a randomized selection of 10 "Calcite" particles
    plot(zep, rowsclass(zep, "Calcite", shuffle=true, n=10), xmax=8.0e3)
"""
function rowsclass(zep::Zeppelin, classname::AbstractString; shuffle=false, n=1000000)
    res = filter(i->zep.data[i,:CLASS]==classname, eachparticle(zep))
    return (shuffle ? Random.shuffle(res) : res)[intersect(eachindex(res),1:n)]
end
rowsclass(zep::Zeppelin, classnames::AbstractVector{<:AbstractString}; shuffle=false, n=1000000) =
    mapreduce(cn->rowsclass(zep, cn, shuffle=shuffle, n=n), union, classnames)



"""
    Base.filter(filt::Function, zep::Zeppelin)

Use a function of the form `filt(row)::Bool` to filter a Zeppelin dataset returning a new Zeppelin dataset
with only the rows for which the function evaluated true.


Example:

    gsr = filter(row->startswith( String(row[:CLASS]), "GSR."), zep)

"""
Base.filter(filt::Function, zep::Zeppelin) = filter(r->filt(zep.data[r,:]), eachparticle(zep))

"""
    multiternary(
        zep::Zeppelin;
        omit = [ n"C", n"O" ],
        palette = TernPalette)

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
    for elm in filter(elm->!(elm in omit), zep.elms)
        sy = convert(Symbol, elm)
        push!(df, [sy, mean(zd[:,sy])])
    end
    sort!(df, :Mean, rev=true)
    elms = df[:,:Elm][1:min(size(df,1),6)]
    NeXLParticle.multiternary(zd, elms, :CLASS, title=zep.header["DESCRIPTION"], palette=palette, norm=100.0, fontsize=fontsize, font=font)
end


const MORPH_COLS = [ :NUMBER, :XABS, :YABS, :DAVG, :DMIN, :DMAX, :DPERP, :PERIMETER, :AREA ]
const CLASS_COLS = [ :CLASS, :VERIFIEDCLASS, :IMPORTANCE ]
const COMP_COLS = [ :FIRSTELM, :FIRSTPCT, :SECONDELM, :SECONDPCT, :THIRDELM,
    :THIRDPCT, :FOURTHELM, :FOURTHPCT, :COUNTS ]
const ALL_ELMS = convert.(Symbol, elements[1:95])

allelms(zep::Zeppelin) = intersect(ALL_ELMS, names(zep.data))

const ALL_COMPOSITIONAL_COLUMNS = ( :FIRST, :FIRSTELM, :SECONDELM, :THIRDELM, :FOURTHELM, :FIRSTELM, :SECONDELM, :THIRDELM,
   :FOURTHELM, :COUNTS1, :COUNTS2, :COUNTS3, :COUNTS4, :FIRSTPCT, :SECONDPCT, :THIRDPCT, :FOURTHPCT, :TYPE4ET, :COUNTS,
   :FITQUAL, :COMPHASH )
const ALL_CLASS_COLS = ( :CLASS, :VERIFIEDCLASS, :IMPORTANCE )

RJLG_ZEPPELIN=format"RJLG Zeppelin"

load(file::File{RJLG_ZEPPELIN}) = Zeppelin(file.filename)
save(f::File{RJLG_ZEPPELIN}, zep) = writeZep(zep, file.filename)

FileIO.add_format(RJLG_ZEPPELIN, iszeppelin, [ ".hdz" ])
