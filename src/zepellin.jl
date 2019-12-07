
struct Zepellin
    headerfile::String
    header::Dict{String,String}
    columns::Vector{Vector{String}}
    data::DataFrame

    function Zepellin(headerfilename::String)
        new(headerfilename, loadZep(headerfilename)...)
    end
    function Zepellin( #
        headerfile::String,
        header::Dict{String,String},
        columns::Vector{Vector{String}},
        data::DataFrame
    )
        new(headerfile, header, columns, data)
    end
end

Base.copy(z::Zepellin) =
    Zepellin(z.headerfile, copy(z.header), copy(z.columns), copy(z.data) )

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
            p = findfirst(c->c=='=', line)
            if !isnothing(p)
                (k, v) = line[1:p-1], line[p+1:end]
                header[k] = v
            end
            hdr = !isequal(uppercase(k), "PARTICLE_PARAMETERS")
        else
            push!(columns, string.(strip.(split(line, "\t"))))
        end
    end
    pxz = CSV.File(replace(headerfilename, ".hdz" => ".pxz"), header = columnnames(columns), normalizenames = true) |> DataFrame

    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    sortedkeys = sort(collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(header))), lt = sortclasses)
    clsnames = map(c -> header[c], sortedkeys)
    pxz=hcat(pxz, DataFrame(CLASSNAME=map(cl->get(clsnames,convert(Int,cl)+1,"####"),pxz[:,:CLASS])))
    categorical!(pxz, :CLASS, compress=true) # Convert class column to pxz
    categorical!(pxz, :CLASSNAME, compress=true) # Convert class column to pxz
    return (header, columns, pxz)
end

function Base.show(io::IO, zep::Zepellin)
    print(io,"Zepellin[$(zep.headerfile),$(size(zep.data))]")
end


function classes(zep::Zepellin)
    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int, c2[6:end]))
    return (c -> get(zep.header, c, "")).(sort(collect(filter(c -> !isnothing(match(r"^CLASS\d+", c)), keys(zep.header))), lt = sortclasses))
end

function NeXLSpectrum.elements(zep::Zepellin)
    sortelms(c1, c2) = isless(parse(Int, c1[5:end]), parse(Int, c2[5:end]))
    elms = (c -> get(zep.header, c, "")).(sort(collect(filter(c -> !isnothing(match(r"^ELEM\d+", c)), keys(zep.header))), lt = sortelms))
    return map(s -> PeriodicTable.elements[parse(Int, split(s, " ")[2])], elms)
end

function header(zep::Zepellin)
    keep(k) = !mapreduce(ty->startswith(k, uppercase(ty)), (a,b)->a||b, ("CLASS","ELEM","MAG"))
    kys = filter(keep, keys(zep.header))
    return DataStructures.SortedDict(filter(kv->kv[1] in kys, zep.header))
end

data(zep::Zepellin) = zep.data
