using CSV
using DataFrames
using PeriodicTable

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
        "PSEM_CLASS" => "CLASS"
    )
    columnnames(zep::Zepellin) = map(cn->get(remapcolumnnames, cn, cn), map(c->c[1], zep.columns))
    header, columns, hdr = Dict{String,String}(), [], true
    for line in readlines(headerfilename)
        if hdr
            (k, v) = strip.(split(line,"="))
            header[k] = v
            hdr = !isequal(uppercase(k),"PARTICLE_PARAMETERS")
        else
            push!(columns, string.(strip.(split(line,"\t"))))
        end
    end
    pxz = CSV.File(replace(headerfile,".hdz"=>".pxz"), header=columnnames(columns), normalizenames=true) |> DataFrame
    categorical!(pxz, :CLASS) # Convert class column to pxz
    return ( header, columns, pxz )
end

struct Zepellin
    headerfile::String
    header::Dict{String,String}
    columns::Vector{NTuple{3, String}}
    data::DataFrame

    function Zepellin(headerfilename::String)
        new(headerfilename, loadZep(headerfilename)...)
    end
end

function classnames(zep::Zepellin)
    sortclasses(c1, c2) = isless(parse(Int, c1[6:end]), parse(Int,c2[6:end]))
    return (c->get(header, c,"")).(sort(collect(filter(c->!isnothing(match(r"^CLASS\d+",c)),keys(zep.header))),lt=sortclasses))
end

function elements(zep::Zepellin)
    sortelms(c1, c2) = isless(parse(Int, c1[5:end]), parse(Int,c2[5:end]))
    elms = (c->get(header, c,"")).(sort(collect(filter(c->!isnothing(match(r"^ELEM\d+",c)),keys(zep.header))),lt=sortelms))
    return map(s->elements[parse(Int,split(s, " ")[2])],elms)
end

data(zep::Zepellin) = zep.data
