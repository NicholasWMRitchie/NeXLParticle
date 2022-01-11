using CategoricalArrays

struct OrderedRuleSet <: ParticleClassifier
    name::String
    rules::Vector{Tuple{String,Function}}

    OrderedRuleSet(name::String, rules::Tuple{String,Function}...) = new(name, collect(rules))

    function OrderedRuleSet(orss::OrderedRuleSet...)
        rules = Vector{Tuple{String,Function}}()
        name = String[]
        for ors in orss
            append!(rules, ors.rules)
            push!(name, ors.name)
        end
        return new(join(name, "+"), rules)
    end
end

Base.show(io::IO, ors::OrderedRuleSet) = print(io, "JORS[$(ors.name)]")
classnames(sr::OrderedRuleSet) = collect(map(r -> r[1], sr.rules))

"""
    classify(zep::Zeppelin, ruleset::AbstractString, classnames::AbstractVector, clsidx::AbstractVector{Int})

Constructs a new Zeppelin item and replaces the "CLASS" column with the classes named in `classnames` and indexed into by `clsidx`.
"""
function classify(zep::Zeppelin, ruleset::AbstractString, classnames::AbstractVector{<:AbstractString}, clsidx::AbstractVector{<:Integer})
    @assert length(clsidx)==nrow(zep.data) "The number of assigned classes must match the number of rows."
    @assert minimum(clsidx)>=-1 "The minimum class index must be -1 or greater."
    @assert maximum(clsidx)<length(classnames) "The maximum class index must less than the number of class names."
    clscol = findfirst(isequal("CLASS"), names(zep.data))
    # Remove "CLASS" and "VERIFIED_CLASS" columns
    rd = copy(zep.data)[:, Cols(x -> !(x in ("CLASS", "VERIFIEDCLASS")))]
    result = Zeppelin("$(zep.headerfile[1:end-4])[$ruleset].hdz", copy(zep.header), rd, classnames)
    if isnothing(clscol)
        tmp = findfirst(isequal("FIRSTELM"), rd)
        clscol = isnothing(tmp) ? ncol(rd) + 1 : tmp
    end
    insertcols!(result.data, clscol, "CLASS" => map(idx->ZepClass(result, idx-1), clsidx))
    insertcols!(result.data, clscol+1, "VERIFIEDCLASS" => fill(ZepClass(result, -1), nrow(rd)))
    return result
end

"""
    classify(zep::Zeppelin, ruleset::AbstractString, catcls::AbstractVector{String})

Create a new Zeppelin file and assign the contents of this CategoricalArray to the "CLASS" column.
"""
classify(zep::Zeppelin, ruleset::AbstractString, catcls::CategoricalVector{String}) =
    classify(zep, ruleset, levels(catcls), catcls.refs)

struct _RowAdapter{T <: DataFrames.DataFrameRow}
    row::T
    default::Dict{Symbol,Any}

    function _RowAdapter(r::DataFrames.DataFrameRow, def::Dict{Symbol,<:Any})
        new{typeof(r)}(r, def)
    end
end

Base.getindex(rw::_RowAdapter, key::Symbol) = #
    get(rw.row, key) do 
        getindex(rw.default, key)
    end

Base.get(rw::_RowAdapter, key::Symbol, def) = #
    get(rw.row, key) do 
        get(rw.default, key, def)
    end

function Base.getproperty(rw::_RowAdapter, key::Symbol)
    if key in ( :row, :default ) 
        getfield(rw, key)
    else
        get(getfield(rw, :row), key) do 
            getindex(getfield(rw, :default), key)
        end
    end
end

Base.sum(rw::_RowAdapter, keys::Symbol...) = sum(key->rw[key], keys)

between(val::AbstractFloat, min::AbstractFloat, max::AbstractFloat) = (val>=min) && (val<=max)
between(rw::_RowAdapter, key::Symbol, min::AbstractFloat, max::AbstractFloat) = between(rw[key], min, max)

"""
    classify(zep::Zeppelin, sr::OrderedRuleSet)::Zeppelin

Classify the rows in `zep`, row-by-row by calling the applying the `OrderedRuleSet` to each row. Returns a new
Zeppelin object.
"""
function classify(zep::Zeppelin, sr::OrderedRuleSet)::Zeppelin
    cat = empty!(categorical(append!(classnames(sr), ["Other", "Missing", "Error"]), compress=false, ordered = true)) # set up the levels
    errcx = 0
    defs = Dict(Symbol(uppercase(e.symbol))=>0.0 for e in elements)
    cres = ThreadsX.map(eachrow(zep.data)) do row
        try
            input = _RowAdapter(row, defs)
            i = findfirst(rule -> rule[2](input), sr.rules)
            isnothing(i) ? "Other" : sr.rules[i][1] 
        catch err
            if (errcx += 1) <= 4
                @info "ERROR[$errcx]: " exception=(err, catch_backtrace())
            end
            "Error"
        end
    end
    append!(cat, cres)
    res = classify(zep, repr(sr), cat)
    res.header["RULE_FILE"]=repr(sr)
    return res
end
