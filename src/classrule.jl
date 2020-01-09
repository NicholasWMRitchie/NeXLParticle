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
    classify(zep::Zeppelin, catcls::CategoricalArray{String})

Assign the CLASS column to the contents of this CategoricalArray.
"""
function classify(zep::Zeppelin, catcls::CategoricalArray{String})::Zeppelin
    @assert length(catcls) == nrow(zep.data)
    # Figure out where to put the column...
    cols = collect(1:ncol(zep.data))
    clscol = findfirst(n -> n == :CLASS, names(zep.data))
    if !isnothing(clscol)
        cols = deleteat!(cols, clscol)
    else
        clscols = findfirst(n -> n == :FIRSTELM, names(zep.data))
        clscols = isnothing(clscols) ? ncols(zep.data) + 1 : clscols
    end
    # Copy the data and insert the column
    result = Zeppelin(zep.headerfile, zep.header, zep.data[:, cols])
    insertcols!(result.data, clscol, :CLASS => catcls)
    return result
end

"""
    classify(zep::Zeppelin, )::Zeppelin

Classify the rows in `zep`, row-by-row by calling the applying the `OrderedRuleSet` to each row. Returns a new
Zeppelin object.
"""
function classify(zep::Zeppelin, sr::OrderedRuleSet)::Zeppelin
    cat = categorical(append!(classnames(sr), ["Other", "Missing"]), false, ordered = true) # set up the levels
    cat = cat[1:0] # delete the data but retain the levels
    errcx = 0
    missingElms = Dict(s => 0.0 for s in filter(
        sy -> !(sy in names(zep.data)),
        map(z -> convert(Symbol, elements[z]), 1:95),
    ))
    for row in eachparticle(zep)
        try
            input = Dict(n => zep.data[row, n] for n in names(zep.data))
            if all(map(v->!ismissing(v),values(input)))
                merge!(input, missingElms)
                push!(cat, classify(sr, input))
            else
                push!(cat, "Missing")
            end
        catch err
            push!(cat, "Other")
            if (errcx += 1) <= 4
                print("ERROR[$errcx]: ")
                showerror(stdout, err)
                println()
                for st in stacktrace(catch_backtrace())[1:3]
                    println(st)
                end
            end
        end
    end
    res = classify(zep, cat)
    res.header["RULE_FILE"]=repr(sr)
    return res
end

function classify(sr::OrderedRuleSet, input::Dict{Symbol,Any})
    i = findfirst(rule -> rule[2](input), sr.rules)
    return isnothing(i) ? "Other" : sr.rules[i][1]
end

const NullRules = OrderedRuleSet("Null", ("Unclassified", inp -> true) )

const BaseRules = OrderedRuleSet(
    "Base",
    ("LowCounts", inp -> inp[:COUNTS] < 2000),
    (
     "Maraging",
     inp -> (inp[:FIRSTELM] == n"Fe") &&
            (inp[:NI] > 5) && (inp[:CO] > 2) && (inp[:FE] + inp[:NI] + inp[:CO] + inp[:MO] > 80.0),
    ),
    ("Quartz-like", inp -> inp[:SI] > 80),
    (
     "Biotite",
     inp -> (inp[:K] > 5) &&
            (inp[:MG] + inp[:FE] > 20) &&
            (inp[:AL] > 10) &&
            (inp[:SI] > 10) && (inp[:K] + inp[:MG] + inp[:FE] + inp[:AL] + inp[:SI] > 80) && (inp[:O] > 5),
    ),
    (
     "Plagioclaise",
     inp -> (inp[:CA] + inp[:NA] > 20) &&
            (inp[:AL] > 10) && (inp[:SI] > 10) && (inp[:NA] + inp[:CA] + inp[:AL] + inp[:SI] > 80) && (inp[:O] > 5),
    ),
    (
     "Othoclaise",
     inp -> (inp[:K] > 5) &&
            (inp[:AL] > 20) && (inp[:SI] > 20) && (inp[:O] > 10) && (inp[:K] + inp[:AL] + inp[:SI] > 80),
    ),
    ("Dolomite", inp -> (inp[:CA] > 40) && (inp[:MG] > 10) && (inp[:CA] + inp[:MG] > 70)),
    ("Olivine", inp -> (inp[:MG] + inp[:FE] > 40) && (inp[:SI] > 20) && (inp[:MG] + inp[:FE] + inp[:SI] > 80)),
    (
     "Feldspar",
     inp -> (inp[:K] + inp[:NA] + inp[:CA] > 20) &&
            (inp[:AL] > 10) && (inp[:SI] > 20) && (inp[:K] + inp[:NA] + inp[:CA] + inp[:AL] + inp[:SI] > 80),
    ),
    ("Willemite", inp -> (inp[:ZN] > 20) && (inp[:SI] > 20) && (inp[:ZN] + inp[:SI] > 80)),
    (
     "Garnet",
     inp -> (inp[:MG] + inp[:FE] + inp[:MN] + inp[:CA] > 10) &&
            (inp[:AL] > 5) && (inp[:SI] > 10) && (inp[:MG] + inp[:FE] + inp[:MN] + inp[:CA] + inp[:AL] + inp[:SI] > 80),
    ),
    (
     "Al-Cu alloy (2XXX)",
     inp -> (inp[:AL] > 60) && (inp[:SECONDELM] == n"Cu") && (inp[:CU] > 3) && (inp[:AL] + inp[:CU] > 90),
    ),
    (
     "Al-Mn alloy (3XXX)",
     inp -> (inp[:AL] > 60) && (inp[:SECONDELM] == n"Mn") && (inp[:MN] > 3) && (inp[:AL] + inp[:MN] > 90),
    ),
    (
     "Al-Si alloy (4XXX)",
     inp -> (inp[:AL] > 80) && (inp[:SECONDELM] == n"Si") && (inp[:SI] > 1) && (inp[:AL] + inp[:SI] > 90),
    ),
    (
     "Al-Mg alloy (5XXX)",
     inp -> (inp[:FIRSTELM] == n"Al") && (inp[:SECONDELM] == n"Mg") && (inp[:MG] > 5) && (inp[:AL] + inp[:MG] > 80),
    ),
    (
     "Al-Si-Mg (6XXX)",
     inp -> (inp[:AL] > 90) &&
            (inp[:SI] > 0) && (inp[:MG] > 0) && (inp[:AL] + inp[:SI] + inp[:MG] + inp[:MN] + inp[:CU] + inp[:FE] > 98),
    ),
    (
     "Al-Zn alloy (7XXX)",
     inp -> (inp[:AL] > 80) && (inp[:SECONDELM] == n"Zn") && (inp[:ZN] > 2) && (inp[:AL] + inp[:ZN] > 90),
    ),
    ("Al-Fe alloy", inp -> (inp[:AL] > 80) && (inp[:SECONDELM] == n"Fe") && (inp[:FE] > 2) && (inp[:AL] + inp[:FE] > 90)),
    ("Aluminosilicate", inp -> (inp[:AL] > 30) && (inp[:SI] > 10) && (inp[:AL] + inp[:SI] > 80)),
    ("Zircon", inp -> (inp[:ZR] > 20) && (inp[:ZR] + inp[:SI] > 80)),
    (
     "Elgiloy",
     inp -> (inp[:CO] > 20) &&
            (inp[:CR] > 10) &&
            (inp[:NI] > 5) &&
            (inp[:FE] > 5) &&
            ((inp[:MO] > 2) || (inp[:S] > 2)) &&
            (inp[:CO] + inp[:CR] + inp[:NI] + inp[:FE] + inp[:MO] + inp[:MN] + inp[:S] > 80),
    ),
    ("Aluminum", inp -> inp[:AL] > 80),
    ("Fe-Al alloy", inp -> (inp[:AL] > 20) && (inp[:FE] > 20) && (inp[:AL] + inp[:FE] > 80)),
    ("Calcium silicate", inp -> (inp[:CA] > 40) && (inp[:SI] > 10) && (inp[:CA] + inp[:SI] > 80)),
    ("Fe-Cr stainless", inp -> (inp[:FE] > 50) && (inp[:CR] > 5) && (inp[:NI] < 3) && (inp[:FE] + inp[:CR] > 90)),
    (
     "Fe-Cr-Ni Stainless",
     inp -> (inp[:FIRSTELM] == n"Fe") && (inp[:CR] > 5) && (inp[:NI] > 5) && (inp[:FE] + inp[:CR] + inp[:NI] > 80),
    ),
    ("Fe-Ni stainless", inp -> (inp[:FE] > 50) && (inp[:NI] > 10) && (inp[:FE] + inp[:NI] > 90)),
    ("Iron oxide", inp -> (inp[:FE] > 80) && (inp[:O] > 10)),
    ("Iron", inp -> inp[:FE] > 80),
    ("Calcite", inp -> (inp[:CA] > 80) && (inp[:O] > 5)),
    ("Zinc", inp -> inp[:ZN] > 80),
    ("Silicate", inp -> inp[:FIRSTELM] == n"Si" ),
    (
     "Anglesite",
     inp -> (inp[:FIRSTELM] == n"Pb") && (inp[:SECONDELM] == n"S") && (inp[:PB] + inp[:S] > 80) && (inp[:O] > 5),
    ),
    ("Hashemite", inp -> (inp[:BA] > 20) && (inp[:CR] > 20) && (inp[:BA] + inp[:CR] > 80) && (inp[:O] > 5)),
    ("Barite", inp -> (inp[:BA] > 40) && (inp[:S] > 10) && (inp[:BA] + inp[:S] > 80) && (inp[:O] > 5)),
    ("Gypsum", inp -> (inp[:CA] > 20) && (inp[:S] > 20) && (inp[:CA] + inp[:S] > 80) && (inp[:O] > 5)),
    ("Lead", inp -> (inp[:PB] > 60) && (inp[:PB] + inp[:S] > 80)),
    ("Salt", inp -> (inp[:K] + inp[:NA] > 20) && (inp[:CL] > 20) && (inp[:K] + inp[:NA] + inp[:CL] > 80)),
    ("Bismuth", inp -> (inp[:BI] > 80)),
    ("Brass", inp -> (inp[:CU] > 20) && (inp[:ZN] > 20) && (inp[:CU] + inp[:ZN] > 80)),
    ("Bronze", inp -> (inp[:CU] > 20) && (inp[:SN] > 20) && (inp[:CU] + inp[:SN] > 80)),
    ("Cupronickel", inp -> (inp[:CU] > 30) && (inp[:NI] > 10) && (inp[:CU] + inp[:NI] > 80)),
    ("Sodium sulfide", inp -> (inp[:NA] > 20) && (inp[:S] > 10) && (inp[:NA] + inp[:S] > 80)),
    ("Wurtzite (ZnS)", inp -> (inp[:ZN] > 20) && (inp[:S] > 10) && (inp[:ZN] + inp[:S] > 80)),
    ("Pyrite (FeS2)", inp -> (inp[:FE] > 20) && (inp[:S] > 20) && (inp[:FE] + inp[:S] > 80)),
    ("Silver sulfide", inp -> (inp[:AG] > 20) && (inp[:S] > 10) && (inp[:AG] + inp[:S] > 80)),
    (
     "Sn62 solder",
     inp -> (inp[:SN] > 30) && (inp[:PB] > 10) && (inp[:AG] > 1) && (inp[:SN] + inp[:PB] + inp[:AG] > 80),
    ),
    (
     "SACZ solder",
     inp -> (inp[:SN] + inp[:AG] + inp[:CU] + inp[:ZN] > 80) &&
            (inp[:SN] > 10) && (inp[:CU] > 10) && (inp[:AG] > 10) && (inp[:ZN] > 2),
    ),
    ("Lead solder", inp -> (inp[:PB] > 20) && (inp[:SN] > 20) && (inp[:PB] + inp[:SN] > 80)),
    (
     "SAC solder",
     inp -> (inp[:SN] + inp[:AG] + inp[:CU] > 80) && (inp[:SN] > 10) && (inp[:CU] > 10) && (inp[:AG] > 10),
    ),
    ("Titanium", inp -> (inp[:TI] > 80)),
    ("Calcium phosphate", inp -> (inp[:CA] > 40) && (inp[:P] > 10) && (inp[:CA] + inp[:P] > 80)),
    ("Calcium fluoride", inp -> (inp[:CA] > 40) && (inp[:F] > 10) && (inp[:CA] + inp[:F] > 80)),
    ("Iron-bearing", inp -> (inp[:FE] > 50)),
    ("Ca-bearing", inp -> inp[:CA] > 50),
    ("Other Silicate", inp -> (inp[:SI] > 20) && (inp[:AL] + inp[:FE] + inp[:CA] + inp[:SI] > 80)),
    ("Gold-Silver alloy", inp -> (inp[:FIRSTELM] == n"Au") && (inp[:SECONDELM] == n"Ag") && (inp[:AU] + inp[:AG] > 80)),
    ("Gold - other", inp -> inp[:AU] > 80),
    ("Tin", inp -> inp[:SN] > 80),
    ("Celestine", inp -> (inp[:SR] > 40) && (inp[:S] > 10) && (inp[:SR] + inp[:S] > 80) && (inp[:O] > 5)),
    ("Barite-Celestine", inp -> (inp[:BA] + inp[:SR] > 40) && (inp[:BA] + inp[:SR] + inp[:S] > 80) && (inp[:O] > 5)),
    ("Cu-rich", inp -> inp[:CU] > 80),
    ("Cr-bearing", inp -> inp[:CR] > 50),
    ("Silver sulfide", inp -> (inp[:AG] > 40) && (inp[:S] > 10) && (inp[:AG] + inp[:S] > 80)),
    ("Silver", inp -> inp[:AG] > 80),
    ("MoS2", inp -> inp[:MO] > 40 && inp[:S] > 30 && (inp[:MO] + inp[:S] > 90)),
    (
     "Monazite",
     inp -> (inp[:LA] + inp[:CE] + inp[:ND] + inp[:TH] + inp[:Y] > 30) &&
            (inp[:P] > 10) && (inp[:LA] + inp[:CE] + inp[:ND] + inp[:TH] + inp[:Y] + inp[:P] + inp[:SI] > 80),
    ),
    ("Lanthanide", inp -> (inp[:LA] + inp[:CE] + inp[:ND] + inp[:Y] > 60)),
    ("Silver chloride", inp -> (inp[:AG] > 40) && (inp[:CL] > 10)),
    ("Magnesium oxide", inp -> (inp[:MG] > 80) && (inp[:O] > 10)),
    ("Nickel", inp -> inp[:NI] > 80),
    ("Ca-Fe silicate", inp -> (inp[:CA] > 20) && (inp[:FE] > 20) && (inp[:SI] > 5)),
    ("Ca-Ti silicate", inp -> (inp[:CA] > 20) && (inp[:TI] > 20) && (inp[:SI] > 10)),
    ("Al-Zr-Cl", inp -> (inp[:AL] > 20) && (inp[:ZR] > 20) && (inp[:CL] > 5)),
    ("Chlorine-rich", inp -> inp[:CL] > 80),
    ("Au-Cu alloy", inp -> (inp[:AU] > 40) && (inp[:CU] > 10) && (inp[:AU] + inp[:CU] + inp[:NI] > 80)),
    ("Uranium-bearing", inp -> inp[:U] > 20),
    ("Thorium-bearing", inp -> inp[:TH] > 20),
    ("Platinum", inp -> inp[:PT] > 80)
);
