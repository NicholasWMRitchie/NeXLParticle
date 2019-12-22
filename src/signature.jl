using NeXLMatrixCorrection
using NeXLSpectrum
using Base.Threads

function askratios(ffr::FilterFitResult)::Vector{KRatio}
    res = KRatio[]
    for lbl in keys(ffr.kratios)
        if lbl isa NeXLSpectrum.CharXRayLabel
            push!(
                res,
                KRatio( #
                    lbl.xrays,
                    ffr.label.spec.properties,
                    lbl.spec.properties,
                    lbl.spec[:Composition],
                    NeXLUncertainties.uncertainvalue(lbl, ffr.kratios),
                ),
            )
        end
    end
    return res
end

"""
    Signature

Signature helps to define how the signature calculation is performed and to cache intermediate calculations for
speed.
"""
struct Signature
    kpure::Dict{CharXRay,<:AbstractFloat}
    matrix::Type{<:MatrixCorrection}
    fluorescence::Type{<:FluorescenceCorrection}
    kratioopt::KRatioOptimizer

    Signature(mc = XPPCorrection, fc = ReedFluorescence, kro = SimpleKRatioOptimizer(1.3)) =
        new(Dict{CharXRay,AbstractFloat}(), mc, fc, kro)
end

function kpure(sig::Signature, std::Material, elm::Element, lines::Vector{CharXRay}, stdE0, unkE0, stdθ, unkθ)::AbstractFloat
    res = get(sig.kpure, brightest(lines), -1.0)
    if res == -1.0
        unkzaf = ZAF(sig.matrix, sig.fluorescence, std, lines, stdE0)
        stdzaf = ZAF(sig.matrix, sig.fluorescence, pure(elm), lines, unkE0)
        res = k(unkzaf, stdzaf, stdθ, unkθ)
        sig.kpure[brightest(lines)] = res
    end
    return res
end

"""
    signature( #
        sig::Signature,
        krs::Vector{KRatio},
        special::Set{Element}
    )::Dict{Element,Float64}

Computes a "particle signature" from the specified set of k-ratios.  A particle signature
is similar to normalized k-ratios relative to pure elements except 1) certain elements
can be removed from the normalization; and 2) these elements can be include in the result as
un-normalized.
"""
function signature( #
    sig::Signature,
    krs::Vector{KRatio},
    special::Set{Element}
)::Dict{Element,<:AbstractFloat}
    kzs = Dict{Element,AbstractFloat}()
    #Scale the k-ratios relative to a pure element...
    for kr in optimizeks(sig.kratioopt, krs)
        # k_{unk,pure} = (I_{unk}/I_{std})*(I_{std}/I_{pure}) = kr.kratio*k(std,pure)
        kstdpure = isapprox(kr.standard[kr.element], 1.0, atol = 0.001) ? 1.0 : #
               kpure(
            sig,
            kr.standard,
            kr.element,
            kr.lines,
            kr.stdProps[:BeamEnergy],
            kr.unkProps[:BeamEnergy],
            kr.stdProps[:TakeOffAngle],
            kr.unkProps[:TakeOffAngle],
        )
        kzs[kr.element] = kr.kratio * kstdpure
    end
    onorm = sum(NeXLUncertainties.value(kzs[elm]) for elm in keys(kzs))
    res = Dict(elm => kzs[elm] / onorm for elm in special)
    notspecial = filter(elm -> !(elm in special), keys(kzs))
    norm = sum(NeXLUncertainties.value(kzs[elm]) for elm in notspecial)
    merge!(res, Dict(elm => kzs[elm] / norm for elm in notspecial))
    return res
end

abstract type CullingRule end

struct NSigmaCulling <: CullingRule
    nsigma::Float64
end

function cull(cr::NSigmaCulling, kr::KRatio)
    k, s = NeXLUncertainties.value(kr.kratio), NeXLUncertainties.σ(kr.kratio)
    return k > cr.nsigma * s ? kr : KRatio(kr.lines, kr.unkProps, kr.stdProps, kr.standard, UncertainValue(0.0, s))
end

function _quant(
    zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum},
    rows::Union{Vector{Int},UnitRange{Int}}=1:1000000;
    strip = Set{Element}( (n"C", ) ), # Elements to fit but not include in the result table.
    special = Set{Element}( ( n"O", ) ), # Element for special treatment in signature
    cullRule::CullingRule = NSigmaCulling(3.0),
    writeResidual::Bool = true,
    withUncertainty::Bool = true
)
    sortedbysig(newcols, row) =
         sort(collect(zip(newcols, row)),lt=(i1,i2)->!isless(i1[1]==n"O" ? 0.0 : i1[2],i2[1]==n"O" ? 0.0 : i2[2]))
    rows = ismissing(rows) ? eachparticle(zep) : rows
    # Build the filtered references
    filt = buildfilter(NeXLSpectrum.GaussianFilter, det)
    filtrefs = FilteredReference[]
    e0, toa = NaN64, NaN64
    # Not really worth threading
    for elm in keys(refs)
        lines = ktransitions
        if elm > n"Ca"
            lines = union(ktransitions, ltransitions)
        elseif elm > n"Ba"
            lines = union(ltransitions, mtransitions)
        end
        cfs = NeXLSpectrum.charFeature(elm, Tuple(lines), maxE = 0.9 * refs[elm][:BeamEnergy])
        fr = filter(refs[elm], det, cfs, filt, 1.0 / dose(refs[elm]))
        append!(filtrefs, fr)
        if isnan(e0)
            e0, toa = refs[elm][:BeamEnergy], refs[elm][:TakeOffAngle]
        end
    end
    # Create a list of columns (by element)
    newcols = sort(collect(filter(elm -> !(elm in strip), keys(refs))))
    quant = Array{Union{UncertainValue,Missing}}(missing, size(zep.data, 1), length(newcols))
    # quantify and tabulate each particle
    if writeResidual
        newdir = joinpath(dirname(zep.headerfile), "ResidualX")
        if !isdir(newdir)
            mkdir(newdir)
        end
    end
    evalrows = intersect(eachparticle(zep),rows)
    ignore = union(strip,special)
    krvs = Array{Union{Missing, Vector{KRatio}}}(missing, length(evalrows))
    counts = Array{Union{Missing,Float64}}(missing,length(evalrows))
    Threads.@threads for row in evalrows
        unk = spectrum(zep, row, false)
        if !ismissing(unk)
            # Particle spectra are often missing critical data items...
            unk[:ProbeCurrent], unk[:LiveTime] = get(unk, :ProbeCurrent, 1.0), get(unk, :LiveTime, 1.0)
            unk[:BeamEnergy], unk[:TakeOffAngle] = get(unk, :BeamEnergy, e0), get(unk, :TakeOffAngle, toa)
            # Fit, cull and then compute the particle signature...
            res = fit(FilteredUnknownW, unk, filt, filtrefs, true)
            counts[row] = NeXLSpectrum.characteristiccounts(res, ignore)
            if writeResidual
                filename = joinpath(dirname(zep.headerfile), "Residual", "00000$(row)"[end-4:end] * ".msa")
                writeEMSA(filename, NeXLSpectrum.residual(res))
            end
            culled = map(kr -> cull(cullRule, kr), askratios(res))
            krvs[row] = filter(kr -> kr.element in newcols, culled)
        end
    end
    nfirst=min(4, length(newcols)-length(special))
    felm = Array{Union{Element,Missing}}(missing, size(zep.data, 1), 4)
    fsig = Array{Union{UncertainValue,Missing}}(missing, size(zep.data, 1), 4)
    sigx = Signature(XPPCorrection, ReedFluorescence, SimpleKRatioOptimizer(1.3))
    for row in evalrows
        krv = krvs[row]
        if !ismissing(krv)
            sig = signature(sigx, krv, special) # Not thread safe...
            quant[row, :] = map(elm -> 100.0 * sig[elm], newcols)
            s2 = copy(sig)
            foreach(i->delete!(s2,i), special) # Special can't be a FIRSTELM, ...
            firstElms = topN(s2, nfirst)
            felm[row, :] = map(i->i<=nfirst ? firstElms[i].first : missing, 1:4)
            fsig[row, :] = map(i->i<=nfirst ? 100.0 * firstElms[i].second : missing, 1:4)
        end
    end
    fourelms = DataFrame( #
        FIRSTELM= felm[:,1], FIRSTPCT = NeXLUncertainties.value.(fsig[:,1]),
        SECONDELM = felm[:,2], SECONDPCT = NeXLUncertainties.value.(fsig[:,2]), #
        THIRDELM = felm[:,3],  THIRDPCT = NeXLUncertainties.value.(fsig[:,3]), #
        FOURTHELM = felm[:,4], FOURTHPCT = NeXLUncertainties.value.(fsig[:,4]), #
        COUNTS = counts) #
        # Return the uncertainties or not...
    cols = Pair{Symbol,AbstractVector{Float64}}[ ]
    for (i, elm) in enumerate(newcols)
        push!(cols, convert(Symbol,elm)=>NeXLUncertainties.value.(quant[:,i]))
        if withUncertainty
            push!(cols, Symbol("U[$(elm.symbol)]")=>NeXLUncertainties.σ.(quant[:,i]))
        end
    end
    quantRes = DataFrame(cols...)
    return hcat(fourelms, quantRes)
end

function quantify(zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum},
    rows::Union{Vector{Int},UnitRange{Int}}=1:1000000;
    strip = Set{Element}( (n"C", ) ), # Elements to fit but not include in the result table.
    special = Set{Element}( ( n"O", ) ), # Element for special treatment in signature
    cullRule::CullingRule = NSigmaCulling(3.0),
    writeResidual::Bool = true,
    withUncertainty::Bool = true)
    qr = _quant(zep,det,refs,rows,strip=strip,special=special,cullRule=cullRule,writeResidual=writeResidual,withUncertainty=withUncertainty)
    # Remove old items...
    removeme = map(elm -> convert(Symbol, elm), zep.elms)
    append!(removeme, map(elm -> Symbol("U[$elm.symbol]"), zep.elms))
    append!(removeme, ALL_COMPOSITIONAL_COLUMNS)
    # append!(removeme, ALL_CLASS_COLUMNS)
    remaining = copy(zep.data[:, filter(f -> !(f in removeme), names(zep.data))])
    data = hcat(remaining, qr)
    # Replace outdated header items
    header = copy(zep.header)
    headerfile =  joinpath(dirname(zep.headerfile),"generated.hdz") # So spectra and other path related items continue to work
    foreach(f->if startswith(f,"ELEM") delete!(header,f) end, keys(header))
    delete!(header, "MAX_PARTICLE")
    delete!(header, "MAX_RESIDUAL")
    # Add/replace header items
    for i in 1:1000
        if !haskey(header, "ANCESTOR[$i]")
            header["ANCESTOR[$i]"]=zep.headerfile
            break
        end
    end
    elmCx = 0
    for elm in sort([keys(refs)...])
        if !(elm in strip)
            elmCx+=1
            header["ELEM[$elmCx]"] = "$(elm.symbol) $(z(elm)) 1"
        end
    end
    header["ELEMENTS"] = "$(elmCx)"
    header["DATAFILES"] = replace(headerfile,r".[h|H][d|D][z|Z]$"=>".*")
    header["VEC_FILE"] = "NeXLParticle"
    return Zeppelin(headerfile, header, data)
end

topN(items::Dict, n::Int) =
    sort( [ items... ],lt=(a1,a2)->isless(a2[2],a1[2]))[1:min(n,length(items))]
