using NeXLMatrixCorrection

using Memoization  # Implements caching of k-ratio relative to pure....

@memoize function kpure( #
    std::Material,
    elm::Element,
    lines::Vector{CharXRay},
    stdE0,
    unkE0,
    stdθ,
    unkθ,
    mc::Type{<:MatrixCorrection},
    fc::Type{<:FluorescenceCorrection},
)::AbstractFloat
    unkzaf = ZAF(mc, fc, std, lines, stdE0)
    stdzaf = ZAF(mc, fc, pure(elm), lines, unkE0)
    kk = k(unkzaf, stdzaf, stdθ, unkθ)
    @info "k[$(elm.symbol) in $(name(std)), $(name(lines))] = $(kk)"
    return kk
end

function askratios(ffr::FilterFitResult)::Vector{KRatio}
    res = KRatio[]
    for lbl in keys(ffr.kratios)
        if lbl isa NeXLSpectrum.CharXRayLabel
            push!(res, KRatio( #
                lbl.xrays,
                ffr.label.spec.properties,
                lbl.spec.properties,
                lbl.spec[:Composition],
                NeXLUncertainties.uncertainvalue(lbl, ffr.kratios),
            ))
        end
    end
    return res
end

"""
    signature( #
        krs::Array{KRatio},
        special::Array{Element} = [n"O"];
        mc::Type{MatrixCorrection} = XPPCorrection,
        fc::Type{FluorescenceCorrection} = ReedFluorescence

    )::Dict{Element,Float64}

Computes a "particle signature" from the specified set of k-ratios.  A particle signature
is similar to normalized k-ratios relative to pure elements except 1) certain elements
can be removed from the normalization; and 2) these elements can be include in the result as
un-normalized.
"""
function signature( #
    krs::Vector{KRatio},
    special::Vector{Element} = Vector[ n"O", ],
    mc::Type{<:MatrixCorrection} = NeXLMatrixCorrection.XPPCorrection,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence
)::Dict{Element,<:AbstractFloat}
    kzs = Dict{Element,AbstractFloat}()
    for kr in optimizeks(SimpleKRatioOptimizer(1.3), krs)
        # k = (Iunk/Istd)/(Istd/Iz) = (Iunk/Iz) = kr.kratio/k(unk,std)
        kstd = isapprox(kr.standard[kr.element], 1.0, atol = 0.001) ? #
               1.0 / kr.standard[kr.element] : #
            # Can't hand in full kr or memoization fails (since kr.kratio changes...)
               kpure(
            kr.standard,
            kr.element,
            kr.lines, #
            kr.stdProps[:BeamEnergy],
            haskey(kr.unkProps, :BeamEnergy) ? kr.unkProps[:BeamEnergy] : kr.stdProps[:BeamEnergy], #
            kr.stdProps[:TakeOffAngle],
            haskey(kr.unkProps, :TakeOffAngle) ? kr.unkProps[:TakeOffAngle] : kr.stdProps[:TakeOffAngle], #
            mc,
            fc,
        )
        kzs[kr.element] = kr.kratio / kstd
    end
    res = Dict{Element,AbstractFloat}()
    notspecial = filter(kr -> !(kr.element in special), krs)
    norm = sum(NeXLUncertainties.value(kr.kratio) for kr in notspecial)
    res = Dict(kr.element => kr.kratio / norm for kr in notspecial)
    merge!(res, Dict(kr.element => kr.kratio for kr in filter(kr -> kr.element in special, krs)))
    return res
end

abstract type CullingRule end

struct NSigmaCulling
    nsigma::Float64
end

function cull(cr::NSigmaCulling, kr::KRatio)
    k, s = NeXLUncertainties.value(kr.kratio), NeXLUncertainties.σ(kr.kratio)
    return k > cr.nsigma*s ? kr :
        KRatio(kr.lines, kr.unkProps, kr.stdProps, kr.standard, UncertainValue(0.0,s))
end

function quantify(
    zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum},
    strip::Vector{Element} = [n"C"], # Element to not include in table
    special::Vector{Element} = [n"O"], # Element for special treatment in signature
    cullRule::NSigmaCulling = NSigmaCulling(3.0)
)
    # Build the filtered references
    filt = buildfilter(NeXLSpectrum.GaussianFilter, det)
    filtrefs = FilteredReference[]
    e0, toa = NaN64, NaN64
    for elm in keys(refs)
        lines = ktransitions
        if elm > n"Ca"
            lines = union(ktransitions, ltransitions)
        elseif elm > n"Ba"
            lines = union(ltransitions, mtransitions)
        end
        cfs = NeXLSpectrum.charFeature(elm, Tuple(lines), maxE=0.9*refs[elm][:BeamEnergy])
        fr = filter(refs[elm], det, cfs, filt, 1.0 / dose(refs[elm]))
        append!(filtrefs,fr)
        e0, toa = refs[elm][:BeamEnergy], refs[elm][:TakeOffAngle]
    end
    # Create a list of columns (by element)
    newcols = sort(collect(filter(z -> !(z in strip), keys(refs))))
    quant = zeros(UncertainValue, size(zep.data,1), length(newcols))
    # quantify and tabulate each particle
    for row in eachparticle(zep)
        unk = spectrum(zep, row, false)
        if !ismissing(unk)
            # Particle spectra are often missing critical data items...
            unk[:ProbeCurrent], unk[:LiveTime] = get(unk, :ProbeCurrent, 1.0), get(unk, :LiveTime, 1.0)
            unk[:BeamEnergy] = get(unk,:BeamEnergy, e0)
            unk[:TakeOffAngle] = get(unk,:TakeOffAngle, toa)
            res = fit(FilteredUnknownG, unk, filt, filtrefs, false)
            culled = map(kr->cull(cullRule, kr), askratios(res))
            krv = filter(kr->!(kr.element in strip), culled)
            sig = signature(krv, special)
            for (i, col) in enumerate(newcols)
                quant[row,i]=100.0*sig[col]
            end
        else
            for (i, col) in enumerate(newcols)
                quant[row,i]=missing
            end
        end
    end
    quantRes = DataFrame( (convert(Symbol,elm)=> quant[:,i] for (i, elm) in enumerate(newcols))...)
    # First, strip out all old elemental data items
    if false
        removeme = map(elm -> convert(Symbol, elm), zep.elms)
        append!(removeme, COMPOSITIONAL_COLUMNS)
        append!(removeme, CLASS_COLUMNS)
        remaining = copy(zep.data[:, filter(f -> !(f in removeme), names(zep.data))])
        return hcat(remaining, quantRes)
    else
        return quantRes
    end
end
