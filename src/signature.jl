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
    special::Vector{Element} = Vector[n"O",],
    mc::Type{<:MatrixCorrection} = NeXLMatrixCorrection.XPPCorrection,
    fc::Type{<:FluorescenceCorrection} = ReedFluorescence,
)::Dict{Element,<:AbstractFloat}
    kzs = Dict{Element,AbstractFloat}()
    #Scale the k-ratios relative to a pure element...
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
            kr.unkProps[:BeamEnergy],
            kr.stdProps[:TakeOffAngle],
            kr.unkProps[:TakeOffAngle], #
            mc,
            fc,
        )
        kzs[kr.element] = kr.kratio * kstd
    end
    res = Dict{Element,AbstractFloat}()
    onorm = sum(NeXLUncertainties.value(kzs[elm]) for elm in keys(kzs))
    res = Dict(elm => (kzs[elm] / onorm) for elm in special)
    notspecial = filter(elm -> !(elm in special), keys(kzs))
    norm = sum(NeXLUncertainties.value(kzs[elm]) for elm in notspecial)
    merge!(res, Dict(elm => (kzs[elm] / norm) for elm in notspecial))
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

function quantify(
    zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum},
    rows::Union{Vector{Int},UnitRange{Int}}=1:1000000;
    strip::Vector{Element} = [n"C"], # Element to not include in table
    special::Vector{Element} = [n"O"], # Element for special treatment in signature
    cullRule::CullingRule = NSigmaCulling(3.0),
)
    rows = ismissing(rows) ? eachparticle(zep) : rows
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
        cfs = NeXLSpectrum.charFeature(elm, Tuple(lines), maxE = 0.9 * refs[elm][:BeamEnergy])
        fr = filter(refs[elm], det, cfs, filt, 1.0 / dose(refs[elm]))
        append!(filtrefs, fr)
        e0, toa = refs[elm][:BeamEnergy], refs[elm][:TakeOffAngle]
    end
    # Create a list of columns (by element)
    newcols = sort(collect(filter(elm -> !(elm in strip), keys(refs))))
    quant = zeros(UncertainValue, size(zep.data, 1), length(newcols))
    # quantify and tabulate each particle
    for row in intersect(eachparticle(zep),rows)
        unk = spectrum(zep, row, false)
        if !ismissing(unk)
            # Particle spectra are often missing critical data items...
            unk[:ProbeCurrent], unk[:LiveTime] = get(unk, :ProbeCurrent, 1.0), get(unk, :LiveTime, 1.0)
            unk[:BeamEnergy], unk[:TakeOffAngle] = get(unk, :BeamEnergy, e0), get(unk, :TakeOffAngle, toa)
            # Fit, cull and then compute the particle signature...
            res = fit(FilteredUnknownW, unk, filt, filtrefs, true)
            culled = map(kr -> cull(cullRule, kr), askratios(res))
            krv = filter(kr -> !(kr.element in strip), culled)
            sig = signature(krv, special)
            quant[row, :] =  collect( 100.0 * sig[elm] for elm in newcols )
        else
            quant[row, :] = collect( missing for _ in eachindex(newcols) )
        end
    end
    quantRes = DataFrame((convert(Symbol, elm) => quant[:, i] for (i, elm) in enumerate(newcols))...)
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
