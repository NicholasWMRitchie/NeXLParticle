using NeXLMatrixCorrection
using NeXLSpectrum

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
        special::Vector{Element}
    )::Dict{Element,Float64}

Computes a "particle signature" from the specified set of k-ratios.  A particle signature
is similar to normalized k-ratios relative to pure elements except 1) certain elements
can be removed from the normalization; and 2) these elements can be include in the result as
un-normalized.
"""
function signature( #
    sig::Signature,
    krs::Vector{KRatio},
    special::Vector{Element}
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
    @assert isequal(special[1],n"O") "$(special[1])!=$(n"O")"
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

function quantify(
    zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum},
    rows::Union{Vector{Int},UnitRange{Int}}=1:1000000;
    strip::Vector{Element} = Element[n"C"], # Element to not include in table
    special::Vector{Element} = Element[n"O"], # Element for special treatment in signature
    cullRule::CullingRule = NSigmaCulling(3.0),
)
    sortedbysig(newcols, row) =
         sort(collect(zip(newcols, row)),lt=(i1,i2)->!isless(i1[1]==n"O" ? 0.0 : i1[2],i2[1]==n"O" ? 0.0 : i2[2]))
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
    quant = Array{Union{UncertainValue,Missing}}(missing, size(zep.data, 1), length(newcols))
    felm = Array{Union{Element,Missing}}(missing, size(zep.data, 1), min(4, length(newcols)))
    fsig = Array{Union{UncertainValue,Missing}}(missing, size(zep.data, 1), min(4, length(newcols)))
    sigx = Signature(XPPCorrection, ReedFluorescence, SimpleKRatioOptimizer(1.3))
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
            sig = signature(sigx, krv, special)
            newrow = collect( 100.0 * sig[elm] for elm in newcols )
            quant[row, :] = newrow
            sbs = sortedbysig(newcols, NeXLUncertainties.value.(newrow))
            felm[row, :] = collect( sbs[i][1] for i in 1:size(fsig,2))
            fsig[row, :] = collect( sbs[i][2] for i in 1:size(fsig,2))
        end
    end
    mss = Vector{Missing}(missing,size(zep.data,1))
    fourelms = DataFrame( #
        FIRSTELM=felm[:,1], FIRSTPCT=fsig[:,1], #
        SECONDELM = length(felm)≥2 ? felm[:,2] : mss, SECONDPCT = length(felm)≥2 ? fsig[:,2] : mss, #
        THIRDELM = length(felm)≥3 ? felm[:,3] : mss,  THIRDPCT = length(felm)≥3 ?  fsig[:,3] : mss, #
        FOURTHELM = length(felm)≥4 ? felm[:,4] : mss, FOURTHPCT = length(felm)≥4 ? fsig[:,4] : mss) #
    quantRes = DataFrame((convert(Symbol, elm) => quant[:, i] for (i, elm) in enumerate(newcols))...)
    return hcat(fourelms, quantRes)
end


function buildNewZep(quant::DataFrame)
    # Remove old items...
    removeme = map(elm -> convert(Symbol, elm), zep.elms)
    append!(removeme, COMPOSITIONAL_COLUMNS)
    append!(removeme, CLASS_COLUMNS)
    remaining = copy(zep.data[:, filter(f -> !(f in removeme), names(zep.data))])
    return hcat(remaining, quant)
end
