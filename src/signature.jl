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
    mc::Type{MatrixCorrection},
    fc::Type{FluorescenceCorrection},
)::AbstractFloat
    unkzaf = ZAF(mc, fc, std, lines, stdE0)
    stdzaf = ZAF(mc, fc, pure(elm), lines, unkE0)
    return k(unkzaf, stdzaf, stdθ, unkθ)
end

"""
    signature( #
        krs::Array{KRatio},
        strip::Array{Element} = [n"O", n"C"],
        includeStrip::Array{Element} = [n"O"];
        mc::Type{MatrixCorrection} = XPPCorrection,
        fc::Type{FluorescenceCorrection} = ReedFluorescence

    )::Dict{Element,Float64}

Computes a "particle signature" from the specified set of k-ratios.  A particle signature
is similar to normalized k-ratios relative to pure elements except 1) certain elements
can be removed from the normalization; and 2) these elements can be include in the result as
un-normalized.
"""
function signature( #
    krs::Array{KRatio},
    strip::Array{Element} = [n"O", n"C"],
    includeStrip::Array{Element} = [n"O"];
    mc::Type{MatrixCorrection} = XPPCorrection,
    fc::Type{FluorescenceCorrection} = ReedFluorescence,

)::Dict{Element,Float64}
    kzs = Dict{Element,AbstractFloat}()
    for kr in krs
        # k = (Iunk/Istd)/(Istd/Iz) = (Iunk/Iz) = kr.kratio/k(unk,std)
        kstd = isapprox(kr.standard[kr.element], 1.0, atol = 0.001) ? #
               kr.standard[kr.element] : #
            # Can't hand in full kr or memoization fails (since kr.kratio changes...)
               kpure(
            kr.standard,
            kr.element,
            kr.lines, #
            kr.stdProps[:BeamEnergy],
            kr.unkProps[:BeamEnergy], #
            kr.stdProps[:TakeOffAngle],
            kr.unkProps[:TakeOffAngle], #
            mc,
            fc,
        )
        kzs[kr.element] = kr.kratio / kstd
    end
    res = Dict{Element,AbstractFloat}()
    notstripped = filter(z -> !(z in strip), keys(krs))
    norm = sum(kzs[elm] for elm in notstripped)
    res = Dict(elm => kzs[elm] / norm for elm in notstripped)
    for elm in intersect(strip, includeStrip)
        res[elm] = kzs[elm]
    end
end

function process(
    zep::Zeppelin,
    det::Detector,
    refs::Dict{Element,Spectrum};
    strip::Array{Element} = [n"O", n"C"],
    includeStrip::Array{Element} = [n"O"],
)
    # First, strip out all old elemental data items
    removeme = map(elm -> convert(Symbol, elm), zep.elms)
    append!(removeme, COMPOSITIONAL_COLUMNS)
    append!(removeme, CLASS_COLUMNS)
    res = copy(zep.data[:, filter(f -> !(f in removeme), names(zep.data))])
    filt = buildfilter(NeXLSpectrum.GaussianFilter, det)
    refs = Reference[]
    for elm in keys(refs)
        lines = ktransitions
        if elm > n"Ca"
            lines = union(ktransitions, ltransitions)
        elseif elm > n"Ba"
            lines = union(ltransitions, mtransitions)
        end
        append!(
            refs,
            map(r -> filter(refs[elm], det, r, filt, 1.0 / dose(refs[elm])), NeXLSpectrum.charFeature(elm, lines)),
        )
    end
    newcols = map(z -> convert(Symbol, z), filter(z -> !((z in strip) && !(z in includeStrip)), keys(refs)))
    quant = DataFrame((sym => [] for sym in newCols))
    for row in eachparticle(zep)
        unk = spectrum(zep, row, withImgs = false)
        if !ismissing(unk)
            res = fit(unk, filt, refs, false)
            sig = signature(res.kratios, strip, includeStrip)
            row = [sig[col] for col in newcols]
            append(!quant, row)
        else
            row = [missing for col in newcols]
        end
        push!(quant, row)
    end
    return hcat(res, quant)
end
