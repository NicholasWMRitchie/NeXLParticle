using NeXLMatrixCorrection

using Memoization  # Implements caching of k-ratio relative to pure....

@memoize function kpure( #
    std::Material,
    elm::Element,
    lines:: Vector{CharXRay},
    stdE0, unkE0,
    stdθ, unkθ,
    mc::Type{MatrixCorrection},
    fc::Type{FluorescenceCorrection}
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
    fc::Type{FluorescenceCorrection} = ReedFluorescence

)::Dict{Element,Float64}
    kzs = Dict{Element,AbstractFloat}()
    for kr in krs
        # k = (Iunk/Istd)/(Istd/Iz) = (Iunk/Iz) = kr.kratio/k(unk,std)
        kstd = isapprox(kr.standard[kr.element], 1.0, atol = 0.001) ? #
            kr.standard[kr.element] : #
            # Can't hand in full kr or memoization fails (since kr.kratio changes...)
            kpure(kr.standard, kr.element, kr.lines, #
                kr.stdProps[:BeamEnergy], kr.unkProps[:BeamEnergy], #
                kr.stdProps[:TakeOffAngle], kr.unkProps[:TakeOffAngle], #
                mc, fc)
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
