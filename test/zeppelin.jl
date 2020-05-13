using NeXLParticle
using Test
using Pkg.Artifacts


zeptest = artifact"zeptest"

@testset "ZepTest" begin
    zep = Zeppelin(joinpath(zeptest,"test.hdz"))
    @test eachparticle(zep)==1:250
    hdz = header(zep)
    @test hdz["MAG0"]=="347 78 250 98.24 5.112"
    @test hdz["INSTRUMENT"]=="NIST's TESCAN MIRA3 in 217 F101"
    es = elms(zep)
    @test length(es)==32
    @test n"Sb" in es
    @test beamenergy(zep)==20.0e3
    @test probecurrent(zep)==0.961832
    @test magdata(zep,0)["Fields"]==78.0
    @test magdata(zep,0)["Area"]==5.112
    @test length(classes(zep))==52
    @test maxparticle(zep,1:10)[123] == 42.0
    relms = ( n"Ag", n"Al", n"Ba", n"Bi", n"Br", n"C", n"Ca", n"Ce", n"Cl", n"Co", n"Cr", n"Cu", n"F", n"Fe", n"K", n"Mg", n"Mn", n"Na", n"Nd", n"Ni",
        n"O", n"P", n"Pb", n"S", n"Sb", n"Si", n"Sn", n"Sr", n"Ti", n"V", n"W", n"Zn", n"Zr" )
    refs = Dict{Element,Spectrum}( elm => readEMSA(joinpath(zeptest,"Standards","$(elm.symbol) std.msa")) for elm in relms)
    det = matching(refs[n"Fe"],132.0)
    res = quantify(zep, det, refs, withUncertainty=true)
    @test isapprox(res[6,:FE], 94.87, atol=0.01)
    @test isapprox(res[26,:CA], 23.48, atol=0.01)
    @test isapprox(res[169, :BA], 39.37, atol=0.01)
    @test res[160,:FIRSTELM]==n"Si"
    @test res[172,:SECONDELM]==n"Ti"
    @test res[179,:THIRDELM]==n"Zn"
    @test ismissing(res[177,:FOURTHELM])
    @test res[195,:FOURTHELM]==n"K"
    NeXLParticle.writeZep(res, res.headerfile)
    res2 = Zeppelin(res.headerfile)
    @test isapprox(res2[6,:FE], 94.87, atol=0.01)
    @test isapprox(res2[26,:CA], 23.48, atol=0.01)
    @test isapprox(res2[169, :BA], 39.37, atol=0.01)
    @test res2[160,:FIRSTELM]==n"Si"
    @test res2[172,:SECONDELM]==n"Ti"
    @test res2[179,:THIRDELM]==n"Zn"
    @test ismissing(res2[177,:FOURTHELM])
    @test res2[195,:FOURTHELM]==n"K"
    #@time res=NeXLParticle.quantify(zep, det, refs, withUncertainty=false);
end
