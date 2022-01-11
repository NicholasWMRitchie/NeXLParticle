using NeXLParticle
using Test
using DataDeps



@testset "ZepTest" begin
    zeptest = datadep"ZepTestArtifact"
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
    @test length(classes(zep))==27
    @test maxparticle(zep,1:10)[123] == 42.0
    refs = references( [
        reference(n"Ag", joinpath(zeptest,"Standards","Ag std.msa")),
        reference(n"Al", joinpath(zeptest,"Standards","Al std.msa")),
        reference(n"Ba", joinpath(zeptest,"Standards","Ba std.msa")),
        reference(n"Bi", joinpath(zeptest,"Standards","Bi std.msa")),
        reference(n"Br", joinpath(zeptest,"Standards","Br std.msa")),
        reference(n"C", joinpath(zeptest,"Standards","C std.msa")),
        reference(n"Ca", joinpath(zeptest,"Standards","Ca std.msa")),
        reference(n"Ce", joinpath(zeptest,"Standards","Ce std.msa")),
        reference(n"Cl", joinpath(zeptest,"Standards","Cl std.msa")),
        reference(n"Co", joinpath(zeptest,"Standards","Co std.msa")),
        reference(n"Cr", joinpath(zeptest,"Standards","Cr std.msa")),
        reference(n"Cu", joinpath(zeptest,"Standards","Cu std.msa")),
        reference(n"F", joinpath(zeptest,"Standards","F std.msa")),
        reference(n"Fe", joinpath(zeptest,"Standards","Fe std.msa")),
        reference(n"K", joinpath(zeptest,"Standards","K std.msa")),
        reference(n"Mg", joinpath(zeptest,"Standards","Mg std.msa")),
        reference(n"Mn", joinpath(zeptest,"Standards","Mn std.msa")),
        reference(n"Na", joinpath(zeptest,"Standards","Na std.msa")),
        reference(n"Nd", joinpath(zeptest,"Standards","Nd std.msa")),
        reference(n"Ni", joinpath(zeptest,"Standards","Ni std.msa")),
        reference(n"O", joinpath(zeptest,"Standards","O std.msa")),
        reference(n"P", joinpath(zeptest,"Standards","P std.msa")),
        reference(n"Pb", joinpath(zeptest,"Standards","Pb std.msa")),
        reference(n"S", joinpath(zeptest,"Standards","S std.msa")),
        reference(n"Sb", joinpath(zeptest,"Standards","Sb std.msa")),
        reference(n"Si", joinpath(zeptest,"Standards","Si std.msa")),
        reference(n"Sn", joinpath(zeptest,"Standards","Sn std.msa")),
        reference(n"Sr", joinpath(zeptest,"Standards","Sr std.msa")),
        reference(n"Ti", joinpath(zeptest,"Standards","Ti std.msa")),
        reference(n"V", joinpath(zeptest,"Standards","V std.msa")),
        reference(n"W", joinpath(zeptest,"Standards","W std.msa")),
        reference(n"Zn", joinpath(zeptest,"Standards","Zn std.msa")),
        reference(n"Zr", joinpath(zeptest,"Standards","Zr std.msa")) ],
        132.0
    )
    @time quantify(zep, refs, withUncertainty=true)
    res = @time quantify(zep, refs, withUncertainty=false)
    @test isapprox(res[6,:FE], 94.87, atol=0.01)
    @test isapprox(res[26,:CA], 23.48, atol=0.01)
    @test isapprox(res[169, :BA], 39.38, atol=0.01)
    @test res[160,:FIRSTELM]==n"Si"
    @test res[172,:SECONDELM]==n"Ti"
    @test res[179,:THIRDELM]==n"Zn"
    @test ismissing(res[177,:FOURTHELM])
    @test res[195,:FOURTHELM]==n"K"
    @test classes(zep)==classes(res)
    outfile = tempname()*".hdz"
    NeXLParticle.writeZep(res, outfile)
    res2 = Zeppelin(outfile)
    @test classes(zep) == classes(res2)
    @test isapprox(res2[6,:FE], 94.87, atol=0.01)
    @test isapprox(res2[26,:CA], 23.48, atol=0.01)
    @test isapprox(res2[169, :BA], 39.38, atol=0.01)
    @test res2[160,:FIRSTELM]==n"Si"
    @test res2[172,:SECONDELM]==n"Ti"
    @test res2[179,:THIRDELM]==n"Zn"
    @test ismissing(res2[177,:FOURTHELM])
    @test res2[195,:FOURTHELM]==n"K"
    @test all(res[:,"CLASS"] .== zep[:,"CLASS"])
    @test all(res[:,"CLASS"] .== res2[:,"CLASS"])
    res3 = @time classify(res, NeXLParticle.BaseRules)
    writeZep(res3, joinpath(homedir(),"Desktop","tmp.hdz"))
    @test repr(res3[2, "CLASS"]) == "Pyrite (FeS2)"
    @test (r->(r.FE > 20) && (r.S > 20) && (r.FE + r.S > 80))(res3[2,:]) 
    @test repr(res3[20, "CLASS"]) == "Dolomite"
    @test (r -> (r.CA > 40) && (r.MG > 10) && (r.CA + r.MG > 70))(res3[20,:]) 
end