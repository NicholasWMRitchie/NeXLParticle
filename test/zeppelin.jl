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
    @test_broken length(classes(zep))==28
    @test maxspectrum(zep,1:10)[123] == 42.0
end
