using NeXLParticle
using Test
using Pkg
using Pkg.Artifacts


# This is the path to the Artifacts.toml we will manipulate
artifacts_toml = joinpath(@__DIR__, "Artifacts.toml")

# Query the `Artifacts.toml` file for the hash bound to the name "iris"
# (returns `nothing` if no such binding exists)
zep_hash = artifact_hash("zeptest", artifacts_toml)

# If the name was not bound, or the hash it was bound to does not exist, create it!
if zep_hash == nothing || !artifact_exists(zep_hash)
    # create_artifact() returns the content-hash of the artifact directory once we're finished creating it
    zep_hash = create_artifact() do artifact_dir
        println("Downloading test data to $artifact_dir")
        tarball = joinpath(artifact_dir, "zeptest.tar.gz")
        download("https://drive.google.com/uc?export=download&id=1TZh4zbw2VY6QFlTfTpiPIG1U5uWepXRh",tarball)
        Pkg.probe_platform_engines!()
        Pkg.unpack(tarball, artifact_dir, verbose=true)
        rm(tarball)
    end
    # Now bind that hash within our `Artifacts.toml`.  `force = true` means that if it already exists,
    # just overwrite with the new content-hash.  Unless the source files change, we do not expect
    # the content hash to change, so this should not cause unnecessary version control churn.
    bind_artifact!(artifacts_toml, "zeptest", zep_hash)
end

# Get the path of the iris dataset, either newly created or previously generated.
# this should be something like `~/.julia/artifacts/dbd04e28be047a54fbe9bf67e934be5b5e0d357a`
zeptest = artifact_path(zep_hash)

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
    refs = Dict{Element,Spectrum}( elm => loadspectrum(joinpath(zeptest,"Standards","$(elm.symbol) std.msa")) for elm in relms)
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
