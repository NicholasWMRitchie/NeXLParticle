using DataDeps

register(DataDep("ZepTestArtifact",
    """
    Dataset:    Zeppelin Test Data Set
    Author:     Nicholas W. M. Ritchie
    License:    CC-SA 3.0
    """,
    "https://drive.google.com/uc?export=download&id=1TZh4zbw2VY6QFlTfTpiPIG1U5uWepXRh",
    "07115ae49fe10b89f0782288d7a023c36909eb77ebba707658b7eef7ed440d22",
    post_fetch_method=DataDeps.unpack
))
ENV["DATADEPS_ALWAYS_ACCEPT"]="true"
# Get the path of the "Zeppelin Test Data Set", downloading it if necessary

include("zeppelin.jl")
include("translate.jl")
include("blob.jl")