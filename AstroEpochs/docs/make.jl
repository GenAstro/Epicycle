using AstroEpochs
using Documenter

DocMeta.setdocmeta!(AstroEpochs, :DocTestSetup, :(using AstroEpochs); recursive=true)

makedocs(;
    modules=[AstroEpochs],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroEpochs.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
