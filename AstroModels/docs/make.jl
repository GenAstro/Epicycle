using AstroModels
using Documenter

DocMeta.setdocmeta!(AstroModels, :DocTestSetup, :(using AstroModels); recursive=true)

makedocs(;
    modules=[AstroModels],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroModels.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
