using AstroCoords
using Documenter

DocMeta.setdocmeta!(AstroCoords, :DocTestSetup, :(using AstroCoords); recursive=true)

makedocs(;
    modules=[AstroCoords],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroCoords.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
