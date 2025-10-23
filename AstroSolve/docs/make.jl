using AstroSolve
using Documenter

DocMeta.setdocmeta!(AstroSolve, :DocTestSetup, :(using AstroSolve); recursive=true)

makedocs(;
    modules=[AstroSolve],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroSolve.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
