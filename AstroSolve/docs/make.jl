using AstroSolve
using Documenter

DocMeta.setdocmeta!(AstroSolve, :DocTestSetup, :(using AstroSolve); recursive=true)

makedocs(;
    modules=[AstroSolve],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroSolve.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroSolve/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "AstroSolve.jl" => "index.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)
