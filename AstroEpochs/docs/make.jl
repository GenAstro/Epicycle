using AstroEpochs
using Documenter

DocMeta.setdocmeta!(AstroEpochs, :DocTestSetup, :(using AstroEpochs); recursive=true)

makedocs(;
    modules=[AstroEpochs],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroEpochs.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroEpochs/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
        "Internal Details" => "internal.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)
