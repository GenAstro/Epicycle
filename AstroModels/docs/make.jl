using AstroModels
using Documenter

DocMeta.setdocmeta!(AstroModels, :DocTestSetup, :(using AstroModels); recursive=true)

makedocs(;
    modules=[AstroModels],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroModels.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroModels/",
        edit_link="main",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=false,
    ),
    pages=[
        "index.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none,       # Skip docstring completeness checks
    doctest=false,         # Skip all doctests completely
    linkcheck=false        # Skip external link checking
)
