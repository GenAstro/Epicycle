using AstroCoords
using Documenter

DocMeta.setdocmeta!(AstroCoords, :DocTestSetup, :(using AstroCoords); recursive=true)

makedocs(;
    modules=[AstroCoords],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroCoords.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroCoords/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference Guide" => "api.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroCoords",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
