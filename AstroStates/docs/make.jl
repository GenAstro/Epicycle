using AstroStates
using Documenter

DocMeta.setdocmeta!(AstroStates, :DocTestSetup, :(using AstroStates); recursive=true)

makedocs(;
    modules=[AstroStates],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroStates.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroStates/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroStates",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)