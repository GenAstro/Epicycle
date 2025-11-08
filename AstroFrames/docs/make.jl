using AstroFrames
using Documenter

DocMeta.setdocmeta!(AstroFrames, :DocTestSetup, :(using AstroFrames); recursive=true)

makedocs(;
    modules=[AstroFrames],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroFrames.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroFrames/",
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
    dirname="AstroFrames",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
