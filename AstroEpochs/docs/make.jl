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
        collapselevel=1,
        sidebar_sitename=false,
    ),
    pages=[
        "index.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroEpochs",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)