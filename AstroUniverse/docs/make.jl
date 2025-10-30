using AstroUniverse
using Documenter

DocMeta.setdocmeta!(AstroUniverse, :DocTestSetup, :(using AstroUniverse); recursive=true)

makedocs(;
    modules=[AstroUniverse],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroUniverse.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroUniverse/",
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
    dirname="AstroUniverse",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
