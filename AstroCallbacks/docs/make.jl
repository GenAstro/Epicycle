using AstroCallbacks
using Documenter

DocMeta.setdocmeta!(AstroCallbacks, :DocTestSetup, :(using AstroCallbacks); recursive=true)

makedocs(;
    modules=[AstroCallbacks],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroCallbacks.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroCallbacks/",
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

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroCallbacks",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
