using AstroMan
using Documenter

DocMeta.setdocmeta!(AstroMan, :DocTestSetup, :(using AstroMan); recursive=true)

makedocs(;
    modules=[AstroMan],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroMan.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroMan/",
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
    dirname="AstroMan",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
