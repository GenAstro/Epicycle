using AstroFun
using Documenter

DocMeta.setdocmeta!(AstroFun, :DocTestSetup, :(using AstroFun); recursive=true)

makedocs(;
    modules=[AstroFun],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroFun.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroFun/",
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
    dirname="AstroFun",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
