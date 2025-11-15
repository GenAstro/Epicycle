using AstroBase
using Documenter

DocMeta.setdocmeta!(AstroBase, :DocTestSetup, :(using AstroBase); recursive=true)

makedocs(;
    modules=[AstroBase],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroBase.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroBase/",
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
    dirname="AstroBase",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)