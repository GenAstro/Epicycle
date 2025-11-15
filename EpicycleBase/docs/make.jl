using EpicycleBase
using Documenter

DocMeta.setdocmeta!(EpicycleBase, :DocTestSetup, :(using EpicycleBase); recursive=true)

makedocs(;
    modules=[EpicycleBase],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="EpicycleBase.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/EpicycleBase/",
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
    dirname="EpicycleBase",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)