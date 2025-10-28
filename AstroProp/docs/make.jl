using AstroProp
using Documenter

DocMeta.setdocmeta!(AstroProp, :DocTestSetup, :(using AstroProp); recursive=true)

makedocs(;
    modules=[AstroProp],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroProp.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroProp/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Public API" => "api.md",
        "Developer API" => "internal.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroProp",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
