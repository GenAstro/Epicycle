using Epicycle
using Documenter

DocMeta.setdocmeta!(Epicycle, :DocTestSetup, :(using Epicycle); recursive=true)

makedocs(;
    modules=[Epicycle],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="Epicycle.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/",
        edit_link="main",
        assets=String[],
        sidebar_sitename=false,
        collapselevel=1,
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => [
            "Component Cheat Sheets" => "unit_examples.md",
            "End-to-End Examples" => "complete_examples.md",
        ],
        "Components" => "components.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="Epicycle",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
