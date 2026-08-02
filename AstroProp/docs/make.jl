using AstroProp
using AstroModels        # SphericalDrag / SphericalSRP are defined here and re-exported by AstroProp
using Documenter

DocMeta.setdocmeta!(AstroProp, :DocTestSetup, :(using AstroProp); recursive=true)

makedocs(;
    # AstroModels is included so @docs can resolve the re-exported spacecraft geometry types.
    # It has no jldoctest blocks, so no foreign doctests run; checkdocs=:none suppresses
    # "missing docstring" noise for the rest of AstroModels.
    modules=[AstroProp, AstroModels],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroProp.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroProp/",
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
    dirname="AstroProp",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
