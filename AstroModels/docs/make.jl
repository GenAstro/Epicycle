using AstroModels
using Documenter

DocMeta.setdocmeta!(AstroModels, :DocTestSetup, :(using AstroModels); recursive=true)

makedocs(;
    modules=[AstroModels],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroModels.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroModels/",
        edit_link="main",
        assets=String[],
        sidebar_sitename=false,
    ),
    pages=[
        "QuickStart" => "index.md",
        "Spacecraft" => [
            "Overview" => "spacecraft.md",
            "Time" => "spacecraft_time.md",
            "State" => "spacecraft_state.md",
            "Coordinate System" => "spacecraft_coord_sys.md",
            "Mass" => "spacecraft_mass.md",
            "CAD Model" => "spacecraft_cad_model.md",
            "History" => "history.md",
        ],
        "Reference" => "reference.md",
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none,       # Skip docstring completeness checks
    doctest=false,         # Skip all doctests completely
    linkcheck=false        # Skip external link checking
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroModels",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
