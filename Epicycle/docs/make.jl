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
    ),
    pages=[
        "Using Epicycle" => [
            "Welcome to Epicycle" => "welcome.md",
            "Getting Started" => [
                "Installation" => "installation.md",
                "Running Examples" => "running_examples.md",
                "Sample Missions" => "sample_missions.md",
                "Getting Help" => "getting_help.md",
            ],
            "Tour of Epicycle" => [
                "System Architecture" => "system_architecture.md",
                "Package Overview" => "package_overview.md",
            ],
        ],
        "Tutorials" => [
            "Unit Examples" => "unit_examples.md",
            "Complete Examples" => "complete_examples.md",
        ],
    ],
    warnonly=true,         # Just warn, don't error
    checkdocs=:none        # Skip docstring completeness checks
)
