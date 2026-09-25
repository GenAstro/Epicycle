# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using AstroProp
using AstroModels        # SphericalDrag / SphericalSRP are defined here and re-exported by AstroProp
using Documenter

DocMeta.setdocmeta!(AstroProp, :DocTestSetup, :(using AstroProp); recursive=true)
DocMeta.setdocmeta!(AstroModels, :DocTestSetup, :(using AstroModels); recursive=true)

makedocs(;
    # AstroModels is included so @docs can resolve the re-exported spacecraft geometry types.
    # Its doctests run here too, so it gets its own DocTestSetup above; checkdocs=:none
    # suppresses "missing docstring" noise for the rest of AstroModels.
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
        "Force Models" => "force_models.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
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
