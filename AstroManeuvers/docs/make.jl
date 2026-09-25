# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using AstroManeuvers
using Documenter

DocMeta.setdocmeta!(AstroManeuvers, :DocTestSetup, :(using AstroManeuvers); recursive=true)

makedocs(;
    modules=[AstroManeuvers],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroManeuvers.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroManeuvers/",
        edit_link="main",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=false,
    ),
    pages=[
        "index.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
    checkdocs=:none,       # Skip docstring completeness checks
    linkcheck=false        # Skip external link checking
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroManeuvers",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
