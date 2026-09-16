# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroRoutines
using Documenter

DocMeta.setdocmeta!(AstroRoutines, :DocTestSetup, :(using AstroRoutines); recursive=true)

makedocs(;
    modules=[AstroRoutines],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroRoutines.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroRoutines/",
        edit_link="main",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=false,
    ),
    pages=[
        "index.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroRoutines",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
