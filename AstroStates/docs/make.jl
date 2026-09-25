# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroStates
using Documenter

DocMeta.setdocmeta!(AstroStates, :DocTestSetup, :(using AstroStates); recursive=true)

makedocs(;
    modules=[AstroStates],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroStates.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroStates/",
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
    dirname="AstroStates",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
