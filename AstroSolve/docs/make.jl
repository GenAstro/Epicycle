# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using AstroSolve
using Documenter

DocMeta.setdocmeta!(AstroSolve, :DocTestSetup, :(using AstroSolve); recursive=true)

makedocs(;
    modules=[AstroSolve],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroSolve.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroSolve/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "AstroSolve.jl" => "index.md",
        "Concepts"      => "concepts.md",
        "Parameter optimization" => "optimization.md",
        "Optimal control"        => "optimal_control.md",
        "Estimation"             => "estimation.md",
        "API reference"          => "api.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroSolve",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
