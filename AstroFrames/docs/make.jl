# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using AstroFrames
using Documenter

DocMeta.setdocmeta!(AstroFrames, :DocTestSetup, :(using AstroFrames); recursive=true)

makedocs(;
    modules=[AstroFrames],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroFrames.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/AstroFrames/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference Guide" => "api.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
    checkdocs=:none        # Skip docstring completeness checks
)

deploydocs(;
    repo="github.com/GenAstro/Epicycle",
    target="build",
    dirname="AstroFrames",
    devbranch="main",
    push_preview=true,
    deploy_config=Documenter.GitHubActions()  # Uses GITHUB_TOKEN
)
