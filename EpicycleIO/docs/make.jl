# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using EpicycleIO
using Documenter

DocMeta.setdocmeta!(EpicycleIO, :DocTestSetup, :(using EpicycleIO); recursive=true)

makedocs(;
    modules=[EpicycleIO],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="EpicycleIO.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/Epicycle/EpicycleIO/",
        edit_link="main",
        assets=String[],
        collapselevel=1,
        sidebar_sitename=false,
    ),
    pages=[
        "Home"          => "index.md",
        "Data plotting" => "data_plotting.md",
        "3D plotting"   => "plotting_3d.md",
        "Dashboard"     => "dashboard.md",
        "Reporting"     => "reporting.md",
        "API"           => "api.md",
    ],
    warnonly=[:missing_docs, :cross_references],   # doctests and examples are fatal
)
