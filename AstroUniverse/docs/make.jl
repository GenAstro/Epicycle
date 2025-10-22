using AstroUniverse
using Documenter

DocMeta.setdocmeta!(AstroUniverse, :DocTestSetup, :(using AstroUniverse); recursive=true)

makedocs(;
    modules=[AstroUniverse],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroUniverse.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
