using AstroProp
using Documenter

DocMeta.setdocmeta!(AstroProp, :DocTestSetup, :(using AstroProp); recursive=true)

makedocs(;
    modules=[AstroProp],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroProp.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
