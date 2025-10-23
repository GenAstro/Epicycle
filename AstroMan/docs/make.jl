using AstroMan
using Documenter

DocMeta.setdocmeta!(AstroMan, :DocTestSetup, :(using AstroMan); recursive=true)

makedocs(;
    modules=[AstroMan],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroMan.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
