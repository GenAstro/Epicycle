using AstroFun
using Documenter

DocMeta.setdocmeta!(AstroFun, :DocTestSetup, :(using AstroFun); recursive=true)

makedocs(;
    modules=[AstroFun],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroFun.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
