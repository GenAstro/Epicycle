using Epicycle
using Documenter

DocMeta.setdocmeta!(Epicycle, :DocTestSetup, :(using Epicycle); recursive=true)

makedocs(;
    modules=[Epicycle],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="Epicycle.jl",
    format=Documenter.HTML(;
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
