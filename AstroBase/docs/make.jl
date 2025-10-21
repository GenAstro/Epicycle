using AstroBase
using Documenter

DocMeta.setdocmeta!(AstroBase, :DocTestSetup, :(using AstroBase); recursive=true)

makedocs(;
    modules=[AstroBase],
    authors="Steve Hughes <steven.hughes@genastro.org>",
    sitename="AstroBase.jl",
    format=Documenter.HTML(;
        canonical="https://GenAstro.github.io/AstroBase.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/GenAstro/AstroBase.jl",
    devbranch="master",
)
