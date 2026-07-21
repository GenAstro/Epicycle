using Pkg
Pkg.add(PackageSpec(
    name = "AstroForceModels",
    url  = "https://github.com/HAMMERHEAD-Space/AstroForceModels.jl",
    rev  = "master",
))
Pkg.status("AstroForceModels")
