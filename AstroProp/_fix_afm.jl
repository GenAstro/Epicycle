using Pkg
# Revert the master pin and re-add from the registry
Pkg.rm(PackageSpec(name="AstroForceModels"))
Pkg.add("AstroForceModels")
Pkg.update()
Pkg.status()
