#!/usr/bin/env julia
"""
CI Build Epicycle Script

Compiles the Epicycle package and all its dependencies.
This is where the heavy compilation work happens.
"""

println("🏗️  Building Epicycle...")

using Pkg
Pkg.activate(".")

# This will trigger compilation of Epicycle and all Astro packages
println("⚡ Loading Epicycle (this will trigger compilation)...")
@time using Epicycle

println("✅ Epicycle build complete!")
println("📊 Loaded packages:")

# Verify all packages are available
packages_to_check = [
    :AstroBase, :AstroStates, :AstroEpochs, :AstroUniverse,
    :AstroCoords, :AstroModels, :AstroMan, :AstroFun, 
    :AstroProp, :AstroSolve
]

for pkg in packages_to_check
    if isdefined(Main, pkg)
        println("  ✅ $pkg loaded successfully")
    else
        println("  ❌ $pkg failed to load")
        exit(1)
    end
end

println("🎉 All packages loaded successfully!")