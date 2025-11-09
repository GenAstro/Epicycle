#!/usr/bin/env julia
"""
CI Setup Environment Script

Sets up the Julia environment for CI by:
1. Activating the project
2. Instantiating dependencies
3. Developing all workspace packages

This script should be run first in the CI pipeline.
"""

println("🔧 Setting up CI environment...")

using Pkg

# Activate the main project
println("📦 Activating project environment...")
Pkg.activate(".")

# Instantiate to get all registered dependencies
println("📥 Installing dependencies...")
Pkg.instantiate()

# Develop all workspace packages
packages = [
    "AstroBase",
    "AstroStates", 
    "AstroEpochs",
    "AstroUniverse",
    "AstroFrames",
    "AstroModels",
    "AstroManeuvers",
    "AstroFun",
    "AstroProp",
    "AstroSolve",
    "Epicycle"
]

println("🔗 Developing workspace packages...")
for pkg in packages
    println("  → Developing $pkg...")
    Pkg.develop(path=pkg)
end

# Resolve any dependency conflicts
println("🎯 Resolving dependencies...")
Pkg.resolve()

println("✅ Environment setup complete!")
println("📋 Package status:")
Pkg.status()

