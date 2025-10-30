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

# Build documentation while everything is hot in memory
println("\n📚 Building documentation...")

# Add Documenter to current environment if needed
try
    using Documenter
    println("  ✅ Documenter already available")
catch
    println("  ➕ Installing Documenter...")
    Pkg.add("Documenter")
    using Documenter
end

# List of packages to build docs for
packages_to_document = [
    "AstroBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroCoords", "AstroModels", "AstroMan", "AstroFun", 
    "AstroProp", "AstroSolve"
]

println("🏗️  Building documentation for $(length(packages_to_document)) packages...")

for pkg_name in packages_to_document
    println("\n📖 Building docs for $pkg_name...")
    
    docs_make_path = joinpath(pkg_name, "docs", "make.jl")
    if !isfile(docs_make_path)
        println("  ⚠️  No docs/make.jl found for $pkg_name, skipping...")
        continue
    end
    
    try
        println("  🔨 Running $docs_make_path...")
        include(joinpath("..", docs_make_path))
        println("  ✅ Documentation built successfully for $pkg_name")
    catch e
        println("  ❌ Failed to build docs for $pkg_name: $e")
        exit(1)
    end
end

println("\n🎉 All documentation built successfully!")