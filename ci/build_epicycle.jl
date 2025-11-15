#!/usr/bin/env julia
"""
CI Build Epicycle Script

Compiles the Epicycle package and all its dependencies.
This is where the heavy compilation work happens.
"""

println("🏗️  Building Epicycle...")

using Pkg
Pkg.activate(".")

# Set coverage environment BEFORE loading any packages
ENV["JULIA_CODE_COVERAGE"] = "user"

# This will trigger compilation of Epicycle and all Astro packages (with coverage)
println("⚡ Loading Epicycle (this will trigger compilation with coverage)...")
@time using Epicycle

println("✅ Epicycle build complete!")
println("📊 Loaded packages:")

# Verify all packages are available
packages_to_check = [
    :EpicycleBase, :AstroStates, :AstroEpochs, :AstroUniverse,
    :AstroFrames, :AstroModels, :AstroManeuvers, :AstroCallbacks,
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
    "EpicycleBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroFrames", "AstroModels", "AstroManeuvers", "AstroCallbacks", 
    "AstroProp", "AstroSolve", "Epicycle"
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
        
        # Add delay to prevent GitHub Pages deployment conflicts
        if pkg_name != packages_to_document[end]  # Don't delay after the last package
            println("  ⏱️  Waiting 30 seconds before next deployment...")
            sleep(30)
        end
    catch e
        println("  ❌ Failed to build docs for $pkg_name: $e")
        exit(1)
    end
end

println("\n🎉 All documentation built successfully!")

# Run tests while everything is hot in memory
println("\n🧪 Running tests with coverage...")

try
    # Path relative to project root, not ci directory
    test_script = joinpath("..", "Epicycle", "util", "test_all_packages.jl")
    include(test_script)
    println("✅ All tests completed successfully!")
catch e
    println("❌ Tests failed: $e")
    exit(1)
end

# Generate coverage
println("\n📈 Generating coverage...")
try
    include("generate_coverage.jl")
    println("✅ Coverage generation completed!")
catch e
    println("⚠️ Coverage generation failed: $e")
end

println("\n🎉 Tests run and coverage generated!")