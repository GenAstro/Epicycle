#!/usr/bin/env julia
"""
CI Test Epicycle Script

Runs comprehensive tests for all packages in the workspace.
Assumes packages are already compiled from build step.
"""

println("🧪 Testing Epicycle...")

using Pkg
Pkg.activate(".")

# Load Epicycle (should be fast since already compiled)
println("📦 Loading Epicycle...")
using Epicycle

# Run the comprehensive test suite
println("🏃 Running comprehensive tests...")

# Get the project root directory (parent of ci directory)
project_root = dirname(@__DIR__)
test_script = joinpath(project_root, "Epicycle", "util", "test_all_packages.jl")
if isfile(test_script)
    println("  → Using test script: $test_script")
    include(test_script)
else
    println("  → Test script not found, running individual package tests...")
    
    # Fallback: run tests for each package individually
    packages = [
        "AstroBase", "AstroStates", "AstroEpochs", "AstroUniverse",
        "AstroCoords", "AstroModels", "AstroMan", "AstroFun", 
        "AstroProp", "AstroSolve", "Epicycle"
    ]
    
    # Get the project root directory
    project_root = dirname(@__DIR__)
    
    for pkg in packages
        pkg_path = joinpath(project_root, pkg)
        if isdir(pkg_path)
            println("  → Testing $pkg...")
            try
                # Set environment variables for coverage
                ENV["JULIA_NUM_THREADS"] = "auto"
                
                # First activate the package to ensure test dependencies are available
                println("    → Activating $pkg environment at $pkg_path")
                Pkg.activate(pkg_path)
                Pkg.instantiate()  # Install test dependencies
                
                # Run tests with coverage enabled - stay in root directory
                # but use the package's test environment
                cd(project_root) do
                    # Enable coverage and run tests
                    withenv("JULIA_CODE_COVERAGE" => "user") do
                        Pkg.test(pkg; coverage=true)
                    end
                end
                
                println("    ✅ $pkg tests passed")
                
                # Return to main project
                Pkg.activate(project_root)
            catch e
                println("    ❌ $pkg tests failed: $e")
                # Return to main project even on failure
                Pkg.activate(project_root)
                exit(1)
            end
        end
    end
end

println("🎉 All tests completed successfully!")

# Generate coverage immediately while .cov files exist
println("\n📈 Generating coverage immediately...")
try
    include(joinpath(dirname(@__DIR__), "ci", "generate_coverage.jl"))
    println("✅ Coverage generation completed!")
catch e
    println("⚠️  Coverage generation failed: $e")
end