#!/usr/bin/env julia
"""
Development Test Script for test_epicycle.jl

Runs the test phase in your current development environment
without the overhead of setting up a pristine environment.

Usage: julia ci/test_epicycle_dev.jl
"""

println("🧪 Testing Epicycle (Development Environment)...")

# Activate the project (should already be active in dev env)
using Pkg
Pkg.activate(".")

println("📦 Loading Epicycle...")
using Epicycle

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
        "EpicycleBase", "AstroStates", "AstroEpochs", "AstroUniverse",
        "AstroFrames", "AstroModels", "AstroManeuvers", "AstroCallbacks", 
        "AstroProp", "AstroSolve", "Epicycle"
    ]
    
    for pkg in packages
        println("  🧪 Testing $pkg...")
        try
            # First activate the package to ensure test dependencies are available
            pkg_path = joinpath(project_root, pkg)
            if isdir(pkg_path)
                println("    → Activating $pkg environment at $pkg_path")
                Pkg.activate(pkg_path)
                Pkg.instantiate()  # Install test dependencies
                Pkg.test()  # Run tests in the package's own environment
                println("  ✅ $pkg tests passed")
                
                # Return to main project
                Pkg.activate(project_root)
            else
                println("    → Package directory not found: $pkg_path")
                Pkg.test(pkg)  # Fallback to global test
                println("  ✅ $pkg tests passed")
            end
        catch e
            println("  ❌ $pkg tests failed: $e")
            # Return to main project even on failure
            Pkg.activate(project_root)
            exit(1)
        end
    end
end

println("🎉 All tests completed successfully!")