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

# Set coverage environment globally for all tests
ENV["JULIA_CODE_COVERAGE"] = "user"
println("🐛 DEBUG: Set JULIA_CODE_COVERAGE globally to: $(get(ENV, "JULIA_CODE_COVERAGE", "NOT SET"))")

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
        "AstroRoutines", "EpicycleBase", "AstroStates", "AstroEpochs", "AstroUniverse",
        "AstroFrames", "AstroModels", "AstroManeuvers", "AstroCallbacks", 
        "AstroProp", "AstroSolve", "Epicycle"
    ]
    
    # Get the project root directory
    project_root = dirname(@__DIR__)
    
    for pkg in packages
        pkg_path = joinpath(project_root, pkg)
        if isdir(pkg_path)
            println("  → Testing $pkg...")
            try
                # Set environment variables for coverage GLOBALLY
                ENV["JULIA_NUM_THREADS"] = "auto"
                ENV["JULIA_CODE_COVERAGE"] = "user"  # Set globally for subprocesses
                
                # DEBUG: Show environment and directory before testing
                println("🐛 DEBUG: About to test $pkg")
                println("🐛 DEBUG: Current working directory: $(pwd())")
                println("🐛 DEBUG: Package path: $pkg_path")
                println("🐛 DEBUG: JULIA_CODE_COVERAGE set to: $(get(ENV, "JULIA_CODE_COVERAGE", "NOT SET"))")
                
                # First activate the package to ensure test dependencies are available
                println("    → Activating $pkg environment at $pkg_path")
                Pkg.activate(pkg_path)
                Pkg.instantiate()  # Install test dependencies
                
                # Run tests with coverage enabled - stay in root directory
                # but use the package's test environment
                cd(project_root) do
                    # Coverage is now enabled globally via ENV
                    println("🐛 DEBUG: About to run tests with coverage for $pkg")
                    println("🐛 DEBUG: Working directory during test: $(pwd())")
                    println("🐛 DEBUG: JULIA_CODE_COVERAGE during test: $(get(ENV, "JULIA_CODE_COVERAGE", "NOT SET"))")
                    
                    # Ensure environment variable is set in this context too
                    withenv("JULIA_CODE_COVERAGE" => "user") do
                        Pkg.test(pkg; coverage=true)
                    end
                    
                    # DEBUG: Check for .cov files immediately after test
                    println("🐛 DEBUG: Checking for .cov files immediately after testing $pkg...")
                    pkg_src = joinpath(project_root, pkg, "src")
                    if isdir(pkg_src)
                        println("🐛 DEBUG: Contents of $pkg_src after test:")
                        for item in readdir(pkg_src)
                            println("    📄 $item")
                            if endswith(item, ".cov")
                                println("      🎯 Found .cov file: $item")
                            end
                        end
                    else
                        println("🐛 DEBUG: $pkg_src does not exist")
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

# DEBUG: Final check before coverage generation
println("🐛 DEBUG: Final state before coverage generation:")
println("🐛 DEBUG: Current working directory: $(pwd())")
println("🐛 DEBUG: JULIA_CODE_COVERAGE: $(get(ENV, "JULIA_CODE_COVERAGE", "NOT SET"))")

# Quick scan for .cov files before coverage generation
println("🐛 DEBUG: Quick scan for .cov files in all package src directories:")
packages_to_check = [
    "AstroRoutines", "EpicycleBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroFrames", "AstroModels", "AstroManeuvers", "AstroCallbacks", 
    "AstroProp", "AstroSolve", "Epicycle"
]

for pkg in packages_to_check
    src_dir = "$pkg/src"
    if isdir(src_dir)
        cov_count = length(filter(f -> endswith(f, ".cov"), readdir(src_dir)))
        println("  📊 $pkg/src: $cov_count .cov files")
    else
        println("  ❌ $pkg/src: directory not found")
    end
end

try
    include(joinpath(dirname(@__DIR__), "ci", "generate_coverage.jl"))
    println("✅ Coverage generation completed!")
catch e
    println("⚠️  Coverage generation failed: $e")
end