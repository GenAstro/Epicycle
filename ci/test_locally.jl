#!/usr/bin/env julia
"""
Local CI Testing Script

Simulates a fresh CI environment locally using Option 2 approach.
Creates an isolated Julia depot and runs the full CI sequence.

Usage: julia ci/test_locally.jl
"""

println("🧹 Setting up isolated CI test environment...")

# Create a temporary depot directory
temp_depot = joinpath(tempdir(), "ci_test_depot_$(rand(1000:9999))")
mkdir(temp_depot)

# Store original depot path
original_depot = get(ENV, "JULIA_DEPOT_PATH", "")

try
    # Set isolated depot
    ENV["JULIA_DEPOT_PATH"] = temp_depot
    println("📁 Using isolated depot: $temp_depot")
    
    # Get current directory
    project_dir = pwd()
    println("📂 Project directory: $project_dir")
    
    println("\n" * "="^50)
    println("🔧 PHASE 1: Setup Environment")
    println("="^50)
    
    # Run setup_environment.jl
    setup_cmd = `julia --project=$project_dir $project_dir/ci/setup_environment.jl`
    println("Running: $setup_cmd")
    run(setup_cmd)
    
    println("\n" * "="^50)
    println("🏗️  PHASE 2: Build Epicycle")
    println("="^50)
    
    # Run build_epicycle.jl
    build_cmd = `julia --project=$project_dir $project_dir/ci/build_epicycle.jl`
    println("Running: $build_cmd")
    run(build_cmd)
    
    println("\n" * "="^50)
    println("🧪 PHASE 3: Test Epicycle")
    println("="^50)
    
    # Run test_epicycle.jl
    test_cmd = `julia --project=$project_dir $project_dir/ci/test_epicycle.jl`
    println("Running: $test_cmd")
    run(test_cmd)
    
    println("\n" * "="^50)
    println("📈 PHASE 4: Generate Coverage")
    println("="^50)
    
    # Run generate_coverage.jl
    coverage_cmd = `julia --project=$project_dir $project_dir/ci/generate_coverage.jl`
    println("Running: $coverage_cmd")
    run(coverage_cmd)
    
    println("\n" * "="^50)
    println("🎉 LOCAL CI TEST COMPLETE!")
    println("="^50)
    println("✅ All phases completed successfully")
    println("🚀 Ready to push to GitHub!")
    
catch e
    println("\n" * "="^50)
    println("❌ LOCAL CI TEST FAILED!")
    println("="^50)
    println("Error: $e")
    println("\n🔍 Fix the issue and run again before pushing to GitHub")
    exit(1)
    
finally
    # Restore original depot
    if !isempty(original_depot)
        ENV["JULIA_DEPOT_PATH"] = original_depot
    else
        delete!(ENV, "JULIA_DEPOT_PATH")
    end
    
    # Clean up temp depot
    try
        rm(temp_depot, recursive=true, force=true)
        println("🧹 Cleaned up temporary depot")
    catch
        println("⚠️  Could not clean up temporary depot: $temp_depot")
    end
end