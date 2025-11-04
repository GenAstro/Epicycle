# Test all packages in the Epicycle monorepo using Pkg.test()
# This properly handles test dependencies from [extras] sections

using Pkg

# Ensure coverage is enabled for all subprocesses
ENV["JULIA_CODE_COVERAGE"] = "user"

# Get the repo root directory (two levels up from this script)
script_dir = dirname(@__FILE__)
repo_root = dirname(dirname(script_dir))

# Store original project to restore later
original_project = Base.active_project()

# List of all packages in the monorepo
packages = [
    "AstroBase",
    "AstroStates",
    "AstroEpochs",
    "AstroUniverse",
    "AstroCoords", 
    "AstroMan",
    "AstroModels",
    "AstroFun",
    "AstroProp",
    "AstroSolve",
]

#packages = ["AstroStates"]

println("Testing all Epicycle packages with proper dependency handling...")
println("Repo root: $repo_root")
println("=" ^ 50)

# PHASE 1: Precompile all test dependencies to avoid noise during testing
println("🔧 PHASE 1: Precompiling test dependencies...")
println("-" ^ 50)

for pkg in packages
    pkg_path = joinpath(repo_root, pkg)
    if isdir(pkg_path)
        println("  → Precompiling $pkg test dependencies...")
        try
            Pkg.activate(pkg_path)
            Pkg.resolve()      # Resolve test environment
            Pkg.instantiate()  # Download and precompile test dependencies
        catch e
            println("    ⚠️  Warning: Could not precompile $pkg dependencies: $e")
        end
    end
end

# Restore original project after precompilation
if original_project !== nothing
    Pkg.activate(original_project)
end

println("✅ Precompilation complete!")
println("\n" * "=" ^ 50)
println("🧪 PHASE 2: Running tests (should be clean now)...")
println("=" ^ 50)

failed_packages = String[]

for pkg in packages
    println("\n🧪 Testing $pkg...")
    pkg_path = joinpath(repo_root, pkg)
    
    if !isdir(pkg_path)
        println("⚠️  Package directory not found: $pkg_path")
        continue
    end
    
    try
        # Activate package environment (deps already precompiled)
        Pkg.activate(pkg_path)
        
        # Run tests with Pkg.test() - should be faster now since coverage already enabled
        println("  → Running Pkg.test() (coverage already enabled globally)...")
        
        # Activate package environment 
        Pkg.activate(pkg_path)
        
        # Run tests - should not recompile since coverage already set
        Pkg.test()
        println("✅ $pkg tests passed")
        
    catch e
        println("❌ $pkg tests failed: $e")
        push!(failed_packages, pkg)
    finally
        # Always restore original project
        if original_project !== nothing
            Pkg.activate(original_project)
        end
    end
end

println("\n" * "=" ^ 50)
println("TEST SUMMARY:")
println("=" ^ 50)

# Ensure we're back in the original project
if original_project !== nothing
    Pkg.activate(original_project)
end

if isempty(failed_packages)
    println("🎉 All packages passed their tests!")
else
    println("❌ Failed packages: $(join(failed_packages, ", "))")
    println("📊 Passed: $(length(packages) - length(failed_packages))/$(length(packages))")
    error("Some packages failed their tests")
end