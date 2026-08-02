#!/usr/bin/env julia
"""
CI Build Epicycle Script

Compiles the Epicycle package and all its dependencies.
This is where the heavy compilation work happens.

Order (changed): load -> TEST -> COVERAGE -> DOCS.
Rationale: tests + coverage now run *before* the docs loop, and a docs
failure no longer aborts the run before coverage is generated. Per-package
docs failures are collected and reported at the end instead of exiting mid-loop,
so one broken package (e.g. Epicycle, built last) can't strand the others or
suppress the coverage upload.
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

# ---------------------------------------------------------------------------
# Track failures across phases and fail at the very end. This guarantees that
# coverage is always generated/uploaded, and that a late docs failure does not
# hide test/coverage results.
# ---------------------------------------------------------------------------
tests_failed = false
docs_failures = String[]

# ---------------------------------------------------------------------------
# PHASE A: Run tests (with coverage). Do NOT exit on failure yet — we still
# want to generate and upload whatever coverage was produced.
# ---------------------------------------------------------------------------
println("\n🧪 Running tests with coverage...")
try
    # Path relative to project root, not ci directory
    test_script = joinpath("..", "Epicycle", "util", "test_all_packages.jl")
    include(test_script)
    println("✅ All tests completed successfully!")
catch e
    tests_failed = true
    println("❌ Tests failed: $e")
    println("   (continuing so coverage is still generated and uploaded)")
end

# ---------------------------------------------------------------------------
# PHASE B: Generate coverage. Runs regardless of test outcome so Codecov always
# receives whatever .cov data exists.
# ---------------------------------------------------------------------------
println("\n📈 Generating coverage...")
try
    include("generate_coverage.jl")
    println("✅ Coverage generation completed!")
catch e
    println("⚠️ Coverage generation failed: $e")
end

# ---------------------------------------------------------------------------
# PHASE C: Build documentation for all packages. Collect per-package failures
# instead of exiting mid-loop, so one broken package can't strand the others.
# ---------------------------------------------------------------------------
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
    catch e
        # Record the failure but keep going so later packages still deploy.
        push!(docs_failures, pkg_name)
        println("  ❌ Failed to build docs for $pkg_name: $e")
    end

    # Add delay to prevent GitHub Pages deployment conflicts
    if pkg_name != packages_to_document[end]  # Don't delay after the last package
        println("  ⏱️  Waiting 30 seconds before next deployment...")
        sleep(30)
    end
end

if isempty(docs_failures)
    println("\n🎉 All documentation built successfully!")
else
    println("\n⚠️  Documentation failed for: $(join(docs_failures, ", "))")
end

# ---------------------------------------------------------------------------
# Final status: fail the job if tests or any docs build failed, but only after
# coverage has been generated and (via the workflow's always() upload step)
# sent to Codecov.
# ---------------------------------------------------------------------------
tests_status = tests_failed ? "❌ FAILED" : "✅ passed"
docs_status = isempty(docs_failures) ? "✅ all built" : "❌ failed: " * join(docs_failures, ", ")

println("\n" * ("=" ^ 50))
println("BUILD SUMMARY")
println("=" ^ 50)
println("  Tests:  $tests_status")
println("  Docs:   $docs_status")

if tests_failed || !isempty(docs_failures)
    docs_list = join(docs_failures, ", ")
    error("CI failed — tests_failed=$tests_failed, docs_failures=[$docs_list]")
end

println("\n🎉 Build, tests, coverage, and docs all completed successfully!")
