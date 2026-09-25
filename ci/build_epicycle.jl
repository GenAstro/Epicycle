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
# The package list is not written here. It is read from the workspace `projects` entry in the
# repo-root Project.toml, which is the one place that has to be right for Julia to resolve at
# all. Eight hand-maintained copies of this list is how AstroRoutines and EpicycleIO came to be
# missing from CI after they moved in, and how CairoMakie survived in the workspace after
# graphics left the umbrella.
using TOML
const WORKSPACE_PACKAGES = TOML.parsefile(
    joinpath(dirname(@__DIR__), "Project.toml"))["workspace"]["projects"]
Pkg.activate(".")

# Set coverage environment BEFORE loading any packages
ENV["JULIA_CODE_COVERAGE"] = "user"

# This will trigger compilation of Epicycle and all Astro packages (with coverage)
println("⚡ Loading Epicycle (this will trigger compilation with coverage)...")
@time using Epicycle

println("✅ Epicycle build complete!")
println("📊 Loaded packages:")

# Verify all packages are available
packages_to_check = Symbol.(filter(!=("Epicycle"), WORKSPACE_PACKAGES))

# Load each package rather than assuming `using Epicycle` pulled it in. AstroRoutines and
# EpicycleIO are workspace members that the umbrella deliberately does not depend on, so
# `isdefined(Main, pkg)` reports them as failures when they are fine.
#
# And record the failure rather than exiting. CI exists to report everything wrong in one run;
# exit(1) here stopped the whole pipeline before a single test ran, and skipped the summary.
load_failures = String[]
for pkg in packages_to_check
    try
        @eval using $pkg
        println("  ✅ $pkg loaded successfully")
    catch e
        println("  ❌ $pkg failed to load: ", sprint(showerror, e))
        push!(load_failures, String(pkg))
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
cdn_failed = false

# ---------------------------------------------------------------------------
# PHASE A0: Are EpicycleIO's third-party browser libraries still reachable?
# Plotly and Cesium are not vendored — dashboard.html pulls both from a CDN at
# exact pinned versions. Nothing else here touches a browser, so a withdrawn
# version would blank every user's plots with the suite still green. Cheap, and
# it runs first because it needs nothing built.
# ---------------------------------------------------------------------------
println("
📡 Checking EpicycleIO browser assets...")
let root = dirname(@__DIR__), script = joinpath(root, "EpicycleIO", "test", "check_cdn_assets.jl")
    try
        run(`$(Base.julia_cmd()) --project=$root $script`)
    catch
        global cdn_failed = true
        @error "EpicycleIO browser assets did not resolve — see the pinned versions in dashboard.html"
    end
end

# ---------------------------------------------------------------------------
# PHASE A: Run tests in a SUBPROCESS with --code-coverage=user.
# Subprocess exit flushes .cov files; the parent process (which runs
# generate_coverage.jl below) then reads them cleanly. Running tests
# in-process would leave coverage counts unflushed and generate_coverage.jl
# would see zeros (root cause of the 0% Codecov upload post 2026-08-07).
# Do NOT exit on failure yet — we still want to generate and upload
# whatever coverage was produced.
# ---------------------------------------------------------------------------
println("\n🧪 Running tests with coverage (subprocess)...")
try
    julia_exe   = Base.julia_cmd().exec[1]
    project_dir = dirname(@__DIR__)                       # repo root
    test_script = joinpath(project_dir, "Epicycle", "util", "test_all_packages.jl")
    test_cmd    = `$julia_exe --project=$project_dir --code-coverage=user $test_script`
    println("   → $test_cmd")
    run(test_cmd)
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
packages_to_document = WORKSPACE_PACKAGES

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
cdn_status = cdn_failed ? "❌ unreachable" : "✅ reachable"
println("Browser assets (Plotly, Cesium): $cdn_status")
tests_status = tests_failed ? "❌ FAILED" : "✅ passed"
docs_status = isempty(docs_failures) ? "✅ all built" : "❌ failed: " * join(docs_failures, ", ")

println("\n" * ("=" ^ 50))
println("BUILD SUMMARY")
println("=" ^ 50)
println("  Tests:  $tests_status")
println("  Docs:   $docs_status")

if tests_failed || !isempty(docs_failures) || cdn_failed || !isempty(load_failures)
    docs_list = join(docs_failures, ", ")
    load_list = join(load_failures, ", ")
    error("CI failed — tests_failed=$tests_failed, docs_failures=[$docs_list], " *
          "load_failures=[$load_list], browser_assets_unreachable=$cdn_failed")
end

println("\n🎉 Build, tests, coverage, and docs all completed successfully!")
