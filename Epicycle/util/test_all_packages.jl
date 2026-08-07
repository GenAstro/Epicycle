# Test all packages in the Epicycle monorepo from a single root Julia session.
#
# Design (see AstroSolve-coverage investigation, 2026-08-06):
#   - Pkg.test spawns a per-package sandbox subprocess whose Manifest can't see
#     dev-only umbrella packages like Epicycle. That masked ~half of AstroSolve's
#     test coverage on CI (three test files that do `using Epicycle` were failing
#     to load).
#   - TestEnv.activate(pkg) layers the package's [extras] onto the current
#     root env instead of building a fresh sandbox, so dev'd deps resolve.
#     No sandbox subprocess, no re-resolve, no re-precompile per package.
#
# Requirements:
#   - Must be called from the root project (workspace env with all sub-packages dev'd).
#   - Julia must be launched with `--code-coverage=user` (or equivalent) for .cov
#     files to be generated; the ENV var alone doesn't turn coverage on mid-process.

using Pkg
using TestEnv

script_dir = dirname(@__FILE__)
repo_root  = dirname(dirname(script_dir))

# Confirm we're running from the root workspace project (has all packages dev'd).
root_project = joinpath(repo_root, "Project.toml")
Base.active_project() == root_project ||
    @warn "test_all_packages.jl: expected root project $root_project active; got $(Base.active_project()). TestEnv layering may fail to resolve dev-only deps."

packages = [
    "EpicycleBase",
    "AstroStates",
    "AstroEpochs",
    "AstroUniverse",
    "AstroFrames",
    "AstroManeuvers",
    "AstroModels",
    "AstroCallbacks",
    "AstroProp",
    "AstroSolve",
    "Epicycle",
]

println("Testing all Epicycle packages via TestEnv from root...")
println("Repo root: $repo_root")
println("=" ^ 60)

failed_packages = String[]

for pkg in packages
    println("\n🧪 Testing $pkg...")
    pkg_path = joinpath(repo_root, pkg)
    if !isdir(pkg_path)
        println("⚠️  Package directory not found: $pkg_path")
        continue
    end

    try
        TestEnv.activate(pkg) do
            include(joinpath(pkg_path, "test", "runtests.jl"))
        end
        println("✅ $pkg tests passed")
    catch e
        println("❌ $pkg tests failed: $(sprint(showerror, e))")
        push!(failed_packages, pkg)
    end
end

println("\n" * "=" ^ 60)
println("TEST SUMMARY:")
println("=" ^ 60)

if isempty(failed_packages)
    println("🎉 All packages passed their tests!")
else
    println("❌ Failed packages: $(join(failed_packages, ", "))")
    println("📊 Passed: $(length(packages) - length(failed_packages))/$(length(packages))")
    error("Some packages failed their tests")
end