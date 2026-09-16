# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# Build one package's documentation, the way CI does.
#
#     julia --startup-file=no --project=. ci/build_docs.jl AstroRoutines
#     julia --startup-file=no --project=. ci/build_docs.jl AstroRoutines AstroStates
#     julia --startup-file=no --project=. ci/build_docs.jl --all
#
# This is `ci/build_epicycle.jl`'s documentation loop for a named package rather than for all
# eleven, so a docs failure can be reproduced and fixed in minutes instead of a full CI run.
#
# Two things it does that a bare `include` of make.jl does not. It develops the package first, so
# an unregistered package resolves from its path rather than from the registry. And it reports the
# failure rather than letting it escape as a stack trace from inside Documenter, which is what
# makes a `@example` failure legible.
#
# Run it from the repository root with the root project active, which is the environment CI builds
# documentation in — not `epicycle-dev`. That difference matters: a package whose docs need a
# dependency the root project lacks builds locally against the warm environment and fails in CI.

using Pkg

const REPO_ROOT = dirname(@__DIR__)

# The packages CI documents, in CI's order.
const ALL_PACKAGES = [
    "AstroRoutines", "EpicycleBase", "AstroStates", "AstroEpochs", "AstroUniverse",
    "AstroFrames", "AstroModels", "AstroManeuvers", "AstroCallbacks",
    "AstroProp", "AstroSolve", "Epicycle",
]

"""Build `pkg`'s documentation. Returns `true` on success."""
function build_one(pkg::AbstractString)
    make_path = joinpath(REPO_ROOT, pkg, "docs", "make.jl")
    if !isfile(make_path)
        println("  ⚠️  $pkg has no docs/make.jl — skipping")
        return true
    end

    println("\n", "="^70)
    println("📖 $pkg")
    println("="^70)

    try
        # A package that is not registered resolves only from its path, so develop before loading —
        # but only when it is not already developed. `Pkg.develop` rewrites Project.toml and can
        # force a re-resolve, and re-resolving the root environment recompiles it, which is minutes
        # every run rather than once.
        if !haskey(Pkg.project().dependencies, pkg)
            println("  🔗 developing $pkg (not yet in this project)")
            Pkg.develop(path = joinpath(REPO_ROOT, pkg))
        end
        include(make_path)
        println("  ✅ $pkg docs built")
        return true
    catch e
        println("  ❌ $pkg docs FAILED")
        println("  ", sprint(showerror, e))
        return false
    end
end

pkgs = isempty(ARGS) ? String[] : (("--all" in ARGS || "-a" in ARGS) ? ALL_PACKAGES : String.(ARGS))

if isempty(pkgs)
    println("usage: julia --startup-file=no --project=. ci/build_docs.jl <Package> [Package...]")
    println("       julia --startup-file=no --project=. ci/build_docs.jl --all")
    println("\nknown packages:\n  ", join(ALL_PACKAGES, "\n  "))
    exit(2)
end

# Instantiate only if something is missing; on a settled environment this is a no-op, and calling
# it unconditionally is what makes a docs rebuild feel like a full CI run.
if !isfile(joinpath(REPO_ROOT, "Manifest.toml"))
    println("🔧 No manifest — resolving the root project (one-off, and slow)...")
    Pkg.instantiate()
end

failed = String[]
for pkg in pkgs
    build_one(pkg) || push!(failed, pkg)
end

println("\n", "="^70)
if isempty(failed)
    println("✅ documentation built for: ", join(pkgs, ", "))
else
    println("❌ documentation failed for: ", join(failed, ", "))
    exit(1)
end
