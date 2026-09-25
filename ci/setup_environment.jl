#!/usr/bin/env julia
"""
CI Setup Environment Script

Sets up the Julia environment for CI by:
1. Activating the project
2. Instantiating dependencies
3. Developing all workspace packages

This script should be run first in the CI pipeline.
"""

println("🔧 Setting up CI environment...")

using Pkg
# The package list is not written here. It is read from the workspace `projects` entry in the
# repo-root Project.toml, which is the one place that has to be right for Julia to resolve at
# all. Eight hand-maintained copies of this list is how AstroRoutines and EpicycleIO came to be
# missing from CI after they moved in, and how CairoMakie survived in the workspace after
# graphics left the umbrella.
using TOML
const WORKSPACE_PACKAGES = TOML.parsefile(
    joinpath(dirname(@__DIR__), "Project.toml"))["workspace"]["projects"]

# Activate the main project
println("📦 Activating project environment...")
Pkg.activate(".")

# Develop the workspace packages first. A package in the root [deps] that is not in any registry
# is looked up in the registries by `instantiate` and fails there — "expected package EpicycleIO
# to be registered" — before the develop loop below could have supplied it from a path. Developing
# first puts the path in the manifest, so the resolve that follows has somewhere to find it.
packages = WORKSPACE_PACKAGES

println("🔗 Developing workspace packages...")
for pkg in packages
    println("  → Developing $pkg...")
    Pkg.develop(path=pkg)
end

# Instantiate to get the registered dependencies
println("📥 Installing dependencies...")
Pkg.instantiate()

# Resolve any dependency conflicts
println("🎯 Resolving dependencies...")
Pkg.resolve()

println("✅ Environment setup complete!")
println("📋 Package status:")
Pkg.status()

