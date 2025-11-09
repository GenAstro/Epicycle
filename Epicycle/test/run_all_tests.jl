# CI test runner: develop local packages in this repo and run each package's tests
using Pkg, Logging

# Go back to repo root from Epicycle/test/
cd(joinpath(@__DIR__, "..", ".."))

# Use a temporary environment to avoid dependency resolution issues
# with the root Project.toml that lists all local packages
temp_env = mktempdir()
Pkg.activate(temp_env)

pkgs = [
    "AstroBase",
    "AstroStates",
    "AstroEpochs",
    "AstroUniverse",
    "AstroFrames",
    "AstroModels",
    "AstroManeuvers",
    "AstroFun",
    "AstroProp",
    "AstroSolve",
]

@info "Developing local packages into active environment"
for p in pkgs
    try
        @info "developing" p
        Pkg.develop(PackageSpec(path=p))
    catch e
        @warn "Pkg.develop failed for $p: $e"
    end
end

@info "Instantiate environment"
Pkg.instantiate()

# Ensure Coverage is available when running with --code-coverage
try
    import Coverage
catch
    @info "Adding Coverage package"
    Pkg.add("Coverage")
    import Coverage
end

# Run tests for each package in the same process. This mirrors your local workflow
succeeded = Dict{String,Bool}()
root = pwd()
for p in pkgs
    try
        @info "Running tests for" p
        pkg_test = joinpath(root, p, "test", "runtests.jl")
        if isfile(pkg_test)
            # Activate the package's environment to get access to [extras] test dependencies
            pkg_path = joinpath(root, p)
            Pkg.activate(pkg_path)
            Pkg.instantiate()  # Ensure test dependencies are available
            
            # cd into package to keep relative paths the same as local runs
            cd(joinpath(root, p)) do
                include(pkg_test)
            end
            succeeded[p] = true
        else
            @warn "No runtests.jl for $p at $pkg_test; marking as skipped"
            succeeded[p] = true
        end
    catch e
        @warn "Tests failed for $p: $e"
        succeeded[p] = false
    end
end

# Switch back to temp environment for coverage generation
Pkg.activate(temp_env)

# Write per-package LCOV files into cov/
mkpath("cov")
using Coverage: process_folder, LCOV
for p in pkgs
    try
        srcdir = joinpath(p, "src")
        if isdir(srcdir)
            outfile = joinpath("cov", string(p, ".lcov"))
            LCOV.writefile(outfile, process_folder(srcdir))
            @info "Wrote coverage for $p to $outfile"
        else
            @warn "No src/ for $p; skipping coverage"
        end
    catch e
        @warn "Coverage generation failed for $p: $e"
    end
end

# Exit with non-zero if any failed
if any(v -> v == false, values(succeeded))
    @error "Some package tests failed"
    exit(1)
else
    @info "All package tests passed"
    exit(0)
end