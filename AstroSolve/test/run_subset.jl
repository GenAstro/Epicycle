# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Runs some of AstroSolve's test files instead of all of them, with the same setup runtests.jl
# gives them. The whole suite in one process can exhaust memory on a machine that also holds a
# REPL, so a change is checked a few files at a time.
#
#   julia --project=<environment> -e 'using TestEnv, AstroSolve;
#       TestEnv.activate("AstroSolve") do; include(joinpath(pkgdir(AstroSolve), "test", "run_subset.jl")); end' \
#       test_correctness_udu.jl test_correctness_rts.jl
#
# File names come from ARGS and are resolved against this directory. The list of files that make
# up the whole suite is runtests.jl, which this does not replace.
#
# Some files use helpers defined in a file runtests.jl includes earlier, so name those first:
#   test_correctness_batch_priors.jl  needs  test_correctness_ekf.jl           (_setup2)
#   test_correctness_oc_manager.jl    needs  test_correctness_brachistochrone.jl (_br_phase)
#                                     and    test_correctness_sims_flanagan.jl   (_sf_phase)

using LinearAlgebra

using EpicycleBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroFrames
using AstroModels
using AstroManeuvers
using AstroCallbacks
using AstroProp
using AstroSolve

# The shared fixtures runtests.jl defines before its includes, which some files read.
posvel = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
sat = Spacecraft(
    state=CartesianState(posvel),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

using Test

isempty(ARGS) && throw(ArgumentError(
    "run_subset.jl needs at least one test file name, such as test_correctness_udu.jl"))

for file in ARGS
    path = joinpath(@__DIR__, file)
    isfile(path) || throw(ArgumentError("run_subset.jl: no test file $(file) in $(@__DIR__)"))
    println("── ", file)
    include(path)
end
