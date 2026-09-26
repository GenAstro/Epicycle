# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using LinearAlgebra
#using OrdinaryDiffEq

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

# Name `solve!` explicitly. Six loaded modules export it — AstroSolve, Epicycle, CommonSolve,
# SciMLBase, DiffEqBase and OrdinaryDiffEqCore — and Julia refuses to choose among implicit
# exports, so an unqualified call raises `UndefVarError: solve! not defined in Main`. Standalone
# the suite is fine, because none of the others is in Main; test_all_packages.jl includes all
# thirteen into one Main, so whether this bites depends on which suite ran first. An explicit
# import wins over every implicit export and makes the order irrelevant.
using AstroSolve: solve!

# Create a spacecraft posvel SolverVariable
#
# There was a top-level `time = Time(...)` here. It was never used — the Spacecraft below
# builds its own — and at column zero it bound `time` in Main, shadowing `Base.time` for every
# suite that runs after this one. test_all_packages.jl includes all thirteen into one Main, so
# EpicycleIO's `_LAST_REQUEST[] = time()` then tried to call a timestamp and errored.
posvel = [7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]
sat = Spacecraft(
    state=CartesianState(posvel), 
    time=Time("2015-09-21T12:23:12", TAI(), ISOT())
    )

using Test

include("runtests_construction_solvervariable.jl")
include("runtests_setget_solvervariable.jl")
include("test_correctness_vary.jl")
include("test_correctness_constraint.jl")
include("runtests_constraint.jl")
include("runtests_events.jl")
include("runtests_events_deltav.jl")
include("runtests_sequence.jl")
include("runtests_sequence_manager.jl")
include("runtests_gaps.jl")
include("runtests_optimize.jl")
include("runtests_sequence_geotransfer.jl")
include("test_correctness_target_block.jl")
include("test_correctness_vary_elements.jl")
include("runtests_show.jl")
include("runtests_sequence_report_enhanced.jl")
include("runtests_solve_history.jl")
include("test_correctness_hermite_simpson.jl")
include("test_correctness_jacobian_sparsity.jl")




include("test_correctness_udu.jl")

include("test_correctness_process_noise.jl")

include("test_correctness_measurements.jl")

include("test_correctness_spring_mass.jl")

include("test_correctness_ekf.jl")

include("test_correctness_rts.jl")

include("test_correctness_od_batch.jl")
include("test_correctness_simulate.jl")

include("test_correctness_od_iterated.jl")

include("test_correctness_tdm_io.jl")

include("test_correctness_batch_priors.jl")

include("test_correctness_brachistochrone.jl")

include("test_correctness_trajectory_tracking.jl")

include("test_correctness_bolza.jl")

include("test_correctness_parameter_id.jl")

include("test_correctness_moon_landing.jl")

include("test_correctness_obstacle_avoidance.jl")

include("test_correctness_sims_flanagan.jl")

include("test_correctness_shooting_sequence.jl")

include("test_correctness_mga_ndsms.jl")

# Last of the collocation checks: it reuses the problems loaded above.
include("test_correctness_jacobian_sparsity_differential.jl")

include("test_correctness_oc_manager.jl")
include("test_correctness_mixed_links.jl")

include("test_correctness_ad_fallback.jl")
include("test_correctness_ad_jacobian_paths.jl")
include("test_correctness_show_collocation.jl")

include("test_correctness_solve_vocabulary.jl")
include("test_correctness_constraint_scale.jl")
include("test_correctness_autodiff_generators.jl")
