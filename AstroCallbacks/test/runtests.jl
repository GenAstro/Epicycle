# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

using Test
using LinearAlgebra

using AstroStates
using AstroEpochs
using EpicycleBase
using AstroUniverse
using AstroFrames
using AstroModels: Spacecraft, get_state 
using AstroManeuvers
using AstroCallbacks 

include("runtests_orbitcalcs.jl")
include("test_correctness_element_setters.jl")
include("test_correctness_quantity_readers.jl")
include("test_correctness_new_quantities.jl")
include("test_correctness_output_partials.jl")
include("test_correctness_history.jl")
include("test_correctness_solver_specs.jl")
include("test_correctness_setters.jl")
include("runtests_maneuvercalcs.jl")
include("runtests_bodycalcs.jl")
include("runtests_inputvalidation.jl")
include("runtests_infrastructure.jl")
include("test_correctness_frame_aware_quantities.jl")

nothing