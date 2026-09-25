# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Runs the Epicycle tour.
#
#   julia --project=<environment> Epicycle/demo/run_demo.jl            # timed pauses
#   julia --project=<environment> Epicycle/demo/run_demo.jl --enter    # Enter advances
#   julia --project=<environment> Epicycle/demo/run_demo.jl --fast     # no pauses
#
# The environment needs the Epicycle packages and EpicycleIO.

include(joinpath(@__DIR__, "epicycle_demo.jl"))

EpicycleDemo.epicycle_demo(pace = "--fast" in ARGS ? 0.0 : 1.0,
                           wait_for_enter = "--enter" in ARGS)
