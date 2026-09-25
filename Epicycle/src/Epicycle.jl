# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

module Epicycle

using Reexport

# Import all the individual packages (now available via Pkg.develop)
@reexport using EpicycleBase
@reexport using AstroRoutines
@reexport using AstroStates
@reexport using AstroEpochs
@reexport using AstroUniverse
@reexport using AstroFrames
@reexport using AstroModels
@reexport using AstroManeuvers
@reexport using AstroCallbacks
@reexport using AstroProp
@reexport using AstroSolve
@reexport using EpicycleIO

# Walking a shipped example a step at a time, for someone learning the interface rather than
# reading it. `Epicycle.example_names()` lists them and `Epicycle.run_example(name)` runs one.
# Qualified rather than exported: a user reaches for this once, not in every script.
include("examples.jl")
using .Examples: example_names, run_example, tutorial_names, run_tutorial

# Re-export the packages a script needs alongside Epicycle. The four OrdinaryDiffEq solver
# packages hold every integrator named anywhere in Epicycle (Vern9, Vern7, Tsit5, DP8, RK4);
# SciMLBase supplies ODEProblem and the callbacks, and CommonSolve supplies `solve`.
# The OrdinaryDiffEq metapackage would bring 176 packages where these bring 91; the eighty-five
# it adds are stiff and specialist solvers, with NonlinearSolve and LinearSolve behind them, and
# nothing here calls any of them.
@reexport using SciMLBase
@reexport using CommonSolve
@reexport using OrdinaryDiffEqTsit5          # Tsit5
@reexport using OrdinaryDiffEqVerner         # Vern7, Vern9
@reexport using OrdinaryDiffEqHighOrderRK    # DP8
@reexport using OrdinaryDiffEqLowOrderRK     # RK4
@reexport using SNOW
@reexport using NLsolve
@reexport using SPICE

# AstroSolve and CommonSolve both export solve!. Re-exporting both leaves the name ambiguous, and
# a script that calls solve! gets an UndefVarError naming two packages. In Epicycle it means the
# solver's, so the umbrella says which.
using AstroSolve: solve!
export solve!

# Drawing is EpicycleIO's: Plotly for data and Cesium for 3D through `orbitview`. It loads with
# the umbrella, because a mission design system whose graphics have to be asked for separately is
# a system whose graphics people do not find. It costs about 0.9 s on top of a 12.8 s
# `using Epicycle`, and its whole dependency tree is PlotlyBase, JSON, Scratch and Sockets, since
# the rendering happens in the browser and there is no graphics stack behind it.
#
# Its trace wrappers are named for what they draw, `xyplot` included, which is why that one is not
# called `plot`. Every plotting package in Julia claims `plot`, `contour`, `heatmap`, `surface`,
# `histogram` and `bar`, and none of them coordinate, so a user who wants Plots.jl alongside this
# renames one side on import: `using Plots: plot as pplot`.
#
# The native GLMakie window that used to live in this package is out of the release. GLMakie was
# a direct dependency of this umbrella until 2026-09-13, the largest single tree in the stack,
# costing about 5 s of every `using Epicycle` and paid by the majority of users who never opened
# a window.

end