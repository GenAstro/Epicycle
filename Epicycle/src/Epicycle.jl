# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

module Epicycle

using Reexport

# Import all the individual packages (now available via Pkg.develop)
@reexport using EpicycleBase
@reexport using AstroStates
@reexport using AstroEpochs
@reexport using AstroUniverse
@reexport using AstroFrames
@reexport using AstroModels
@reexport using AstroManeuvers
@reexport using AstroCallbacks
@reexport using AstroProp
@reexport using AstroSolve

# Rexport packages required to use Epicycle
@reexport using OrdinaryDiffEq
@reexport using SNOW
@reexport using NLsolve
@reexport using SPICE

# Graphics module (3D visualization)
using GLMakie
using FileIO
using MeshIO
using GeometryBasics
using Random
using LinearAlgebra
using Colors: RGB, N0f8

# Include graphics source files
include("graphics/view3d.jl")
include("graphics/trajectories.jl")
include("graphics/bodies.jl")
include("graphics/stars.jl")
include("graphics/planes.jl")
include("graphics/labels.jl")
include("graphics/models.jl")

# Export graphics API
export View3D, add_spacecraft!, display_view

end