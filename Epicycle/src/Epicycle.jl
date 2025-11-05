# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

module Epicycle

using Reexport

# Import all the individual packages (now available via Pkg.develop)
@reexport using AstroBase
@reexport using AstroStates
@reexport using AstroEpochs
@reexport using AstroUniverse
@reexport using AstroCoords
@reexport using AstroModels
@reexport using AstroMan
@reexport using AstroFun
@reexport using AstroProp
@reexport using AstroSolve

# Rexport packages required to use Epicycle
@reexport using OrdinaryDiffEq
@reexport using SNOW
@reexport using NLsolve
@reexport using SPICE

end