# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

""" 
Module containing physical models such as spacecraft. 
"""
module AstroModels

using AstroEpochs
using AstroStates
using AstroFrames
using AstroUniverse

# Import commonly used types to avoid qualification
using AstroFrames: ICRFAxes, CoordinateSystem
using AstroUniverse: earth

export Spacecraft, get_state, to_posvel, set_posvel!
export CADModel
export HistorySegment, SpacecraftHistory
export push_state!, push_segment!, new_segment!
export to_float64

import EpicycleBase: AbstractPoint

include("cadmodel.jl")
include("spacecraft_history.jl")
include("spacecraft.jl")

end