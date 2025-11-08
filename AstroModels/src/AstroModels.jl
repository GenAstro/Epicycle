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
export push_history_segment!

import AstroBase: AbstractPoint

include("spacecraft.jl")

end