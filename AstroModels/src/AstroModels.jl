# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

__precompile__()

""" 
Module containing physical models such as spacecraft. 
"""
module AstroModels

using AstroEpochs
using AstroStates
using AstroCoords
using AstroUniverse

# Import commonly used types to avoid qualification
using AstroCoords: ICRFAxes, CoordinateSystem
using AstroUniverse: earth

export Spacecraft, get_state, to_posvel, set_posvel!
export push_history_segment!

import AstroBase: AbstractPoint

include("spacecraft.jl")

end