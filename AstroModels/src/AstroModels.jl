# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

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

include("spacecraft.jl")

end