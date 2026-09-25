# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

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
using AstroFrames: ICRF, CoordinateSystem
using AstroFrames: AbstractAxes, GCRF, ITRF, Coordinate
using AstroUniverse: CelestialBody
using AstroStates: CartesianState, to_vector
using LinearAlgebra: norm, dot
using StaticArrays: SVector
using AstroUniverse: earth

export Spacecraft, get_state, to_posvel, set_posvel!, total_mass
export AbstractDragGeometry, SphericalDrag
export AbstractSRPGeometry, SphericalSRP
export CADModel
export HistorySegment, SpacecraftHistory
export push_segment!
export to_float64
export Mass
export GroundStation, AbstractGeodeticReference, Ellipsoid, is_visible

import EpicycleBase: AbstractPoint, AbstractParamTag, get_field, set_field!

include("cadmodel.jl")
include("spacecraft_history.jl")
include("drag_geometry.jl")
include("srp_geometry.jl")
include("spacecraft.jl")
include("ground_station.jl")
include("tags.jl")
include("frame_conversion.jl")

end