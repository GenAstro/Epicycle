# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

__precompile__()
module AstroCoords

using AstroBase
using AstroUniverse

abstract type AbstractCoordinateSystem end
abstract type AbstractAxes end

export AbstractAxes

# TODO Rename getting rid of "Axes" 
export ICRFAxes, MJ2000Axes, VNB, Inertial

export AbstractCoordinateSystem
export CoordinateSystem

mutable struct CoordinateSystem{O<:AbstractPoint, A<:AbstractAxes} <: AbstractCoordinateSystem
    origin::O
    axes::A
end

end