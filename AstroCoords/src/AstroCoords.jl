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

struct ICRFAxes <: AbstractAxes end
struct MJ2000Axes <: AbstractAxes end   
struct VNB <: AbstractAxes end
struct Inertial <: AbstractAxes end

mutable struct CoordinateSystem{O<:AbstractPoint, A<:AbstractAxes} <: AbstractCoordinateSystem
    origin::O
    axes::A
end

# Safe property accessor for show
@inline _maybe_get(x, s::Symbol) = (s in propertynames(x)) ? getfield(x, s) : nothing

"""
    Base.show(io::IO, ::MIME"text/plain", cs::CoordinateSystem)

Pretty-print a CoordinateSystem with safe handling of arbitrary origin types.
Prints origin name/type, axes type, and μ only if available on the origin.
"""
function Base.show(io::IO, ::MIME"text/plain", cs::CoordinateSystem)
    o = cs.origin
    a = cs.axes
    oname = _maybe_get(o, :name)
    omu   = _maybe_get(o, :mu)

    println(io, "CoordinateSystem:")
    print(io, "  origin = ")
    if oname !== nothing
        println(io, oname)
    end
    println(io, "  axes   = ", typeof(a))
end

# Delegate generic show to text/plain
function Base.show(io::IO, cs::CoordinateSystem)
    show(io, MIME"text/plain"(), cs)
end

end