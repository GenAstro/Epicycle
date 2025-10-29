# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    AstroCoords

Coordinate systems and reference frames for astrodynamics applications.

Provides types for defining coordinate systems with customizable origins and axes orientations,
supporting common reference frames like ICRF, J2000, and VNB for example .
"""
module AstroCoords

using AstroBase
using AstroUniverse

export AbstractAxes
export ICRFAxes, MJ2000Axes, VNB, Inertial
export AbstractCoordinateSystem
export CoordinateSystem

"""
    AbstractCoordinateSystem

Abstract type for coordinate system implementations.

All coordinate systems must define an origin and axes orientation.
See also: [`CoordinateSystem`](@ref).
"""
abstract type AbstractCoordinateSystem end

"""
    AbstractAxes

Abstract type for coordinate system axes definitions.

Concrete types specify the orientation of coordinate axes.
See also: [`ICRFAxes`](@ref), [`VNB`](@ref), [`Inertial`](@ref).
"""
abstract type AbstractAxes end

"""
    ICRFAxes

International Celestial Reference Frame axes.
"""
struct ICRFAxes <: AbstractAxes end

"""
    MJ2000Axes

Mean J2000.0 coordinate axes.
"""
struct MJ2000Axes <: AbstractAxes end   

"""
    VNB

Velocity-Normal-Binormal coordinate axes.

Trajectory-relative coordinate system where V is along velocity vector, 
N is normal to the orbital plane, and B completes the right-handed frame.
"""
struct VNB <: AbstractAxes end

"""
    Inertial

Inertial coordinate axes.
"""
struct Inertial <: AbstractAxes end

"""
    CoordinateSystem{O,A}(origin, axes)
    CoordinateSystem(origin, axes)

Coordinate system defined by an origin point and axes orientation.

# Arguments
- `origin::AbstractPoint`: Origin point of the coordinate system
- `axes::AbstractAxes`: Axes orientation specification

# Examples
```julia
# Earth-centered inertial coordinate system
cs = CoordinateSystem(earth, ICRFAxes())

# Spacecraft-relative VNB frame  
cs_vnb = CoordinateSystem(spacecraft_position, VNB())
```

See also: [`AbstractAxes`](@ref), [`ICRFAxes`](@ref), [`VNB`](@ref).
"""
mutable struct CoordinateSystem{O<:AbstractPoint, A<:AbstractAxes} <: AbstractCoordinateSystem
    origin::O
    axes::A        # COV_EXCL_LINE tested but inlined away so cov does not detect
end

# Safe property accessor for show
@inline _maybe_get(x, s::Symbol) = (s in propertynames(x)) ? getfield(x, s) : nothing

"""
    Base.show(io::IO, ::MIME"text/plain", cs::CoordinateSystem)

Display a CoordinateSystem showing origin name and axes type. 
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

""" 
    function Base.show(io::IO, cs::CoordinateSystem)

Delegate generic show to text/plain
"""
function Base.show(io::IO, cs::CoordinateSystem)
    show(io, MIME"text/plain"(), cs)
end

end