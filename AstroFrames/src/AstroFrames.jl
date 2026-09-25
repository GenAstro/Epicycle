# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    AstroFrames

Coordinate systems and reference frames for astrodynamics applications.

Provides types for defining coordinate systems with customizable origins and
axes orientations. Axes types follow community naming (`ICRF`, `GCRF`, `ITRF`,
etc.); some are restricted to specific origins per the origin-coupling rules
declared in `origin_coupling.jl`.

See also:
- [`CoordinateSystem`](@ref) — pairs an origin and axes.
- [`ICRF`](@ref), [`GCRF`](@ref), [`ITRF`](@ref), [`CelestialBodyFixed`](@ref) — axes types.
"""
module AstroFrames

using EpicycleBase
using AstroUniverse
using StaticArrays
using LinearAlgebra: I, norm, cross, dot
using AstroStates: to_vector, AbstractOrbitState, OrbitState, state_tag_to_type
import AstroStates: CartesianState   # methods added in coordinate.jl

export AbstractAxes
export ICRF, GCRF, CIRS, TIRS, ITRF
export ICRFAxes, MJ2000Axes   # deprecated; see axes_types.jl
export MJ2000Eq, MODEq, TODEq, MODEc, TODEc, PEF, TEME
export MJ2000Ec
export MoonPA, MoonME
export CelestialBodyFixed
export RIC, LVLH, VNB
export Inertial
export AbstractCoordinateSystem, CoordinateSystem
export axes_rotation, origin_translation, edge_theory

# The extension contract - what you write a frame of your own against.
export EpochScales, epoch_tdb, epoch_tt, epoch_utc
export hub_axes, valid_origin, needs_reference_orbit

# Conversion for callers with no Spacecraft - a state, its frame, its epoch.
export Coordinate, state_of, frame_of, epoch_of

# Named frames — see named_frames.jl.
export EarthMJ2000Eq, EarthMJ2000Ec, EarthFixed, EarthTODEq, EarthICRF
export MoonFixed, MoonPrincipalAxes, SunMJ2000Ec

"""
    AbstractCoordinateSystem

Abstract type for coordinate system implementations. All coordinate systems
must define an origin and axes orientation.

See also: [`CoordinateSystem`](@ref).

# Example

```jldoctest
EarthMJ2000Eq isa AbstractCoordinateSystem

# output
true
```
"""
abstract type AbstractCoordinateSystem end

"""
    AbstractAxes

Abstract type for coordinate system axes definitions. Concrete subtypes
specify the orientation of the coordinate axes.

See also: [`ICRF`](@ref), [`GCRF`](@ref), [`ITRF`](@ref), [`CelestialBodyFixed`](@ref).

# Example

```jldoctest
ICRF() isa AbstractAxes, ITRF() isa AbstractAxes

# output
(true, true)
```
"""
abstract type AbstractAxes end

include("axes_types.jl")
include("origin_coupling.jl")
include("rotation_math.jl")
include("axes_rotation.jl")
include("orbit_relative.jl")
include("origin_translation.jl")

"""
    CoordinateSystem{O,A}(origin, axes)
    CoordinateSystem(origin, axes)

Coordinate system defined by an origin point and axes orientation.

# Fields
- `origin::AbstractPoint`: Origin point (e.g. a `CelestialBody`).
- `axes::AbstractAxes`: Axes orientation specification.

Both are also the constructor's arguments, in that order.

# Origin coupling

Some axes types are only physically meaningful at specific origins:

- Earth-restricted (`GCRF`, `CIRS`, `TIRS`, `ITRF`, `MODEq`, `TODEq`, `MODEc`, `TODEc`, `PEF`)
  — require an Earth origin.
- Moon-restricted (`MoonPA`, `MoonME`) — require a Moon origin.
- `CelestialBodyFixed{OT}` — requires an origin of type `OT` (the body encoded in
  the axes type parameter).

Invalid combinations throw an `ArgumentError` at construction with a
domain-language message.

# Examples
```julia
# Earth-centered ICRF coordinate system.
cs = CoordinateSystem(earth, ICRF())

# Earth-fixed frame with EOP.
cs_itrs = CoordinateSystem(earth, ITRF())

# Mars body-fixed (IAU 2015) — origin fills in the body:
cs_mars = CoordinateSystem(mars, CelestialBodyFixed())

# Invalid — CoordinateSystem(sun, ITRF()) throws:
#   "ITRF axes require an Earth origin (use `earth`); got Sun."
```
"""
mutable struct CoordinateSystem{O<:AbstractPoint, A<:AbstractAxes} <: AbstractCoordinateSystem
    origin::O
    axes::A

    function CoordinateSystem{O,A}(origin::O, axes::A) where {O<:AbstractPoint, A<:AbstractAxes}
        valid_origin(axes, origin) || throw(ArgumentError(
            "$(A) axes require $(_axes_family_hint(A)); got origin $(_origin_display(origin))."))
        return new{O,A}(origin, axes)
    end
end

CoordinateSystem(origin::O, axes::A) where {O<:AbstractPoint, A<:AbstractAxes} =
    CoordinateSystem{O,A}(origin, axes)

# `CelestialBodyFixed()` (unresolved sentinel) — resolve NAIF ID from the origin.
#
# Any body the universe knows how to orient is allowed, not a fixed list: a
# user who registers an orientation model for their own body gets body-fixed
# axes for it with no further ceremony.
function CoordinateSystem(origin::AbstractPoint, ::CelestialBodyFixed{0})
    n = _naifid_of(origin)
    n === nothing && throw(ArgumentError(
        "CelestialBodyFixed axes need an origin with a NAIF ID; got $(_origin_display(origin))."))
    AstroUniverse.orientation_model(n)   # raises, naming what to do, if there is none
    return CoordinateSystem(origin, CelestialBodyFixed{n}())
end

_origin_display(origin) = (:name in propertynames(origin)) ? origin.name : string(origin)

# Safe property accessor for `show`
@inline _maybe_get(x, s::Symbol) = (s in propertynames(x)) ? getfield(x, s) : nothing

"""
    Base.show(io::IO, ::MIME"text/plain", cs::CoordinateSystem)

Display a CoordinateSystem showing origin name and axes type.
"""
function Base.show(io::IO, ::MIME"text/plain", cs::CoordinateSystem)
    o = cs.origin
    a = cs.axes
    oname = _maybe_get(o, :name)

    println(io, "CoordinateSystem:")
    print(io, "  origin = ")
    if oname !== nothing
        println(io, oname)
    else
        println(io, typeof(o))
    end
    print(io, "  axes   = ", typeof(a))
end

"""
    Base.show(io::IO, cs::CoordinateSystem)

Delegate generic show to text/plain.
"""
Base.show(io::IO, cs::CoordinateSystem) = show(io, MIME"text/plain"(), cs)

# Needs `CoordinateSystem`, defined above.
include("coordinate.jl")
include("named_frames.jl")

end
