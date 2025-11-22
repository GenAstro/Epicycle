# Frame Transform Prototype
# Prototyping frame transformation design before integration into AstroFrames

using StaticArrays
using LinearAlgebra
using SatelliteToolboxTransformations
using AstroModels
using AstroUniverse
using AstroEpochs

# =============================================================================
# Core Type Hierarchy
# =============================================================================

abstract type AbstractReferenceSystem end
struct FK5 <: AbstractReferenceSystem end
struct IAU2006 <: AbstractReferenceSystem end
struct Inherited <: AbstractReferenceSystem end

abstract type AbstractOrientationType end
struct Inertial <: AbstractOrientationType end
struct Fixed <: AbstractOrientationType end
struct Derived <: AbstractOrientationType end

abstract type AbstractAxes end

# =============================================================================
# Axes Types
# =============================================================================

struct GCRF <: AbstractAxes end
struct ITRF <: AbstractAxes end
struct LVLH <: AbstractAxes end

# =============================================================================
# CoordinateFrame and Coordinate Structs
# =============================================================================

"""
    CoordinateFrame

Defines a reference frame with an origin point and orientation axes.
"""
struct CoordinateFrame
    origin::AbstractPoint
    axes::AbstractAxes
end

"""
    Coordinate{T<:Real}

Represents a state in space-time with position, velocity, acceleration,
time, and coordinate frame. Parametric type T supports Float64, ForwardDiff.Dual, etc.
for automatic differentiation compatibility.
"""
struct Coordinate{T<:Real}
    pos::SVector{3,T}
    vel::SVector{3,T}
    acc::SVector{3,T}
    time::Time
    frame::CoordinateFrame
end

# Convenience constructors
function Coordinate(pos::SVector{3,T}, vel::SVector{3,T}, time::Time, frame::CoordinateFrame) where T<:Real
    Coordinate{T}(
        pos,
        vel,
        SVector{3,T}(zero(T), zero(T), zero(T)),
        time,
        frame
    )
end

# Allow construction from any vector-like input
function Coordinate(pos, vel, time::Time, frame::CoordinateFrame)
    T = promote_type(eltype(pos), eltype(vel))
    Coordinate(
        SVector{3,T}(pos),
        SVector{3,T}(vel),
        SVector{3,T}(zero(T), zero(T), zero(T)),
        time,
        frame
    )
end
