# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

"""
Module containing models for spacecraft maneuvers. 
"""
module AstroManeuvers

using LinearAlgebra

using EpicycleBase
using AstroEpochs
using AstroStates
using AstroFrames
using AstroModels: Spacecraft, to_posvel, set_posvel!
using AstroModels: HistorySegment, SpacecraftHistory, push_segment!

export ImpulsiveManeuver, maneuver!
import AstroFrames: AbstractAxes, VNB, Inertial

"""
    MANEUVER_AXES

Allowed maneuver axes types for ImpulsiveManeuver.
"""
const MANEUVER_AXES = (VNB, Inertial)

"""
    ImpulsiveManeuver{A<:AbstractAxes, T<:Real}

Represents an impulsive delta-v command defined in a local maneuver frame.

Fields
- axes::A        local maneuver frame (e.g., VNB() or Inertial())
- g0::T          standard gravity in m/s^2 (internally converted to km/s^2 where needed)
- Isp::T         specific impulse in seconds
- element1::T    delta-v component along axis 1 of `axes` (km/s)
- element2::T    delta-v component along axis 2 of `axes` (km/s)
- element3::T    delta-v component along axis 3 of `axes` (km/s)

# Notes:
- Use keyword constructor for convenience to avoid setting all values. 
- Delta-v components are interpreted in the provided `axes` frame and expressed in km/s.
- Numeric types are generic over Real to support AD and high-precision arithmetic.

# Examples
```jldoctest
using AstroManeuvers, AstroFrames

m = ImpulsiveManeuver(axes=VNB(), 
                      g0=9.81, Isp=300.0,
                      element1=0.010, 
                      element2=0.005, 
                      element3=-0.002)

# output
ImpulsiveManeuver(
  axes = VNB()
  g0 = 9.81
  Isp = 300.0
  element1 = 0.01
  element2 = 0.005
  element3 = -0.002)
```
"""
mutable struct ImpulsiveManeuver{A<:AbstractAxes, T<:Real}
    axes::A
    g0::T
    Isp::T
    element1::T
    element2::T
    element3::T
    function ImpulsiveManeuver(axes::A, g0::T, Isp::T, element1::T, element2::T, element3::T) where {A<:AbstractAxes, T<:Real}
        _validate_maneuver_inputs(axes, g0, Isp)
        new{A,T}(axes, g0, Isp, element1, element2, element3)
    end
end

""" 
    function _validate_maneuver_inputs(axes::AbstractAxes, g0::Real, Isp::Real)

Validates inputs for constructing an ImpulsiveManeuver.
"""
@inline function _validate_maneuver_inputs(axes::AbstractAxes, g0::Real, Isp::Real)
    if !any(T -> axes isa T, MANEUVER_AXES)
        allowed = join(string.(MANEUVER_AXES), ", ")
        throw(ArgumentError("ImpulsiveManeuver: axes must be one of ($allowed); got $(typeof(axes))"))
    end
    if !(g0 > 0)
        throw(ArgumentError("ImpulsiveManeuver: g0 must be > 0 m/s^2; got $(g0)"))
    end
    if !(Isp > 0)
        throw(ArgumentError("ImpulsiveManeuver: Isp must be > 0 s; got $(Isp)"))
    end
    return nothing
end

"""
    ImpulsiveManeuver(; axes=VNB(), g0=9.81, Isp=220.0, element1=0.0, element2=0.0, element3=0.0)

Construct an impulsive maneuver from keyword arguments with defaults.
"""
function ImpulsiveManeuver(;
    axes::AbstractAxes = VNB(),
    g0::Real = 9.81,
    Isp::Real = 220.0,
    element1::Real = 0.0,
    element2::Real = 0.0,
    element3::Real = 0.0,
)
    # Promote all numeric arguments to common type
    T = promote_type(typeof(g0), typeof(Isp), typeof(element1), typeof(element2), typeof(element3))
    ImpulsiveManeuver(axes, T(g0), T(Isp), T(element1), T(element2), T(element3))
end

"""
    Base.show(io::IO, sc::ImpulsiveManeuver)

Pretty-print an ImpulsiveManeuver in a human-readable, multi-line summary.
"""
function Base.show(io::IO, m::ImpulsiveManeuver)
    println(io, "ImpulsiveManeuver(")
    println(io, "  axes = ", m.axes)
    println(io, "  g0 = ", m.g0)
    println(io, "  Isp = ", m.Isp)
    println(io, "  element1 = ", m.element1)
    println(io, "  element2 = ", m.element2)
    println(io, "  element3 = ", m.element3, ")")
end

"""
    get_deltav_elements(m::ImpulsiveManeuver) -> (dv1, dv2, dv3)

Return the delta-v components in the maneuver axes (as defined by `m.axes`).

Arguments
- m::ImpulsiveManeuver    maneuver definition

# Returns
- (dv1, dv2, dv3)::Tuple  components in the maneuver frame axes
"""
get_deltav_elements(m::ImpulsiveManeuver) = (m.element1, m.element2, m.element3)

"""
    rot_mat_axes_to_inertial(::VNB, posvel) -> 3x3

Return the DCM from the maneuver frame to inertial.
"""
function rot_mat_axes_to_inertial(::VNB, posvel::AbstractVector)
    rot_mat_vnb_to_inertial(posvel)
end

"""
    rot_mat_axes_to_inertial(::Inertial, posvel) -> 3x3

Return the DCM from the maneuver frame to inertial.
"""
function rot_mat_axes_to_inertial(::Inertial, posvel::AbstractVector)
    T = eltype(posvel)
    return Matrix{T}(I, 3, 3)
end

"""
    get_deltav_inertial(m::ImpulsiveManeuver, sc) -> (dvi, dvj, dvk)

Returns maneuver components in the inertial frame.
"""
function get_deltav_inertial(m::ImpulsiveManeuver, sc)
    posvel = to_posvel(sc)
    R = rot_mat_axes_to_inertial(m.axes, posvel)
    T = eltype(posvel)
    dv_vec = T[m.element1, m.element2, m.element3]
    dv_inertial = R * dv_vec
    return Tuple(dv_inertial)
end

"""
    rot_mat_vnb_to_inertial(posvel::AbstractVector) -> 3x3

Build the VNB→Inertial DCM from position/velocity (km, km/s). 

Notes:
- V is along v,
- N is orbital normal r×v,
- B completes right-handed triad.
"""
function rot_mat_vnb_to_inertial(posvel::AbstractVector)
    r̄ = posvel[1:3]
    v̄ = posvel[4:6]
    nv = norm(v̄)
    nr = norm(r̄)
    if nv == 0 || nr == 0
        error("Cannot form VNB frame: zero-norm r or v")
    end
    v̂ = v̄ / nv
    n̄ = cross(r̄, v̄)
    nn = norm(n̄)
    nn == 0 && error("Cannot form VNB frame: r×v is zero (degenerate orbit)")
    n̂ = n̄ / nn
    b̂ = cross(v̂, n̂)
    return hcat(v̂, n̂, b̂)
end

"""
    compute_mass_used(m::ImpulsiveManeuver, initial_mass::Real, Isp::Real) -> Real

Computes mass consumed from impulsive maneuver using the rocket equation.
"""
function compute_mass_used(m::ImpulsiveManeuver, initial_mass::Real, Isp::Real)
    T = promote_type(typeof(initial_mass), typeof(m.g0), typeof(Isp), typeof(m.element1))
    dv = sqrt(T(m.element1)^2 + T(m.element2)^2 + T(m.element3)^2)
    g0_km = T(m.g0) / T(1000)
    denom = g0_km * T(Isp)
    if denom == zero(T)
        # Input validation should make this impossible to reach. Extra caution here.
        error("Invalid g0 or Isp; cannot compute mass used")  # COV_EXCL_LINE
    end
    return T(initial_mass) * (one(T) - exp(-dv / denom))
end

"""
    consume_fuel!(m::ImpulsiveManeuver, sc)

Subtract mass used by maneuver from the spacecraft mass.
"""
function consume_fuel!(m::ImpulsiveManeuver, sc::Spacecraft)
    total_mass = sc.mass
    sc.mass -= compute_mass_used(m, total_mass, m.Isp)
    return sc
end

"""
    maneuver!(sc::Spacecraft, m::ImpulsiveManeuver) -> Spacecraft

Apply an impulsive delta-v maneuver to a spacecraft, updating velocity, mass, and history.

Arguments
- sc::Spacecraft           Spacecraft being maneuvered
- m::ImpulsiveManeuver     Maneuver definition (axes, Δv components, g0, Isp)

# Notes:
- Model updates spacecraft velocity, mass, and history. 
- In places updates to spacecraft are performed.

# Returns
- sc::Spacecraft  The same spacecraft instance with updated state, mass, and history

# Examples
```jldoctest
using AstroManeuvers, AstroFrames, AstroModels

m = ImpulsiveManeuver(axes=Inertial(), 
                      Isp=300.0, 
                      element1=0.01, 
                      element2=0.0, 
                      element3=0.0)

sc = Spacecraft()
maneuver!(sc, m)

# output
Spacecraft: unnamed
  AstroEpochs.Time
    value  = 2015-09-21T12:23:12.000
    scale  = UTC()
    format = ISOT()
  OrbitState:
    statetype: AstroStates.Cartesian
  CartesianState:
    x   =  7000.00000000
    y   =     0.00000000
    z   =     0.00000000
    vx  =     0.01000000
    vy  =     7.50000000
    vz  =     0.00000000
  CoordinateSystem:
    origin: Earth
    axes: ICRFAxes
  Total Mass = 996.6078730003628 kg
```
"""
function maneuver!(sc::Spacecraft, m::ImpulsiveManeuver)
    Δv_inertial = get_deltav_inertial(m, sc)
    state = to_posvel(sc)
    state[4:6] .+= Δv_inertial
    set_posvel!(sc, state)
    consume_fuel!(m, sc)
    push_history_segment!(sc, m)
    return sc
end

"""
    push_history_segment!(sc::Spacecraft, m::ImpulsiveManeuver)

Append a new segment with maneuver to spacecraft history.
"""
function push_history_segment!(sc::Spacecraft, m::ImpulsiveManeuver)
    # Create single-point history segment for maneuver
    times = [sc.time]
    states = [CartesianState(to_posvel(sc))]
    segment = HistorySegment(times, states, sc.coord_sys, name="maneuver")
    push_segment!(sc.history, segment)
    return sc
end

"""
    promote(m::ImpulsiveManeuver, ::Type{T}) where T<:Real -> ImpulsiveManeuver{A,T}

Promote all Real fields of an ImpulsiveManeuver to type T for automatic differentiation.

# Arguments
- m::ImpulsiveManeuver: The maneuver to promote
- ::Type{T}: Target numeric type (e.g., ForwardDiff.Dual, BigFloat)

# Returns  
- ImpulsiveManeuver{A,T}: New maneuver with all Real fields promoted to type T

# Examples
```julia
using ForwardDiff
m = ImpulsiveManeuver(axes=VNB(), g0=9.81, Isp=300.0, element1=0.01, element2=0.02, element3=0.03)
m_dual = promote(m, ForwardDiff.Dual{Float64})
```
"""
function Base.promote(m::ImpulsiveManeuver, ::Type{T}) where T<:Real
    ImpulsiveManeuver(
        m.axes,           
        T(m.g0),          
        T(m.Isp),           
        T(m.element1),    
        T(m.element2),    
        T(m.element3),    
    )
end

end
