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

"""
Module containing orbital state representations and conversions.
"""
module AstroStates

# TODO validate state type inputs and call from OrbitState and AbstractOrbitState constructors

using AstroBase

abstract type AbstractOrbitState <: AbstractState end

import Base: show
import AstroBase: AbstractOrbitStateType

using LinearAlgebra
using Printf

# Export public types 
export AbstractState, AbstractOrbitState, AbstractOrbitStateType 
export CartesianState, KeplerianState, SphericalRADECState, ModifiedEquinoctialState
export OutGoingAsymptoteState, IncomingAsymptoteState, ModifiedKeplerianState, SphericalAZIFPAState
export EquinoctialState, AlternateEquinoctialState

# Export conversion functions 
export cart_to_kep, kep_to_cart, sphradec_to_cart, cart_to_sphradec, cart_to_mee, mee_to_cart
export outasymptote_to_kep, cart_to_outasymptote, to_vector, cart_to_inasymptote
export inasymptote_to_kep, sphazfpa_to_cart, cart_to_sphazfpa, equinoctial_to_alt_equinoctial
export alt_equinoctial_to_equinoctial, equinoctial_to_cart, cart_to_equinoctial
export modkep_to_kep, kep_to_modkep
export _marker_type
export state_tag_to_type, state_type_to_tag

struct Cartesian <: AbstractOrbitStateType end
struct Keplerian <: AbstractOrbitStateType end
struct Equinoctial <: AbstractOrbitStateType end
struct SphericalRADEC <: AbstractOrbitStateType end
struct SphericalAZIFPA <: AbstractOrbitStateType end
struct ModifiedEquinoctial <: AbstractOrbitStateType end
struct OutGoingAsymptote <: AbstractOrbitStateType end
struct IncomingAsymptote <: AbstractOrbitStateType end
struct ModifiedKeplerian <: AbstractOrbitStateType end
struct AlternateEquinoctial <: AbstractOrbitStateType end

export AbstractOrbitStateType, Cartesian, Keplerian, Equinoctial, SphericalRADEC, SphericalAZIFPA
export ModifiedEquinoctial, OutGoingAsymptote, IncomingAsymptote,ModifiedKeplerian, AlternateEquinoctial
export OrbitState, to_state

include("cart_to_kep.jl")
include("kep_to_cart.jl")

include("sphazfpa_to_cart.jl")
include("cart_to_sphazfpa.jl")

include("sphradec_to_cart.jl")
include("cart_to_sphradec.jl")

include("cart_to_outasymptote.jl")
include("outasymptote_to_kep.jl")

include("cart_to_inasymptote.jl")
include("inasymptote_to_kep.jl")

include("cart_to_mee.jl")
include("mee_to_cart.jl")

include("equinoctial_to_alt_equinoctial.jl")
include("alt_equinoctial_to_equinoctial.jl")

include("equinoctial_to_cart.jl")
include("cart_to_equinoctial.jl")

include("modkep_to_kep.jl")
include("kep_to_modkep.jl")

# =============================================================================
# State Type Definitions
# =============================================================================
"""
    CartesianState(posvel)

Mutable Cartesian state type. Parameterized for AD / alternate precisions.

# Fields
- `posvel::Vector{T}`: 6-element position/velocity vector (x,y,z,vx,vy,vz)

# Notes
- Distance/time units must be consistent with the gravitational parameter `μ` used elsewhere.
"""
mutable struct CartesianState{T} <: AbstractOrbitState
    posvel::Vector{T}   
    numvars::Int

    function CartesianState(state::Vector{T}) where {T<:Real}
        @assert length(state) == 6 "CartesianState must be a 6-element vector"
        new{T}(state, 6)
    end
end

function show(io::IO, state::CartesianState)
    posvel = state.posvel
    println(io, "CartesianState:")
    println(io, @sprintf("  x   = %14.8f", posvel[1]))
    println(io, @sprintf("  y   = %14.8f", posvel[2]))
    println(io, @sprintf("  z   = %14.8f", posvel[3]))
    println(io, @sprintf("  vx  = %14.8f", posvel[4]))
    println(io, @sprintf("  vy  = %14.8f", posvel[5]))
    println(io, @sprintf("  vz  = %14.8f", posvel[6]))
end

"""
    KeplerianState(a, e, i, raan, aop, ta)

Keplerian orbital elements representation using classical osculating elements.

# Units
- Distance/time units must be consistent with the gravitational parameter `μ` used in conversions.
- All angular quantities are in **radians**.

# Fields (all `::T` where `T<:Real`)
- `sma`: Semi-major axis (must be nonzero). Defines orbit size and energy.
  * If sma > 0: elliptic orbit (bound). If sma < 0: hyperbolic orbit (unbound).
  * For elliptic orbits, sma is half the major axis length.
- `ecc`: Eccentricity. Defines orbit shape.
  * ecc = 0: circular orbit. 0 < ecc < 1: elliptical orbit.
  * ecc = 1: parabolic orbit. ecc > 1: hyperbolic orbit.
  * Valid range: [0, ∞), but ecc ≈ 1 results in infinite sma.
- `inc`: Inclination (rad). Angle between orbit plane and reference xy-plane.
  * Range: [0, π]. If inc < π/2: prograde orbit. If inc > π/2: retrograde orbit.
- `raan`: Right ascension of ascending node (rad). Orients the orbit plane.
  * Range: [0, 2π). Angle from +x axis to ascending node, measured in xy-plane.
  * Defines where orbit plane intersects reference plane (ascending crossing).
  * Undefined for equatorial orbits (inc ≈ 0 or π).
- `aop`: Argument of periapsis (rad). Orients the orbit within its plane.
  * Range: [0, 2π). Angle from ascending node to periapsis point.
  * Defines orientation of orbit ellipse within the orbital plane.
  * Undefined for circular orbits (ecc ≈ 0).
- `ta`: True anomaly (rad). Spacecraft position within the orbit.
  * Range: [0, 2π). Angle from periapsis to current spacecraft position.
  * ta = 0: at periapsis. ta = π: at apoapsis (for elliptic orbits).
  * Undefined for circular orbits (ecc ≈ 0).

# Notes
- Parametric so automatic differentiation and high-precision types are supported.
- Classical Keplerian elements have well-known singularities for special cases.

# Examples
```julia
k = KeplerianState(7000.0, 0.01, π/4, 0.0, 0.0, π/3)
```
"""
struct KeplerianState{T<:Real} <: AbstractOrbitState
    sma::T
    ecc::T
    inc::T
    raan::T
    aop::T
    ta::T
end

function show(io::IO, state::KeplerianState)
    println(io, "KeplerianState:")
    println(io, @sprintf("  sma        = %14.6f", state.sma))
    println(io, @sprintf("  ecc        = %14.8f", state.ecc))
    println(io, @sprintf("  inc  (deg) = %14.6f", rad2deg(state.inc)))
    println(io, @sprintf("  raan (deg) = %14.6f", rad2deg(state.raan)))
    println(io, @sprintf("  aop  (deg) = %14.6f", rad2deg(state.aop)))
    println(io, @sprintf("  ta   (deg) = %14.6f", rad2deg(state.ta)))
end

"""
    SphericalRADECState(r, dec, ra, v, decv, rav)

Spherical coordinates state with right ascension and declination components.

# Units
- Distance/time units must be consistent with the gravitational parameter `μ` used in conversions.
- All angular quantities are in **radians**.

# Fields (all `::T` where `T<:Real`)
- `r`: Radial distance. Magnitude of position vector from origin to spacecraft.
  * Range: r > 0. Typically r ≥ 1e-10 for numerical stability.
- `dec`: Declination (rad). Elevation angle of position above/below xy-plane.
  * Range: [-π/2, π/2]. dec = 0: position in xy-plane. dec = ±π/2: at poles.
- `ra`: Right ascension (rad). Azimuthal angle of position in xy-plane.
  * Range: [0, 2π) or (-∞, ∞). Measured counterclockwise from +x axis.
- `v`: Velocity magnitude. Speed of spacecraft motion.
  * Range: v ≥ 0. Typically v ≥ 1e-10 for numerical stability.
- `decv`: Declination of velocity (rad). Elevation angle of velocity vector.
  * Range: [-π/2, π/2]. Angle between velocity vector and xy-plane.
- `rav`: Right ascension of velocity (rad). Azimuthal angle of velocity.
  * Range: [0, 2π) or (-∞, ∞). Direction of velocity in xy-plane.

# Notes
- Parametric for automatic differentiation and arbitrary precision.
- No inherent singularities, but precision loss occurs when r or v approach zero.
- Useful for astronomy applications and telescope pointing.
"""
struct SphericalRADECState{T<:Real} <: AbstractOrbitState
    r::T
    dec::T
    ra::T
    v::T
    decv::T
    rav::T
end

function show(io::IO, state::SphericalRADECState)
    println(io, "SphericalRADECState:")
    println(io, @sprintf("  r         = %14.6f", state.r))
    println(io, @sprintf("  dec (deg) = %11.6f", rad2deg(state.dec)))
    println(io, @sprintf("  ra  (deg) = %11.6f", rad2deg(state.ra)))
    println(io, @sprintf("  v         = %11.6f", state.v))
    println(io, @sprintf("  decv(deg) = %11.6f", rad2deg(state.decv)))
    println(io, @sprintf("  rav (deg) = %11.6f", rad2deg(state.rav)))
end

"""
    SphericalAZIFPAState(r, ra, dec, v, vazi, fpa)

Spherical azimuth / flight-path-angle representation of an orbital state.

# Units
- Distance/time units must be consistent with the gravitational parameter `μ` in use.
- All angular values are radians.

# Fields (all `::T` where `T<:Real`)
- `r`: Radial distance from the central body center
- `ra`: Right ascension
- `dec`: Declination
- `v`: Velocity magnitude
- `vazi`: Velocity azimuth (angle east of north) [rad]
- `fpa`: Flight path angle (above local horizontal) [rad]

# Notes
- Parametric for AD / arbitrary precision. 
"""
struct SphericalAZIFPAState{T<:Real} <: AbstractOrbitState
    r::T
    ra::T
    dec::T
    v::T
    vazi::T
    fpa::T
end

function show(io::IO, state::SphericalAZIFPAState)
    println(io, "SphericalAZIFPAState:")
    println(io, @sprintf("  r          = %14.8f", state.r))
    println(io, @sprintf("  ra   (deg) = %14.8f", rad2deg(state.ra)))
    println(io, @sprintf("  dec  (deg) = %14.8f", rad2deg(state.dec)))
    println(io, @sprintf("  v          = %14.8f", state.v))
    println(io, @sprintf("  vazi (deg) = %14.8f", rad2deg(state.vazi)))
    println(io, @sprintf("  fpa  (deg) = %14.8f", rad2deg(state.fpa)))
end

"""
    ModifiedEquinoctialState(p, f, g, h, k, L)

Modified equinoctial elements representation.

# Units
- Distance/time units must be consistent with the gravitational parameter `μ` used in conversions.
- Angles in radians.

# Fields (all `::T` where `T<:Real`)
- `p`: Semi-latus rectum
- `f`: Eccentricity vector component
- `g`: Eccentricity vector component
- `h`: Inclination vector component
- `k`: Inclination vector component
- `L`: True longitude (rad)

# Notes
- Parametric for AD and high precision.
"""
struct ModifiedEquinoctialState{T<:Real} <: AbstractOrbitState
    p::T
    f::T
    g::T
    h::T
    k::T
    L::T
end

function show(io::IO, state::ModifiedEquinoctialState)
    println(io, "ModifiedEquinoctialState:")
    println(io, @sprintf("  p       = %14.8f", state.p))
    println(io, @sprintf("  f       = %14.8f", state.f))
    println(io, @sprintf("  g       = %14.8f", state.g))
    println(io, @sprintf("  h       = %14.8f", state.h))
    println(io, @sprintf("  k       = %14.8f", state.k))
    println(io, @sprintf("  L (deg) = %14.8f", rad2deg(state.L)))
end

"""
    OutGoingAsymptoteState(rp, c3, rla, dla, bpa, ta)

Outgoing asymptote parameters for hyperbolic trajectories.

# Units
- Distance and time units are arbitrary but must be consistent with the gravitational parameter `μ` used in the simulation.

# Fields (all `::T` where `T<:Real`)
- `rp: Periapsis radius
- `c3: Characteristic energy
- `rla: Right ascension of asymptote (rad)
- `dla: Declination of asymptote (rad)
- `bpa: B-plane angle (rad)
- `ta: True anomaly at asymptote (rad)
"""
struct OutGoingAsymptoteState{T<:Real} <: AbstractOrbitState
    rp::T
    c3::T
    rla::T
    dla::T
    bpa::T
    ta::T
end

function show(io::IO, state::OutGoingAsymptoteState)
    println(io, "OutGoingAsymptoteState:")
    println(io, @sprintf("  rp        = %14.8f", state.rp))
    println(io, @sprintf("  c3        = %14.8f", state.c3))
    println(io, @sprintf("  rla (deg) = %14.8f", rad2deg(state.rla)))
    println(io, @sprintf("  dla (deg) = %14.8f", rad2deg(state.dla)))
    println(io, @sprintf("  bpa (deg) = %14.8f", rad2deg(state.bpa)))
    println(io, @sprintf("  ta  (deg) = %14.8f", rad2deg(state.ta)))
end

"""
    IncomingAsymptoteState(rp, c3, rla, dla, bpa, ta)

Incoming asymptote parameters for hyperbolic trajectories.

# Units
- Distance and time units are arbitrary but must be consistent with the gravitational parameter `μ` used in the simulation.

# Fields (all `::T` where `T<:Real`)
- `rp`: Periapsis radius
- `c3`: Characteristic energy
- `rla`: Right ascension of asymptote (rad)
- `dla`: Declination of asymptote (rad)
- `bpa`: B-plane angle (rad)
- `ta`: True anomaly at asymptote (rad)
"""
struct IncomingAsymptoteState{T<:Real} <: AbstractOrbitState
    rp::T
    c3::T
    rla::T
    dla::T
    bpa::T
    ta::T
end

function show(io::IO, state::IncomingAsymptoteState)
    println(io, "OutGoingAsymptoteState:")
    println(io, @sprintf("  rp        = %14.8f", state.rp))
    println(io, @sprintf("  c3        = %14.8f", state.c3))
    println(io, @sprintf("  rla (deg) = %14.8f", rad2deg(state.rla)))
    println(io, @sprintf("  dla (deg) = %14.8f", rad2deg(state.dla)))
    println(io, @sprintf("  bpa (deg) = %14.8f", rad2deg(state.bpa)))
    println(io, @sprintf("  ta  (deg) = %14.8f", rad2deg(state.ta)))
end

"""
    ModifiedKeplerianState(rp, ra, inc, raan, aop, ta)

Modified Keplerian orbital elements using periapsis and apoapsis radii instead of semi-major axis and eccentricity.

# Units
- Distance and time units must be consistent with the gravitational parameter `μ` used in the simulation.
- All angular quantities are in **radians**.

# Fields (all `::T` where `T<:Real`)
- `rp`: Radius of periapsis. Closest approach distance to central body.
  * Range: rp > 0. Typically rp > central body radius for physical orbits.
- `ra`: Radius of apoapsis. Farthest distance from central body (elliptic orbits only).
  * Range: ra ≥ rp for elliptic orbits. For hyperbolic orbits: ra = ∞ (not used).
- `inc`: Inclination (rad). Angle between orbit plane and reference xy-plane.
  * Range: [0, π]. If inc < π/2: prograde orbit. If inc > π/2: retrograde orbit.
  * inc = 0: equatorial orbit in xy-plane. inc = π/2: polar orbit.
- `raan`: Right ascension of ascending node (rad). Orients the orbit plane.
  * Range: [0, 2π). Angle from +x axis to ascending node, measured in xy-plane.
  * Defines where orbit plane intersects reference plane (ascending crossing).
  * Undefined for equatorial orbits (inc ≈ 0 or π).
- `aop`: Argument of periapsis (rad). Orients the orbit within its plane.
  * Range: [0, 2π). Angle from ascending node to periapsis point.
  * Defines orientation of orbit ellipse within the orbital plane.
  * Undefined for circular orbits (rp ≈ ra).
- `ta`: True anomaly (rad). Spacecraft position within the orbit.
  * Range: [0, 2π). Angle from periapsis to current spacecraft position.
  * ta = 0: at periapsis. ta = π: at apoapsis (for elliptic orbits).
  * Undefined for circular orbits (rp ≈ ra).

# Notes
- Parametric so automatic differentiation and high-precision types are supported.
- Alternative to classical Keplerian elements using radius parameters instead of sma/ecc.
- Shares same singularities as classical Keplerian elements for circular and equatorial orbits.

# Examples
```julia
# 400 km x 35,786 km orbit (GTO-like)
modkep = ModifiedKeplerianState(6778.0, 42164.0, π/6, 0.0, 0.0, 0.0)
```
"""
struct ModifiedKeplerianState{T<:Real} <: AbstractOrbitState
    rp::T
    ra::T
    inc::T
    raan::T
    aop::T
    ta::T
end

function show(io::IO, state::ModifiedKeplerianState)
    println(io, "ModifiedKeplerianState:")
    println(io, @sprintf("  rp         = %14.8f", state.rp))
    println(io, @sprintf("  ra         = %14.8f", state.ra))
    println(io, @sprintf("  inc  (deg) = %14.8f", rad2deg(state.inc)))
    println(io, @sprintf("  raan (deg) = %14.8f", rad2deg(state.raan)))
    println(io, @sprintf("  aop  (deg) = %14.8f", rad2deg(state.aop)))
    println(io, @sprintf("  ta   (deg) = %14.8f", rad2deg(state.ta)))
end

"""
    EquinoctialState(a, h, k, p, q, λ)

Equinoctial orbital elements representation.

# Fields (all `::T` where `T<:Real`)
- `a`: Semi-major axis [must be ≠ 0]
- `h`: e⋅g component of eccentricity vector
- `k`: e⋅f component of eccentricity vector
- `p`: tan(i/2)⋅cos(Ω)
- `q`: tan(i/2)⋅sin(Ω)
- `mlong`: Mean longitude (rad), normalized to [0, 2π)

# Notes
- Distance and time units must be consistent with the gravitational parameter `μ` used elsewhere.
- Angles are always in radians.
"""
struct EquinoctialState{T<:Real} <: AbstractOrbitState
    a::T
    h::T
    k::T
    p::T
    q::T
    mlong::T
end

function show(io::IO, s::EquinoctialState)
    println(io, "EquinoctialState:")
    println(io, @sprintf("  a      = %14.8f", s.a))
    println(io, @sprintf("  h      = %11.8f", s.h))
    println(io, @sprintf("  k      = %11.8f", s.k))
    println(io, @sprintf("  p      = %11.8f", s.p))
    println(io, @sprintf("  q      = %11.8f", s.q))
    println(io, @sprintf("  mlong (deg.) = %11.8f", rad2deg(s.mlong)))
end

"""
    AlternateEquinoctialState(a, h, k, altp, altq, λ)

Alternate equinoctial elements representation.

# Units
- Distance and time units are arbitrary but must be consistent with the gravitational parameter `μ` used in the simulation.
- Angles are in radians.

# Fields (all `::T` where `T<:Real`)
- `a`: Semi-major axis
- `h`: e⋅g component of eccentricity vector
- `k`: e⋅f component of eccentricity vector
- `altp`: sin(i/2)⋅cos(Ω)
- `altq`: sin(i/2)⋅sin(Ω)
- `mlong`: Mean longitude [rad]
"""
struct AlternateEquinoctialState{T<:Real} <: AbstractOrbitState
    a::T
    h::T
    k::T
    altp::T
    altq::T
    mlong::T
end

function show(io::IO, s::AlternateEquinoctialState)
    println(io, "AlternateEquinoctialState:")
    println(io, @sprintf("  a      = %14.8f", s.a))
    println(io, @sprintf("  h      = %11.8f", s.h))
    println(io, @sprintf("  k      = %11.8f", s.k))
    println(io, @sprintf("  altp   = %11.8f", s.altp))
    println(io, @sprintf("  altq   = %11.8f", s.altq))
    println(io, @sprintf("  mlong (deg.)  = %11.8f", rad2deg(s.mlong)))
end

@inline state_tag_to_type(::Cartesian)            = CartesianState            # COV_EXCL_LINE
@inline state_tag_to_type(::Keplerian)            = KeplerianState            # COV_EXCL_LINE
@inline state_tag_to_type(::ModifiedKeplerian)    = ModifiedKeplerianState    # COV_EXCL_LINE
@inline state_tag_to_type(::Equinoctial)          = EquinoctialState          # COV_EXCL_LINE
@inline state_tag_to_type(::AlternateEquinoctial) = AlternateEquinoctialState # COV_EXCL_LINE
@inline state_tag_to_type(::ModifiedEquinoctial)  = ModifiedEquinoctialState  # COV_EXCL_LINE
@inline state_tag_to_type(::SphericalRADEC)       = SphericalRADECState       # COV_EXCL_LINE
@inline state_tag_to_type(::SphericalAZIFPA)      = SphericalAZIFPAState      # COV_EXCL_LINE
@inline state_tag_to_type(::IncomingAsymptote)    = IncomingAsymptoteState    # COV_EXCL_LINE
@inline state_tag_to_type(::OutGoingAsymptote)    = OutGoingAsymptoteState    # COV_EXCL_LINE

@inline state_type_to_tag(::Type{CartesianState})            = Cartesian()
@inline state_type_to_tag(::Type{KeplerianState})            = Keplerian()
@inline state_type_to_tag(::Type{ModifiedKeplerianState})    = ModifiedKeplerian()
@inline state_type_to_tag(::Type{EquinoctialState})          = Equinoctial()
@inline state_type_to_tag(::Type{AlternateEquinoctialState}) = AlternateEquinoctial()
@inline state_type_to_tag(::Type{ModifiedEquinoctialState})  = ModifiedEquinoctial()
@inline state_type_to_tag(::Type{SphericalRADECState})       = SphericalRADEC()
@inline state_type_to_tag(::Type{SphericalAZIFPAState})      = SphericalAZIFPA()
@inline state_type_to_tag(::Type{IncomingAsymptoteState})    = IncomingAsymptote()
@inline state_type_to_tag(::Type{OutGoingAsymptoteState})    = OutGoingAsymptote()
#@inline state_type_to_tag(st::AbstractOrbitState)            = state_type_to_tag(typeof(st))

@inline state_type_to_tag(::Type{<:CartesianState})            = Cartesian()           
@inline state_type_to_tag(::Type{<:KeplerianState})            = Keplerian()           
@inline state_type_to_tag(::Type{<:EquinoctialState})          = Equinoctial()         
@inline state_type_to_tag(::Type{<:SphericalRADECState})       = SphericalRADEC()      
@inline state_type_to_tag(::Type{<:SphericalAZIFPAState})      = SphericalAZIFPA()     
@inline state_type_to_tag(::Type{<:ModifiedEquinoctialState})  = ModifiedEquinoctial() 
@inline state_type_to_tag(::Type{<:OutGoingAsymptoteState})    = OutGoingAsymptote()   
@inline state_type_to_tag(::Type{<:IncomingAsymptoteState})    = IncomingAsymptote()   
@inline state_type_to_tag(::Type{<:ModifiedKeplerianState})    = ModifiedKeplerian()   
@inline state_type_to_tag(::Type{<:AlternateEquinoctialState}) = AlternateEquinoctial()
@inline state_type_to_tag(x::AbstractOrbitState) = state_type_to_tag(typeof(x))        

# Bring in the OrbitState wrapper and helpers now that concrete types exist
include("orbit_state.jl")

# =============================================================================
# Vector Conversion Interface
# =============================================================================

function to_vector(state::CartesianState)
    state.posvel
end
function to_vector(state::KeplerianState)
    [state.sma, state.ecc, state.inc, state.raan, state.aop, state.ta]
end
function to_vector(state::SphericalRADECState)
    [state.r, state.dec, state.ra, state.v, state.decv, state.rav]
end
function to_vector(state::ModifiedEquinoctialState)
    [state.p, state.f, state.g, state.h, state.k, state.L]
end
function to_vector(state::OutGoingAsymptoteState)
    [state.rp, state.c3, state.rla, state.dla, state.bpa, state.ta]
end
function to_vector(state::IncomingAsymptoteState)
    [state.rp, state.c3, state.rla, state.dla, state.bpa, state.ta]
end
function to_vector(state::ModifiedKeplerianState)
    [state.rp, state.ra, state.inc, state.raan, state.aop, state.ta]
end
function to_vector(state::SphericalAZIFPAState)
    [state.r, state.ra, state.dec, state.v, state.vazi, state.fpa]
end
function to_vector(state::EquinoctialState)
    [state.a, state.h, state.k, state.p, state.q, state.mlong]
end
function to_vector(state::AlternateEquinoctialState)
    [state.a, state.h, state.k, state.altp, state.altq, state.mlong]
end
function to_vector(state::OrbitState)
    state.state
end
function to_vector(v::Vector{<:Real})
    error("Cannot call `to_vector` on raw Vector — a state type was expected 
    but got Vector{Float64}. This means one of your conversions returned a 
    vector instead of a struct.")
end

# =============================================================================
# Overloaded Constructors
# =============================================================================

# Identity constructors
CartesianState(state::CartesianState, μ::Real) = state
CartesianState(state::CartesianState) = state
KeplerianState(state::KeplerianState, μ::Real) = state
SphericalRADECState(state::SphericalRADECState, μ::Real) = state
SphericalRADECState(state::SphericalRADECState) = state
ModifiedEquinoctialState(state::ModifiedEquinoctialState, μ::Real) = state
OutGoingAsymptoteState(state::OutGoingAsymptoteState, μ::Real) = state
IncomingAsymptoteState(state::IncomingAsymptoteState, μ::Real) = state
ModifiedKeplerianState(state::ModifiedKeplerianState, μ::Real) = state
SphericalAZIFPAState(state::SphericalAZIFPAState, μ::Real) = state
SphericalAZIFPAState(state::SphericalAZIFPAState) = state
EquinoctialState(state::EquinoctialState, μ::Real) = state
EquinoctialState(state::EquinoctialState) = state
AlternateEquinoctialState(state::AlternateEquinoctialState, μ::Real) = state
AlternateEquinoctialState(state::AlternateEquinoctialState) = state

# Conversions to CartesianState
CartesianState(s::Vector{<:Real}, μ::Real) = CartesianState(s)
CartesianState(keplerian::KeplerianState, μ::Real) = CartesianState(kep_to_cart(to_vector(keplerian), μ))
CartesianState(equinoct::EquinoctialState, μ::Real) = CartesianState(equinoctial_to_cart(to_vector(equinoct),μ))
CartesianState(equinoct::AlternateEquinoctialState, μ::Real) = CartesianState(
    equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ))
CartesianState(modkeplerian::ModifiedKeplerianState, μ::Real) = 
               CartesianState(kep_to_cart(modkep_to_kep(to_vector(modkeplerian)),μ))
CartesianState(sph::SphericalRADECState) = CartesianState(sphradec_to_cart(to_vector(sph)))
CartesianState(sph::SphericalRADECState, μ::Real) = CartesianState(sphradec_to_cart(to_vector(sph)))
CartesianState(sph::SphericalAZIFPAState) = CartesianState(sphazfpa_to_cart(to_vector(sph)))
CartesianState(sph::SphericalAZIFPAState, μ::Real) = CartesianState(sphazfpa_to_cart(to_vector(sph)))
CartesianState(mee::ModifiedEquinoctialState, μ::Real) = CartesianState(mee_to_cart(to_vector(mee), μ))
CartesianState(outasymptote::OutGoingAsymptoteState, μ::Real) = CartesianState(
                   kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ))
CartesianState(inasymptote::IncomingAsymptoteState, μ::Real) = CartesianState(
                   kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ))

# Conversions to KeplerianState
#KeplerianState(s::Vector{<:Real}, μ::Real) = KeplerianState(s...)
#KeplerianState(s::Vector{<:Real}) = KeplerianState(s...)
KeplerianState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    KeplerianState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
KeplerianState(s::AbstractVector{T}) where {T<:Real} = KeplerianState(s, zero(T))

KeplerianState(modkep::ModifiedKeplerianState, μ::Real) = KeplerianState(modkep_to_kep(to_vector(modkep))...)
KeplerianState(cart::CartesianState, μ::Real) = KeplerianState(cart_to_kep(to_vector(cart), μ)...)
KeplerianState(sph::SphericalRADECState, μ::Real) = KeplerianState(CartesianState(sph), μ)
KeplerianState(equinoct::EquinoctialState, μ::Real) = KeplerianState(
                   CartesianState(equinoctial_to_cart(to_vector(equinoct),μ)), μ)
KeplerianState(equinoct::AlternateEquinoctialState, μ::Real) = KeplerianState(
                   CartesianState(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ)), μ)
KeplerianState(sph::SphericalAZIFPAState, μ::Real) = KeplerianState(CartesianState(sph), μ)
KeplerianState(mee::ModifiedEquinoctialState, μ::Real) = KeplerianState(CartesianState(mee, μ), μ)
KeplerianState(outasymptote::OutGoingAsymptoteState, μ::Real) = KeplerianState(
                   outasymptote_to_kep(to_vector(outasymptote), μ)...)
KeplerianState(inasymptote::IncomingAsymptoteState, μ::Real) = KeplerianState(
                   inasymptote_to_kep(to_vector(inasymptote), μ)...)

# Conversions to SphericalRADECState
#SphericalRADECState(s::Vector{<:Real}, μ::Real) = SphericalRADECState(s...)
#SphericalRADECState(s::Vector{<:Real}) = SphericalRADECState(s...)
SphericalRADECState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    SphericalRADECState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
SphericalRADECState(s::AbstractVector{T}) where {T<:Real} = SphericalRADECState(s, zero(T))
SphericalRADECState(cart::CartesianState) = SphericalRADECState(cart_to_sphradec(to_vector(cart))...)
SphericalRADECState(cart::CartesianState,μ::Real) = SphericalRADECState(cart_to_sphradec(to_vector(cart))...)
SphericalRADECState(sph::SphericalAZIFPAState) = SphericalRADECState( cart_to_sphradec(sphazfpa_to_cart(to_vector(sph))))
SphericalRADECState(sph::SphericalAZIFPAState,μ::Real) = SphericalRADECState( cart_to_sphradec(sphazfpa_to_cart(to_vector(sph))))
SphericalRADECState(equinoct::EquinoctialState,μ::Real) = SphericalRADECState(
                                  cart_to_sphradec(equinoctial_to_cart(to_vector(equinoct),μ)))
SphericalRADECState(equinoct::AlternateEquinoctialState,μ::Real) = SphericalRADECState(
                                  cart_to_sphradec(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ)))
SphericalRADECState(keplerian::KeplerianState, μ::Real) = SphericalRADECState(CartesianState(keplerian, μ))
SphericalRADECState(modkep::ModifiedKeplerianState, μ::Real) = SphericalRADECState(
    CartesianState(kep_to_cart(modkep_to_kep(to_vector(modkep)), μ)))
SphericalRADECState(mee::ModifiedEquinoctialState, μ::Real) = SphericalRADECState(CartesianState(mee, μ))
SphericalRADECState(outasymptote::OutGoingAsymptoteState, μ::Real) = SphericalRADECState(
       cart_to_sphradec(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ))...)
SphericalRADECState(inasymptote::IncomingAsymptoteState, μ::Real) = SphericalRADECState(
       cart_to_sphradec(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ))...)

# Conversions to SphericalAZIFPAState
#SphericalAZIFPAState(s::Vector{<:Real}, μ::Real) = SphericalAZIFPAState(s...)
#SphericalAZIFPAState(s::Vector{<:Real}) = SphericalAZIFPAState(s...)
SphericalAZIFPAState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    SphericalAZIFPAState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
SphericalAZIFPAState(s::AbstractVector{T}) where {T<:Real} = SphericalAZIFPAState(s, zero(T))
SphericalAZIFPAState(cart::CartesianState) = SphericalAZIFPAState(cart_to_sphazfpa(to_vector(cart))...)
SphericalAZIFPAState(cart::CartesianState,μ::Real) = SphericalAZIFPAState(cart_to_sphazfpa(to_vector(cart))...)
SphericalAZIFPAState(sph::SphericalRADECState) = SphericalAZIFPAState( cart_to_sphazfpa(sphradec_to_cart(to_vector(sph))))
SphericalAZIFPAState(sph::SphericalRADECState,μ::Real) = SphericalAZIFPAState( cart_to_sphazfpa(sphradec_to_cart(to_vector(sph))))
SphericalAZIFPAState(equinoct::EquinoctialState,μ::Real) = SphericalAZIFPAState( 
                                cart_to_sphazfpa(equinoctial_to_cart(to_vector(equinoct),μ)))
SphericalAZIFPAState(equinoct::AlternateEquinoctialState,μ::Real) = SphericalAZIFPAState( 
                                cart_to_sphazfpa(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ)))
SphericalAZIFPAState(keplerian::KeplerianState, μ::Real) = SphericalAZIFPAState(CartesianState(keplerian, μ))
SphericalAZIFPAState(modkep::ModifiedKeplerianState, μ::Real) = SphericalAZIFPAState(
    CartesianState(kep_to_cart(modkep_to_kep(to_vector(modkep)), μ)))
SphericalAZIFPAState(mee::ModifiedEquinoctialState, μ::Real) = SphericalAZIFPAState(CartesianState(mee, μ))
SphericalAZIFPAState(outasymptote::OutGoingAsymptoteState, μ::Real) = SphericalAZIFPAState(
       cart_to_sphazfpa(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ))...)
SphericalAZIFPAState(inasymptote::IncomingAsymptoteState, μ::Real) = SphericalAZIFPAState(
       cart_to_sphazfpa(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ))...)

# Conversions to ModifiedEquinoctialState
#ModifiedEquinoctialState(s::Vector{<:Real}, μ::Real) = ModifiedEquinoctialState(s...)
#ModifiedEquinoctialState(s::Vector{<:Real}) = ModifiedEquinoctialState(s...)
ModifiedEquinoctialState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    ModifiedEquinoctialState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
ModifiedEquinoctialState(s::AbstractVector{T}) where {T<:Real} =
    ModifiedEquinoctialState(s, zero(T))
ModifiedEquinoctialState(cart::CartesianState, μ::Real) = ModifiedEquinoctialState(cart_to_mee(to_vector(cart), μ)...)
ModifiedEquinoctialState(equinoct::EquinoctialState, μ::Real) = ModifiedEquinoctialState(
                           cart_to_mee(equinoctial_to_cart(to_vector(equinoct),μ), μ)...)
ModifiedEquinoctialState(equinoct::AlternateEquinoctialState, μ::Real) = ModifiedEquinoctialState(
                           cart_to_mee(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ), μ)...)
ModifiedEquinoctialState(keplerian::KeplerianState, μ::Real) = ModifiedEquinoctialState(CartesianState(keplerian, μ), μ)
ModifiedEquinoctialState(modkep::ModifiedKeplerianState, μ::Real) = ModifiedEquinoctialState(
    CartesianState( kep_to_cart(modkep_to_kep(to_vector(modkep)), μ)) ,μ)
ModifiedEquinoctialState(sph::SphericalRADECState, μ::Real) = ModifiedEquinoctialState(CartesianState(sph, μ), μ)
ModifiedEquinoctialState(sph::SphericalAZIFPAState, μ::Real) = ModifiedEquinoctialState(CartesianState(sph, μ), μ)
ModifiedEquinoctialState(outasymptote::OutGoingAsymptoteState, μ::Real) = ModifiedEquinoctialState(
       cart_to_mee(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ),μ)...)
ModifiedEquinoctialState(inasymptote::IncomingAsymptoteState, μ::Real) = ModifiedEquinoctialState(
       cart_to_mee(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ),μ)...)

# Conversions to OutGoingAsymptoteState
#OutGoingAsymptoteState(s::Vector{<:Real}, μ::Real) = OutGoingAsymptoteState(s...)
#OutGoingAsymptoteState(s::Vector{<:Real}) = OutGoingAsymptoteState(s...)
OutGoingAsymptoteState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    OutGoingAsymptoteState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
OutGoingAsymptoteState(s::AbstractVector{T}) where {T<:Real} = OutGoingAsymptoteState(s, zero(T))
OutGoingAsymptoteState(cart::CartesianState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(to_vector(cart),μ)...)
OutGoingAsymptoteState(equinoct::EquinoctialState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(equinoctial_to_cart(to_vector(equinoct),μ),μ)...)
OutGoingAsymptoteState(equinoct::AlternateEquinoctialState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ),μ)...)
OutGoingAsymptoteState(kep::KeplerianState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(kep_to_cart(to_vector(kep),μ),μ)...)
OutGoingAsymptoteState(modkep::ModifiedKeplerianState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(kep_to_cart(modkep_to_kep(to_vector(modkep)),μ),μ)...)
OutGoingAsymptoteState(sph::SphericalRADECState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(sphradec_to_cart(to_vector(sph)),μ)...)
OutGoingAsymptoteState(sph::SphericalAZIFPAState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(sphazfpa_to_cart(to_vector(sph)),μ)...)
OutGoingAsymptoteState(mee::ModifiedEquinoctialState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(mee_to_cart(to_vector(mee),μ),μ)...)
OutGoingAsymptoteState(inasymptote::IncomingAsymptoteState, μ::Real) = OutGoingAsymptoteState(
       cart_to_outasymptote(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote),μ),μ),μ)...)

# Conversions to IncomingAsymptoteState
#IncomingAsymptoteState(s::Vector{<:Real}, μ::Real) = IncomingAsymptoteState(s...)
#IncomingAsymptoteState(s::Vector{<:Real}) = IncomingAsymptoteState(s...)
IncomingAsymptoteState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    IncomingAsymptoteState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
IncomingAsymptoteState(s::AbstractVector{T}) where {T<:Real} = IncomingAsymptoteState(s, zero(T))

IncomingAsymptoteState(cart::CartesianState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(to_vector(cart),μ)...)
IncomingAsymptoteState(equinoct::EquinoctialState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(equinoctial_to_cart(to_vector(equinoct),μ),μ)...)
IncomingAsymptoteState(equinoct::AlternateEquinoctialState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ),μ)...)
IncomingAsymptoteState(kep::KeplerianState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(kep_to_cart(to_vector(kep),μ),μ)...)
IncomingAsymptoteState(modkep::ModifiedKeplerianState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(kep_to_cart(modkep_to_kep(to_vector(modkep)),μ),μ)...)
IncomingAsymptoteState(sph::SphericalRADECState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(sphradec_to_cart(to_vector(sph)),μ)...)
IncomingAsymptoteState(sph::SphericalAZIFPAState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(sphazfpa_to_cart(to_vector(sph)),μ)...)
IncomingAsymptoteState(mee::ModifiedEquinoctialState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(mee_to_cart(to_vector(mee),μ),μ)...)
IncomingAsymptoteState(outasymptote::OutGoingAsymptoteState, μ::Real) = IncomingAsymptoteState(
       cart_to_inasymptote(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote),μ),μ),μ)...)

# Conversions to Modified Keplerian State
#ModifiedKeplerianState(s::Vector{<:Real}, μ::Real) = ModifiedKeplerianState(s...)
#ModifiedKeplerianState(s::Vector{<:Real}) = ModifiedKeplerianState(s...)
ModifiedKeplerianState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    ModifiedKeplerianState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
ModifiedKeplerianState(s::AbstractVector{T}) where {T<:Real} = ModifiedKeplerianState(s, zero(T))
ModifiedKeplerianState(kep::KeplerianState, μ::Real) = ModifiedKeplerianState(kep_to_modkep(to_vector(kep))...)
ModifiedKeplerianState(cart::CartesianState, μ::Real) = ModifiedKeplerianState(kep_to_modkep(cart_to_kep(to_vector(cart), μ))...)
ModifiedKeplerianState(equinoct::EquinoctialState, μ::Real) = ModifiedKeplerianState(
                       kep_to_modkep(cart_to_kep(equinoctial_to_cart(to_vector(equinoct), μ), μ))...)
ModifiedKeplerianState(equinoct::AlternateEquinoctialState, μ::Real) = ModifiedKeplerianState(
                       kep_to_modkep(cart_to_kep(equinoctial_to_cart(alt_equinoctial_to_equinoctial(to_vector(equinoct)), μ), μ))...)
ModifiedKeplerianState(sph::SphericalRADECState, μ::Real) = ModifiedKeplerianState(CartesianState(sph), μ)
ModifiedKeplerianState(sph::SphericalAZIFPAState, μ::Real) = ModifiedKeplerianState(CartesianState(sph), μ)
ModifiedKeplerianState(mee::ModifiedEquinoctialState, μ::Real) = ModifiedKeplerianState(CartesianState(mee, μ), μ)
ModifiedKeplerianState(outasymptote::OutGoingAsymptoteState, μ::Real) = ModifiedKeplerianState(
                   kep_to_modkep(outasymptote_to_kep(to_vector(outasymptote), μ))...)
ModifiedKeplerianState(inasymptote::IncomingAsymptoteState, μ::Real) = ModifiedKeplerianState(
                   kep_to_modkep(inasymptote_to_kep(to_vector(inasymptote), μ))...)

# Conversions to EquinoctialState
#EquinoctialState(s::Vector{<:Real}, μ::Real) = EquinoctialState(s...)
#EquinoctialState(s::Vector{<:Real}) = EquinoctialState(s...)
EquinoctialState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    EquinoctialState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
EquinoctialState(s::AbstractVector{T}) where {T<:Real} =
    EquinoctialState(s, zero(T))
EquinoctialState(equinoct::AlternateEquinoctialState,μ::Real) = EquinoctialState(alt_equinoctial_to_equinoctial(to_vector(equinoct)),μ)
EquinoctialState(cart::CartesianState, μ::Real) = EquinoctialState(cart_to_equinoctial(to_vector(cart), μ)...)
EquinoctialState(mee::ModifiedEquinoctialState, μ::Real) = EquinoctialState(
                           cart_to_equinoctial(mee_to_cart(to_vector(mee),μ), μ)...)
EquinoctialState(keplerian::KeplerianState, μ::Real) = EquinoctialState(CartesianState(keplerian, μ), μ)
EquinoctialState(modkep::ModifiedKeplerianState, μ::Real) = EquinoctialState(
    CartesianState( kep_to_cart(modkep_to_kep(to_vector(modkep)), μ)) ,μ)
EquinoctialState(sph::SphericalRADECState, μ::Real) = EquinoctialState(CartesianState(sph, μ), μ)
EquinoctialState(sph::SphericalAZIFPAState, μ::Real) = EquinoctialState(CartesianState(sph, μ), μ)
EquinoctialState(outasymptote::OutGoingAsymptoteState, μ::Real) = EquinoctialState(
       cart_to_equinoctial(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ),μ)...)
EquinoctialState(inasymptote::IncomingAsymptoteState, μ::Real) = EquinoctialState(
       cart_to_equinoctial(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ),μ)...)

# Conversions to AlternateEquinoctialState
#AlternateEquinoctialState(s::Vector{<:Real}, μ::Real) = AlternateEquinoctialState(s...)
#AlternateEquinoctialState(s::Vector{<:Real}) = AlternateEquinoctialState(s...)
AlternateEquinoctialState(s::AbstractVector{T}, μ::Real) where {T<:Real} = begin
    @assert length(s) == 6
    AlternateEquinoctialState{T}(s[1], s[2], s[3], s[4], s[5], s[6])
end
AlternateEquinoctialState(s::AbstractVector{T}) where {T<:Real} =
    AlternateEquinoctialState(s, zero(T))
AlternateEquinoctialState(equinoct::EquinoctialState,μ::Real) = AlternateEquinoctialState(equinoctial_to_alt_equinoctial(to_vector(equinoct)),μ)
AlternateEquinoctialState(cart::CartesianState, μ::Real) = AlternateEquinoctialState(equinoctial_to_alt_equinoctial(cart_to_equinoctial(to_vector(cart), μ))...)
AlternateEquinoctialState(mee::ModifiedEquinoctialState, μ::Real) = AlternateEquinoctialState(
                           equinoctial_to_alt_equinoctial(cart_to_equinoctial(mee_to_cart(to_vector(mee),μ), μ))...)
AlternateEquinoctialState(kep::KeplerianState, μ::Real) = AlternateEquinoctialState(
    CartesianState( kep_to_cart(to_vector(kep), μ)) ,μ)
AlternateEquinoctialState(modkep::ModifiedKeplerianState, μ::Real) = AlternateEquinoctialState(
    CartesianState( kep_to_cart(modkep_to_kep(to_vector(modkep)), μ)) ,μ)
AlternateEquinoctialState(sph::SphericalRADECState, μ::Real) = AlternateEquinoctialState(CartesianState(sph, μ), μ)
AlternateEquinoctialState(sph::SphericalAZIFPAState, μ::Real) = AlternateEquinoctialState(CartesianState(sph, μ), μ)
AlternateEquinoctialState(outasymptote::OutGoingAsymptoteState, μ::Real) = AlternateEquinoctialState(
       equinoctial_to_alt_equinoctial(cart_to_equinoctial(kep_to_cart(outasymptote_to_kep(to_vector(outasymptote), μ),μ),μ))...)
AlternateEquinoctialState(inasymptote::IncomingAsymptoteState, μ::Real) = AlternateEquinoctialState(
       equinoctial_to_alt_equinoctial(cart_to_equinoctial(kep_to_cart(inasymptote_to_kep(to_vector(inasymptote), μ),μ),μ))...)

# =============================================================================
# Fallbacks for unsupported conversions
# =============================================================================

KeplerianState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to KeplerianState")
SphericalRADECState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to SphericalRADECState")
ModifiedEquinoctialState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to ModifiedEquinoctialState")
CartesianState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to CartesianState")
SphericalAZIFPAState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to SphericalAZIFPAState")
ModifiedKeplerianState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to ModifiedKeplerianState")
EquinoctialState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to EquinoctialState")
AlternateEquinoctialState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to AlternateEquinoctialState")
OutGoingAsymptoteState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to OutGoingAsymptoteState")
IncomingAsymptoteState(state::AbstractOrbitState, μ::Real) = error("No conversion defined from $(typeof(state)) to IncomingAsymptoteState")

end
