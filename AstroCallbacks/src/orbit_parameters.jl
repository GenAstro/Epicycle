# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbit parameters of the current state (#11, sections 1, 2 and 6).
#
# The mathematics lives below this package. A value a state representation already carries is
# read from that representation: the apsis radii from ModifiedKeplerianState, the semi-latus
# rectum and true longitude from ModifiedEquinoctialState, and the B-plane from
# IncomingAsymptoteState. Two-body formulas of the state and the anomaly conversions are
# AstroRoutines functions; the two-body ones are the same functions, given a method here that
# takes a spacecraft. What stays here is the quantity interface: the frame, the gravitational
# parameter of the frame's origin, the domain errors, labels, units and setters.
#
# Each reads the subject in `cs`, or in its own frame without one. A parameter undefined for the
# orbit throws an ArgumentError saying why, rather than returning NaN, because a NaN in a
# constraint or a stopping condition fails somewhere far from its cause.
#
# Every reader is generic in the state's element type, so ForwardDiff passes through.
# =============================================================================

using AstroStates: KeplerianState, ModifiedKeplerianState, ModifiedEquinoctialState,
                   IncomingAsymptoteState, CartesianState, to_vector
using AstroRoutines: true_to_mean_anomaly, true_to_eccentric_anomaly, mean_to_true_anomaly,
                     true_to_hyperbolic_anomaly, hyperbolic_to_mean_anomaly,
                     mean_to_hyperbolic_anomaly, hyperbolic_to_true_anomaly
import AstroRoutines: specific_energy, hyperbolic_excess_speed, apoapsis_speed, periapsis_speed,
                      angular_momentum, angular_momentum_magnitude, orbit_normal,
                      orbital_angular_velocity, eccentricity_vector, orbital_period, mean_motion

# --- shared reads ------------------------------------------------------------

"Position, velocity and gravitational parameter of `subject` in `cs`."
function _rvμ(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    x = to_vector(CartesianState(Coordinate(subject, cs, params)))
    return x[1:3], x[4:6], get_gravparam(cs.origin)
end

"The state of `subject` in `cs` as representation `S`, about the origin's gravitational parameter."
function _state(S, subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    return S(CartesianState(Coordinate(subject, cs, params)), get_gravparam(cs.origin))
end

"Keplerian elements and gravitational parameter of `subject` in `cs`."
function _keplerian(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    return _state(KeplerianState, subject, cs, params), get_gravparam(cs.origin)
end

_wrap2π(x) = mod(x, 2π)

function _require_closed(e, what)
    e < 1 || throw(ArgumentError(
        "$what is defined only for a closed orbit (eccentricity below 1); this orbit has " *
        "eccentricity $e."))
    return nothing
end

function _require_hyperbolic(e, what)
    e > 1 || throw(ArgumentError(
        "$what is defined only for a hyperbolic orbit (eccentricity above 1); this orbit has " *
        "eccentricity $e."))
    return nothing
end

# --- anomalies ---------------------------------------------------------------

"""
    mean_anomaly(subject[, cs][, params])

Osculating mean anomaly of `subject`, in radians, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
- On a closed orbit, M = E − e sin E, in [0, 2π).
- On a hyperbola, M = e sinh H − H, negative before periapsis and positive after. It is not an
  angle and does not wrap.
- A parabola (eccentricity exactly 1) has no mean anomaly in either form and throws.
- Settable with [`mean_anomaly!`](@ref).

# Returns
The mean anomaly in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
mean_anomaly(sat)
mean_anomaly(sat, EarthMJ2000Ec)
```
"""
function mean_anomaly(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    e, ν = kep.ecc, kep.ta
    e < 1 && return _wrap2π(true_to_mean_anomaly(ν, e))
    e > 1 && return hyperbolic_to_mean_anomaly(true_to_hyperbolic_anomaly(ν, e), e)
    throw(ArgumentError("a parabolic orbit (eccentricity 1) has no mean anomaly."))
end

"""
    eccentric_anomaly(subject[, cs][, params])

Eccentric anomaly of `subject`, in radians in [0, 2π), measured in `cs` or the subject's own
frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only for a closed orbit; a hyperbola throws. Its counterpart there is
[`hyperbolic_anomaly`](@ref).

# Returns
The eccentric anomaly in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
eccentric_anomaly(sat)
```
"""
function eccentric_anomaly(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    _require_closed(kep.ecc, "Eccentric anomaly")
    return _wrap2π(true_to_eccentric_anomaly(kep.ta, kep.ecc))
end

"""
    hyperbolic_anomaly(subject[, cs][, params])

Hyperbolic anomaly H of `subject`, in radians, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only for a hyperbolic orbit; a closed orbit throws. H is negative before periapsis and
positive after, with tanh(H/2) = √((e−1)/(e+1)) tan(ν/2).

# Returns
The hyperbolic anomaly in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
hyperbolic_anomaly(sat)
```
"""
function hyperbolic_anomaly(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    _require_hyperbolic(kep.ecc, "Hyperbolic anomaly")
    return true_to_hyperbolic_anomaly(kep.ta, kep.ecc)
end

# --- energy and speed --------------------------------------------------------

"""
    specific_energy(subject[, cs][, params])

Specific orbital energy of `subject`, v²/2 − μ/r, in km²/s², measured in `cs` or the subject's
own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Negative on a closed orbit, positive on a hyperbola. In a rotating frame the velocity is the
rotating-frame velocity, so the value is not the inertial energy.

# Returns
The specific energy in km²/s².

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
specific_energy(sat)
```
"""
specific_energy(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    specific_energy(_rvμ(subject, cs, params)...)

"""
    hyperbolic_excess_speed(subject[, cs][, params])

Hyperbolic excess speed v∞ of `subject`, √(2ε) = √(−μ/a), in km/s, measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only on an escape orbit (positive specific energy); a closed orbit throws. Its square is
the characteristic energy C3.

# Returns
The hyperbolic excess speed in km/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
hyperbolic_excess_speed(sat)
```
"""
hyperbolic_excess_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    hyperbolic_excess_speed(_rvμ(subject, cs, params)...)

"""
    apoapsis_speed(subject[, cs][, params])

Speed of `subject` at apoapsis, μ(1 − e)/h, in km/s, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only for a closed orbit; a hyperbola has no apoapsis and throws. The value is the
two-body speed the current osculating orbit has at its apoapsis, not a propagated one.

# Returns
The apoapsis speed in km/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
apoapsis_speed(sat)
```
"""
apoapsis_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    apoapsis_speed(_rvμ(subject, cs, params)...)

"""
    periapsis_speed(subject[, cs][, params])

Speed of `subject` at periapsis, μ(1 + e)/h, in km/s, measured in `cs` or the subject's own
frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined for every conic. The value is the two-body speed the current osculating orbit has at its
periapsis, not a propagated one.

# Returns
The periapsis speed in km/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
periapsis_speed(sat)
```
"""
periapsis_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    periapsis_speed(_rvμ(subject, cs, params)...)

# --- angular momentum and the orbit plane ------------------------------------

"""
    angular_momentum(subject[, cs][, params])

Specific angular momentum vector of `subject`, r × v, in km²/s, expressed in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Returns
A three-element vector in km²/s, in the axes of `cs`.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
angular_momentum(sat)
angular_momentum(sat, EarthMJ2000Ec)
```
"""
function angular_momentum(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, _ = _rvμ(subject, cs, params)
    return angular_momentum(r, v)
end

"""
    angular_momentum_magnitude(subject[, cs][, params])

Magnitude of the specific angular momentum of `subject`, ‖r × v‖, in km²/s, measured in `cs` or
the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Returns
The angular momentum magnitude in km²/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
angular_momentum_magnitude(sat)
```
"""
function angular_momentum_magnitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, _ = _rvμ(subject, cs, params)
    return angular_momentum_magnitude(r, v)
end

"""
    orbit_normal(subject[, cs][, params])

Unit vector along the angular momentum of `subject`, (r × v)/‖r × v‖, expressed in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Undefined for rectilinear motion, where r × v is zero, and throws there.

# Returns
A three-element unit vector, in the axes of `cs`.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
orbit_normal(sat)
```
"""
function orbit_normal(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, _ = _rvμ(subject, cs, params)
    return orbit_normal(r, v)
end

"""
    orbital_angular_velocity(subject[, cs][, params])

Orbital angular velocity vector of `subject`, (r × v)/r², in rad/s, expressed in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
The instantaneous rate at which the position vector turns, along the orbit normal. Its
magnitude equals the mean motion only on a circular orbit.

# Returns
A three-element vector in rad/s, in the axes of `cs`.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
orbital_angular_velocity(sat)
```
"""
function orbital_angular_velocity(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, _ = _rvμ(subject, cs, params)
    return orbital_angular_velocity(r, v)
end

"""
    eccentricity_vector(subject[, cs][, params])

Eccentricity vector of `subject`, pointing to periapsis with magnitude e, expressed in `cs` or
the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
e = ((v² − μ/r) r − (r·v) v)/μ. On a circular orbit it is the zero vector and gives no
periapsis direction.

# Returns
A three-element dimensionless vector, in the axes of `cs`.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
eccentricity_vector(sat)
```
"""
eccentricity_vector(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    eccentricity_vector(_rvμ(subject, cs, params)...)

# --- period and motion -------------------------------------------------------

"""
    orbital_period(subject[, cs][, params])

Orbital period of `subject`, 2π√(a³/μ), in seconds, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only for a closed orbit (negative specific energy); an escape orbit throws. It is the
period of the osculating two-body orbit, not a measured one.

# Returns
The period in seconds.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
orbital_period(sat)
```
"""
orbital_period(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    orbital_period(_rvμ(subject, cs, params)...)

"""
    mean_motion(subject[, cs][, params])

Mean motion of `subject`, √(μ/|a|³), in rad/s, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined on closed orbits and hyperbolas, where it is the rate of the hyperbolic mean anomaly.
A parabola has no semi-major axis and throws.

# Returns
The mean motion in rad/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
mean_motion(sat)
```
"""
mean_motion(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    mean_motion(_rvμ(subject, cs, params)...)

# --- angles in the orbit plane -----------------------------------------------

"""
    argument_of_latitude(subject[, cs][, params])

Argument of latitude of `subject`, ω + ν, in radians in [0, 2π), measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
The angle from the ascending node to the position, in the orbit plane. It stays defined on a
circular orbit, where ω and ν separately do not.

# Returns
The argument of latitude in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
argument_of_latitude(sat)
```
"""
function argument_of_latitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    return _wrap2π(kep.aop + kep.ta)
end

"""
    longitude_of_periapsis(subject[, cs][, params])

Longitude of periapsis of `subject`, Ω + ω, in radians in [0, 2π), measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Stays defined on an equatorial orbit, where Ω and ω separately do not.

# Returns
The longitude of periapsis in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
longitude_of_periapsis(sat)
```
"""
function longitude_of_periapsis(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    return _wrap2π(kep.raan + kep.aop)
end

"""
    descending_node_right_ascension(subject[, cs][, params])

Right ascension of the descending node of `subject`, Ω + π, in radians in [0, 2π), measured in
`cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Returns
The right ascension of the descending node in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
descending_node_right_ascension(sat)
```
"""
function descending_node_right_ascension(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    return _wrap2π(kep.raan + π)
end

# --- named shortcuts (#11, section 2) ----------------------------------------
#
# Each is a field of a state representation, read from it; a name makes it usable directly in a
# constraint or a report.

"""
    apoapsis_radius(subject[, cs][, params])

Apoapsis radius of `subject`, a(1 + e), in km, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined only for a closed orbit; a hyperbola throws. Settable with
[`apoapsis_radius!`](@ref).

# Returns
The apoapsis radius in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
apoapsis_radius(sat)
```
"""
function apoapsis_radius(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    radii = _state(ModifiedKeplerianState, subject, cs, params)
    _require_closed((radii.ra - radii.rp) / (radii.ra + radii.rp), "Apoapsis radius")
    return radii.ra
end

"""
    periapsis_radius(subject[, cs][, params])

Periapsis radius of `subject`, p/(1 + e), in km, measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Defined for every conic. Settable with [`periapsis_radius!`](@ref).

# Returns
The periapsis radius in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
periapsis_radius(sat)
```
"""
periapsis_radius(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _state(ModifiedKeplerianState, subject, cs, params).rp

"""
    semilatus_rectum(subject[, cs][, params])

Semi-latus rectum of `subject`, p = h²/μ = a(1 − e²), in km, measured in `cs` or the subject's
own frame.

# Arguments
See [`position_vector`](@ref).

# Returns
The semi-latus rectum in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
semilatus_rectum(sat)
```
"""
semilatus_rectum(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _state(ModifiedEquinoctialState, subject, cs, params).p

"""
    c3(subject[, cs][, params])

Characteristic energy of `subject`, v² − 2μ/r = −μ/a, in km²/s², measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Twice the specific energy. Positive on an escape orbit, where it is v∞², and negative on a
closed orbit.

# Returns
C3 in km²/s².

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
c3(sat)
```
"""
# The asymptote states carry C3 too, but their conversion fails on a circular orbit, where C3 is
# still defined; twice the specific energy is defined everywhere.
c3(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    2 * specific_energy(_rvμ(subject, cs, params)...)

"""
    true_longitude(subject[, cs][, params])

True longitude of `subject`, Ω + ω + ν, in radians in [0, 2π), measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Stays defined on circular and equatorial orbits, where its three terms separately do not.

# Returns
The true longitude in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
true_longitude(sat)
```
"""
true_longitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _wrap2π(_state(ModifiedEquinoctialState, subject, cs, params).L)

for (fn, field, what, other) in
        ((:apoapsis_radius,  :ra, "apoapsis radius",  "periapsis radius"),
         (:periapsis_radius, :rp, "periapsis radius", "apoapsis radius"))
    fn! = Symbol(fn, "!")
    doc = """
        $(fn!)(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to` km, measured in `cs` or the spacecraft's own frame. The
    $other, inclination, RAAN, argument of periapsis and true anomaly are held in that frame.

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when the radii that result
    do not describe an orbit.

    # Examples
    ```julia
    sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
    $(fn!)(sat; to = $(fn === :apoapsis_radius ? "7500.0" : "6800.0"))
    ```
    """
    @eval begin
        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::Real)
            isfinite(to) && to > 0 || throw(ArgumentError(
                "$($what) must be a positive number of km; got $to. The state was not changed."))
            return set_element!(sat, ModifiedKeplerianState, cs, $(QuoteNode(field)), params; to = to)
        end
        @doc $doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)
    end
end

# --- B-plane (#11, section 6) --------------------------------------------------
#
# The B-plane is normal to the incoming asymptote Ŝ and passes through the body's centre. T̂ lies
# along Ŝ × Ẑ, with Ẑ the coordinate system's Z axis, and R̂ = Ŝ × T̂. B points from the centre to
# where the incoming asymptote pierces the plane; its length is the impact parameter h/v∞.
#
# IncomingAsymptoteState carries the B-plane: its periapsis radius and C3 give |B|, and its
# B-vector azimuth, measured from R̂ toward T̂ (GMAT's IncomingBVAZI, counter-clockwise from
# south), gives the direction, so B·T = |B| sin(azimuth) and B·R = |B| cos(azimuth).

"|B|, B·T and B·R of a hyperbolic state, from its incoming asymptote, in km."
function _b_plane(subject, cs::AbstractCoordinateSystem, params::NamedTuple, what)
    asymptote = _state(IncomingAsymptoteState, subject, cs, params)
    μ = get_gravparam(cs.origin)
    any(isnan, to_vector(asymptote)) && throw(ArgumentError(
        "$what is undefined for this state: its incoming asymptote could not be formed (the " *
        "orbit is circular or parabolic, or the asymptote lies along the Z axis of the " *
        "coordinate system)."))
    rp, c3, azimuth = asymptote.rp, asymptote.c3, asymptote.bpa
    _require_hyperbolic(1 + rp * c3 / μ, what)
    magnitude = sqrt(rp^2 + 2μ * rp / c3)
    return magnitude, magnitude * sin(azimuth), magnitude * cos(azimuth)
end

"""
    b_dot_t(subject[, cs][, params])

Component of the B vector of `subject` along T̂, in km, measured in `cs` or the subject's own
frame.

# Arguments
See [`position_vector`](@ref).

# Notes
T̂ lies along Ŝ × Ẑ, with Ŝ the incoming asymptote and Ẑ the Z axis of `cs`, so B·T depends on the
coordinate system: an equatorial and an ecliptic frame give different values. Defined only for a
hyperbolic orbit; a closed orbit throws.

# Returns
B·T in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
b_dot_t(sat)
b_dot_t(sat, EarthMJ2000Ec)
```
"""
b_dot_t(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _b_plane(subject, cs, params, "B·T")[2]

"""
    b_dot_r(subject[, cs][, params])

Component of the B vector of `subject` along R̂ = Ŝ × T̂, in km, measured in `cs` or the subject's
own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Depends on the coordinate system, as [`b_dot_t`](@ref) does. Defined only for a hyperbolic orbit.

# Returns
B·R in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
b_dot_r(sat)
```
"""
b_dot_r(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _b_plane(subject, cs, params, "B·R")[3]

"""
    b_magnitude(subject[, cs][, params])

Length of the B vector of `subject`, the impact parameter h/v∞, in km, measured in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
The miss distance of the incoming asymptote from the body's centre. Unlike its components it
does not depend on the coordinate system's orientation. Defined only for a hyperbolic orbit.

# Returns
The B-vector magnitude in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
b_magnitude(sat)
```
"""
b_magnitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    _b_plane(subject, cs, params, "B-vector magnitude")[1]

"""
    b_plane_angle(subject[, cs][, params])

Angle of the B vector of `subject` from T̂ toward R̂, atan(B·R, B·T), in radians in [0, 2π),
measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
- GMAT's `BVectorAngle`. It is 90° less the B-vector azimuth that `IncomingAsymptoteState` carries
  (GMAT's `IncomingBVAZI`), which is measured from R̂ toward T̂.
- Depends on the coordinate system, as [`b_dot_t`](@ref) does. Defined only for a hyperbolic
  orbit.

# Returns
The B-plane angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
b_plane_angle(sat)
```
"""
function b_plane_angle(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    _, bt, br = _b_plane(subject, cs, params, "B-plane angle")
    return _wrap2π(atan(br, bt))
end

# --- the shorter methods -----------------------------------------------------

const _ORBIT_PARAMETERS = (
    :mean_anomaly, :eccentric_anomaly, :hyperbolic_anomaly,
    :specific_energy, :hyperbolic_excess_speed, :apoapsis_speed, :periapsis_speed,
    :angular_momentum, :angular_momentum_magnitude, :orbit_normal, :orbital_angular_velocity,
    :eccentricity_vector, :orbital_period, :mean_motion,
    :argument_of_latitude, :longitude_of_periapsis, :descending_node_right_ascension,
    :apoapsis_radius, :periapsis_radius, :semilatus_rectum, :c3, :true_longitude,
    :b_dot_t, :b_dot_r, :b_magnitude, :b_plane_angle)

for fn in _ORBIT_PARAMETERS
    @eval begin
        $fn(subject, cs::AbstractCoordinateSystem) = $fn(subject, cs, NamedTuple())
        $fn(subject) = $fn(subject, frame_of(subject))
    end
end

# --- setting the mean anomaly ------------------------------------------------

"True anomaly for mean anomaly `M` and eccentricity `e`, from Kepler's equation."
_true_from_mean(M, e) = e < 1 ? _wrap2π(mean_to_true_anomaly(rem2pi(M, RoundNearest), e)) :
                        hyperbolic_to_true_anomaly(mean_to_hyperbolic_anomaly(M, e), e)

"""
    mean_anomaly!(sat[, cs][, params]; to) -> sat

Set the mean anomaly of `sat` to `to` radians, measured in `cs` or the spacecraft's own frame.
The other Keplerian elements are held in that frame.

# Notes
On a closed orbit any value is accepted and taken modulo 2π. On a hyperbola, `to` is the
hyperbolic mean anomaly e sinh H − H and is not wrapped.

# Returns
`sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when `to` is not finite.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
mean_anomaly!(sat; to = 1.0)
mean_anomaly(sat)
```
"""
function mean_anomaly!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                       params::NamedTuple = NamedTuple(); to::Real)
    isfinite(to) || throw(ArgumentError(
        "mean anomaly must be finite; got $to. The state was not changed."))
    kep, _ = _keplerian(sat, cs, params)
    kep.ecc == 1 && throw(ArgumentError(
        "a parabolic orbit (eccentricity 1) has no mean anomaly. The state was not changed."))
    return set_element!(sat, KeplerianState, cs, :ta, params; to = _true_from_mean(to, kep.ecc))
end

EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof(mean_anomaly),
                           cs::AbstractCoordinateSystem = frame_of(sat),
                           params::NamedTuple = NamedTuple(); to) =
    mean_anomaly!(sat, cs, params; to = to)

# --- traits ------------------------------------------------------------------

for (f, l, u) in (
        (:mean_anomaly,                    "Mean anomaly",                 "rad"),
        (:eccentric_anomaly,               "Eccentric anomaly",            "rad"),
        (:hyperbolic_anomaly,              "Hyperbolic anomaly",           "rad"),
        (:specific_energy,                 "Specific energy",              "km²/s²"),
        (:hyperbolic_excess_speed,         "Hyperbolic excess speed",      "km/s"),
        (:apoapsis_speed,                  "Apoapsis speed",               "km/s"),
        (:periapsis_speed,                 "Periapsis speed",              "km/s"),
        (:angular_momentum,                "Angular momentum",             "km²/s"),
        (:angular_momentum_magnitude,      "Angular momentum magnitude",   "km²/s"),
        (:orbit_normal,                    "Orbit normal",                 ""),
        (:orbital_angular_velocity,        "Orbital angular velocity",     "rad/s"),
        (:eccentricity_vector,             "Eccentricity vector",          ""),
        (:orbital_period,                  "Orbital period",               "s"),
        (:mean_motion,                     "Mean motion",                  "rad/s"),
        (:argument_of_latitude,            "Argument of latitude",         "rad"),
        (:longitude_of_periapsis,          "Longitude of periapsis",       "rad"),
        (:descending_node_right_ascension, "RA of descending node",        "rad"),
        (:apoapsis_radius,                 "Apoapsis radius",              "km"),
        (:periapsis_radius,                "Periapsis radius",             "km"),
        (:semilatus_rectum,                "Semi-latus rectum",            "km"),
        (:c3,                              "C3",                           "km²/s²"),
        (:true_longitude,                  "True longitude",               "rad"),
        (:b_dot_t,                         "B·T",                          "km"),
        (:b_dot_r,                         "B·R",                          "km"),
        (:b_magnitude,                     "B-vector magnitude",           "km"),
        (:b_plane_angle,                   "B-plane angle",                "rad"))
    @eval EpicycleBase.label(::typeof($f)) = $l
    isempty(u) || @eval EpicycleBase.unit(::typeof($f)) = $u
    @eval subject_type(::typeof($f)) = Union{Spacecraft, Coordinate}
end

# Angles that wrap at 2π. The hyperbolic anomaly and a hyperbola's mean anomaly do not wrap, but a
# closed orbit's mean anomaly does, and that is the case a stopping condition meets.
for f in (:mean_anomaly, :eccentric_anomaly, :argument_of_latitude, :longitude_of_periapsis,
          :descending_node_right_ascension, :true_longitude, :b_plane_angle)
    @eval EpicycleBase.is_cyclic(::typeof($f)) = true
    @eval EpicycleBase.cycle(::typeof($f))     = 2π
end
