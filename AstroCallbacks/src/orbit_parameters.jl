# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbit parameters that are two-body formulas of the current state (#11, sections 1 and 2).
#
# Each reads the subject in `cs`, or in its own frame without one, and takes the gravitational
# parameter from the coordinate system's origin, as `raan` does. A parameter undefined for the
# orbit throws an ArgumentError saying why, rather than returning NaN, because a NaN in a
# constraint or a stopping condition fails somewhere far from its cause.
#
# Every reader is generic in the state's element type, so ForwardDiff passes through.
# =============================================================================

using AstroStates: KeplerianState, ModifiedKeplerianState, CartesianState, to_vector
using LinearAlgebra: cross, dot, norm

# --- shared reads ------------------------------------------------------------

"Position, velocity and gravitational parameter of `subject` in `cs`."
function _rvμ(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    x = to_vector(CartesianState(Coordinate(subject, cs, params)))
    return x[1:3], x[4:6], get_gravparam(cs.origin)
end

"Keplerian elements and gravitational parameter of `subject` in `cs`."
function _keplerian(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    μ = get_gravparam(cs.origin)
    return KeplerianState(CartesianState(Coordinate(subject, cs, params)), μ), μ
end

_specific_energy(r, v, μ) = dot(v, v) / 2 - μ / norm(r)

_wrap2π(x) = mod(x, 2π)

# True anomaly is stored in [0, 2π); a hyperbola's lies in (-π, π), negative on the way in.
_signed_true_anomaly(ν) = ν > π ? ν - 2π : ν

_eccentric_from_true(ν, e)  = atan(sqrt(1 - e^2) * sin(ν), e + cos(ν))
_hyperbolic_from_true(ν, e) = asinh(sqrt(e^2 - 1) * sin(ν) / (1 + e * cos(ν)))

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
    if e < 1
        E = _eccentric_from_true(ν, e)
        return _wrap2π(E - e * sin(E))
    elseif e > 1
        H = _hyperbolic_from_true(_signed_true_anomaly(ν), e)
        return e * sinh(H) - H
    end
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
    return _wrap2π(_eccentric_from_true(kep.ta, kep.ecc))
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
    return _hyperbolic_from_true(_signed_true_anomaly(kep.ta), kep.ecc)
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
function specific_energy(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return _specific_energy(r, v, μ)
end

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
function hyperbolic_excess_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    ε = _specific_energy(r, v, μ)
    ε > 0 || throw(ArgumentError(
        "Hyperbolic excess speed is defined only on an escape orbit (positive specific energy); " *
        "this orbit has specific energy $ε km²/s²."))
    return sqrt(2ε)
end

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
function apoapsis_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    e = norm(_eccentricity_vector(r, v, μ))
    _require_closed(e, "Apoapsis speed")
    return μ * (1 - e) / norm(cross(r, v))
end

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
function periapsis_speed(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return μ * (1 + norm(_eccentricity_vector(r, v, μ))) / norm(cross(r, v))
end

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
    return cross(r, v)
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
    return norm(cross(r, v))
end

"""
    orbit_normal(subject[, cs][, params])

Unit vector along the angular momentum of `subject`, (r × v)/‖r × v‖, expressed in `cs` or the
subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Undefined for rectilinear motion, where r × v is zero; the result is then not finite.

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
    h = cross(r, v)
    return h ./ norm(h)
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
    return cross(r, v) ./ dot(r, r)
end

_eccentricity_vector(r, v, μ) = ((dot(v, v) - μ / norm(r)) .* r .- dot(r, v) .* v) ./ μ

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
function eccentricity_vector(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return _eccentricity_vector(r, v, μ)
end

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
function orbital_period(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    ε = _specific_energy(r, v, μ)
    ε < 0 || throw(ArgumentError(
        "Orbital period is defined only for a closed orbit (negative specific energy); this " *
        "orbit has specific energy $ε km²/s²."))
    a = -μ / (2ε)
    return 2π * sqrt(a^3 / μ)
end

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
function mean_motion(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    ε = _specific_energy(r, v, μ)
    ε == 0 && throw(ArgumentError("a parabolic orbit (zero specific energy) has no mean motion."))
    return sqrt(μ / abs(μ / (2ε))^3)
end

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
# Each is reachable through `element` on some representation today; a name makes it usable
# directly in a constraint or a report.

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
    r, v, μ = _rvμ(subject, cs, params)
    e = norm(_eccentricity_vector(r, v, μ))
    _require_closed(e, "Apoapsis radius")
    return dot(cross(r, v), cross(r, v)) / μ / (1 - e)
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
function periapsis_radius(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    e = norm(_eccentricity_vector(r, v, μ))
    return dot(cross(r, v), cross(r, v)) / μ / (1 + e)
end

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
function semilatus_rectum(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    h = cross(r, v)
    return dot(h, h) / μ
end

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
function c3(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return 2 * _specific_energy(r, v, μ)
end

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
function true_longitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    kep, _ = _keplerian(subject, cs, params)
    return _wrap2π(kep.raan + kep.aop + kep.ta)
end

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

"B·T and B·R of a hyperbolic state, in km."
function _b_plane(r, v, μ, what)
    ev = _eccentricity_vector(r, v, μ)
    e  = norm(ev)
    _require_hyperbolic(e, what)
    h  = cross(r, v)
    ĥ  = h ./ norm(h)
    ê  = ev ./ e
    Ŝ  = ê ./ e .+ sqrt(1 - 1 / e^2) .* cross(ĥ, ê)
    T  = cross(Ŝ, [zero(e), zero(e), one(e)])
    T̂  = T ./ norm(T)
    R̂  = cross(Ŝ, T̂)
    B  = (norm(h) / sqrt(2 * _specific_energy(r, v, μ))) .* cross(Ŝ, ĥ)
    return dot(B, T̂), dot(B, R̂)
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
function b_dot_t(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return _b_plane(r, v, μ, "B·T")[1]
end

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
function b_dot_r(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    return _b_plane(r, v, μ, "B·R")[2]
end

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
function b_magnitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    bt, br = _b_plane(r, v, μ, "B-vector magnitude")
    return sqrt(bt^2 + br^2)
end

"""
    b_plane_angle(subject[, cs][, params])

Angle of the B vector of `subject` from T̂ toward R̂, atan(B·R, B·T), in radians in [0, 2π),
measured in `cs` or the subject's own frame.

# Arguments
See [`position_vector`](@ref).

# Notes
Depends on the coordinate system, as [`b_dot_t`](@ref) does. Defined only for a hyperbolic orbit.

# Returns
The B-plane angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))
b_plane_angle(sat)
```
"""
function b_plane_angle(subject, cs::AbstractCoordinateSystem, params::NamedTuple)
    r, v, μ = _rvμ(subject, cs, params)
    bt, br = _b_plane(r, v, μ, "B-plane angle")
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
function _true_from_mean(M, e)
    if e < 1
        E = AstroStates._eccentric_anomaly(M, e)        # the bracketed solver the state conversions use
        return _wrap2π(2 * atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2)))
    end
    H = asinh(M / e)
    for _ in 1:50
        δ = (e * sinh(H) - H - M) / (e * cosh(H) - 1)
        H -= δ
        abs(δ) < 1e-14 * max(1, abs(H)) && break
    end
    return 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(H / 2))
end

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
