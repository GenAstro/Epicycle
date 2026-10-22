# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbit parameters relative to a body's shape and rotation, and to the Sun and Moon
# (#11, sections 3 and 4).
#
# A parameter tied to a body takes the body, not a coordinate system: `altitude(sat, moon)`.
# Without one it uses the central body of the subject's own frame. Positions relative to the
# body's surface are read in its body-fixed frame: ITRF for the Earth, the Moon's mean-Earth
# axes, and the IAU body-fixed axes for the rest. Sun and Moon positions come from the loaded
# planetary ephemeris at the subject's epoch.
#
# Every parameter is generic in the state's element type, so ForwardDiff passes through. The Sun
# and Moon positions depend on the epoch only.
# =============================================================================

using AstroFrames: CoordinateSystem, ICRF, CelestialBodyFixed, EarthFixed, MODEq, MoonFixed, epoch_of
using AstroUniverse: CelestialBody, translate, earth, moon, sun

# --- bodies and frames ---------------------------------------------------------

function _central_body(subject)
    body = frame_of(subject).origin
    body isa CelestialBody || throw(ArgumentError(
        "the subject's frame is centred on $(body), not a celestial body; pass the body as the " *
        "second argument."))
    return body
end

_fixed_frame(body::CelestialBody) =
    body.naifid == 399 ? EarthFixed :
    body.naifid == 301 ? MoonFixed  : CoordinateSystem(body, CelestialBodyFixed())

_inertial_frame(body::CelestialBody) = CoordinateSystem(body, ICRF())

_jd_tdb(subject) = epoch_of(subject).tdb.jd

"Position and velocity of `subject` in the body-fixed frame of `body`."
function _fixed_rv(subject, body::CelestialBody)
    x = to_vector(CartesianState(Coordinate(subject, _fixed_frame(body), NamedTuple())))
    return x[1:3], x[4:6]
end

"A body-centred ICRF position `r` at `subject`'s epoch, re-expressed in the body-fixed axes."
function _to_fixed(r, subject, body::CelestialBody)
    c = Coordinate([r; zero(r)], _inertial_frame(body), epoch_of(subject))
    return to_vector(CartesianState(Coordinate(c, _fixed_frame(body), NamedTuple())))[1:3]
end

"A body-fixed position `r` at `subject`'s epoch, re-expressed in body-centred ICRF axes."
function _from_fixed(r, subject, body::CelestialBody)
    c = Coordinate([r; zero(r)], _fixed_frame(body), epoch_of(subject))
    return to_vector(CartesianState(Coordinate(c, _inertial_frame(body), NamedTuple())))[1:3]
end

_angle_between(a, b) = atan(norm(cross(a, b)), dot(a, b))

# --- geodetic position ---------------------------------------------------------

"""
Geodetic latitude and height of body-fixed position `r` on an ellipsoid of equatorial radius `R`
and flattening `f`. Fixed-point iteration on the latitude; eight passes reach well below a
millimetre for any orbit above the surface, and a fixed count keeps it differentiable.
"""
function _geodetic(r, R, f)
    e2 = f * (2 - f)
    p  = sqrt(r[1]^2 + r[2]^2)
    φ  = atan(r[3], p * (1 - e2))
    for _ in 1:8
        s = sin(φ)
        N = R / sqrt(1 - e2 * s^2)
        φ = atan(r[3] + e2 * N * s, p)
    end
    s = sin(φ)
    h = p * cos(φ) + r[3] * s - R * sqrt(1 - e2 * s^2)
    return φ, h
end

"""
    altitude(subject[, body])

Height of `subject` above the reference ellipsoid of `body`, along the ellipsoid normal, in km.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.
- `body`: the `CelestialBody` whose ellipsoid the height is measured from. Without it, the
  central body of the subject's own frame.

# Notes
The ellipsoid is the body's equatorial radius and flattening in AstroUniverse; a body with zero
flattening is a sphere and the altitude is r − R.

# Returns
The altitude in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
altitude(sat)
altitude(sat, earth)
```
"""
function altitude(subject, body::CelestialBody)
    r, _ = _fixed_rv(subject, body)
    return _geodetic(r, body.equatorial_radius, body.flattening)[2]
end

"""
    geodetic_latitude(subject[, body])

Geodetic latitude of the point below `subject` on the reference ellipsoid of `body`, in radians
in [−π/2, π/2].

# Arguments
See [`altitude`](@ref).

# Notes
The angle between the ellipsoid normal through the subject and the body's equatorial plane. On a
spherical body it equals the geocentric latitude.

# Returns
The geodetic latitude in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
geodetic_latitude(sat)
```
"""
function geodetic_latitude(subject, body::CelestialBody)
    r, _ = _fixed_rv(subject, body)
    return _geodetic(r, body.equatorial_radius, body.flattening)[1]
end

"""
    longitude(subject[, body])

East longitude of the point below `subject`, in the body-fixed frame of `body`, in radians in
(−π, π].

# Arguments
See [`altitude`](@ref).

# Returns
The longitude in radians, positive east of the prime meridian.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
longitude(sat)
```
"""
function longitude(subject, body::CelestialBody)
    r, _ = _fixed_rv(subject, body)
    return atan(r[2], r[1])
end

"""
    velocity_azimuth(subject[, body])

Azimuth of the body-fixed velocity of `subject`, measured clockwise from local north toward
east, in radians in [0, 2π).

# Arguments
See [`altitude`](@ref).

# Notes
North and east are those of the ellipsoid at the point below the subject. The velocity is
relative to the rotating body, so it includes the body's rotation.

# Returns
The velocity azimuth in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
velocity_azimuth(sat)
```
"""
function velocity_azimuth(subject, body::CelestialBody)
    r, v = _fixed_rv(subject, body)
    φ, _ = _geodetic(r, body.equatorial_radius, body.flattening)
    λ = atan(r[2], r[1])
    east  = (-sin(λ) * v[1] + cos(λ) * v[2])
    north = (-sin(φ) * cos(λ) * v[1] - sin(φ) * sin(λ) * v[2] + cos(φ) * v[3])
    return _wrap2π(atan(east, north))
end

# --- the apsides above the ellipsoid -------------------------------------------

function _apsis_altitude(subject, body::CelestialBody, which::Symbol)
    r, v, μ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    ev = eccentricity_vector(r, v, μ)
    e  = norm(ev)
    p  = dot(cross(r, v), cross(r, v)) / μ
    if which === :apoapsis
        _require_closed(e, "Apoapsis altitude")
        radius, dir = p / (1 - e), -ev ./ e
    else
        radius, dir = p / (1 + e), ev ./ e
    end
    rf = _to_fixed(radius .* dir, subject, body)
    return _geodetic(rf, body.equatorial_radius, body.flattening)[2]
end

"""
    apoapsis_altitude(subject[, body])

Height of the osculating apoapsis of `subject` above the reference ellipsoid of `body`, in km.

# Arguments
See [`altitude`](@ref).

# Notes
The apoapsis point is placed on the ellipsoid using the body's equator at the subject's epoch,
so on an oblate body the result is the apoapsis radius less the ellipsoid radius at the
apoapsis latitude, not less the equatorial radius. Defined only for a closed orbit; a hyperbola
throws, as does a circular orbit, which has no apoapsis direction.

# Returns
The apoapsis altitude in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
apoapsis_altitude(sat)
```
"""
apoapsis_altitude(subject, body::CelestialBody) = _apsis_altitude(subject, body, :apoapsis)

"""
    periapsis_altitude(subject[, body])

Height of the osculating periapsis of `subject` above the reference ellipsoid of `body`, in km.

# Arguments
See [`altitude`](@ref).

# Notes
The periapsis point is placed on the ellipsoid using the body's equator at the subject's epoch,
as [`apoapsis_altitude`](@ref) does. Defined for every conic except a circular orbit, which has
no periapsis direction.

# Returns
The periapsis altitude in km.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
periapsis_altitude(sat)
```
"""
periapsis_altitude(subject, body::CelestialBody) = _apsis_altitude(subject, body, :periapsis)

# --- rotation ------------------------------------------------------------------

# The body's equator frame, in which the meridian hour angle and the mean geographic longitude are
# both measured, so they share one plane and one zero. Its Z axis is the body's pole. Its X axis
# is the mean equinox of date for the Earth, so the hour angle is Greenwich mean sidereal
# time; for any other body it is the ascending node of the body's equator on the ICRF equator,
# the zero of the IAU prime-meridian angle W.
"Axes of the equator frame of `body` at `subject`'s epoch, as unit vectors in body-centred ICRF."
function _equator_axes(subject, body::CelestialBody)
    if body.naifid == 399
        mod_(e) = to_vector(CartesianState(Coordinate(
            Coordinate([e; 0.0; 0.0; 0.0], CoordinateSystem(earth, MODEq()), epoch_of(subject)),
            _inertial_frame(earth), NamedTuple())))[1:3]
        return mod_([1.0, 0.0, 0.0]), mod_([0.0, 1.0, 0.0]), mod_([0.0, 0.0, 1.0])
    end
    ẑ = _from_fixed([0.0, 0.0, 1.0], subject, body)
    n = cross([0.0, 0.0, 1.0], ẑ)
    x̂ = norm(n) > 1e-12 ? n ./ norm(n) : [1.0, 0.0, 0.0]   # pole on the ICRF pole: no node
    return x̂, cross(ẑ, x̂), ẑ
end

"""
    meridian_hour_angle(subject[, body])

Angle of the prime meridian of `body` measured eastward in its equator, at the epoch of
`subject`, in radians in [0, 2π).

# Arguments
See [`altitude`](@ref). Only the subject's epoch is used.

# Notes
- For the Earth it is Greenwich mean sidereal time, measured in the mean equator from the mean
  equinox of date. It differs from apparent sidereal time by the equation of the equinoxes, at
  most about 1 s; [`longitude`](@ref) and [`true_local_time`](@ref) use the true orientation.
- For any other body it is the IAU prime-meridian angle W, measured from the ascending node of
  the body's equator on the ICRF equator.
- [`local_sidereal_time`](@ref) and [`mean_geographic_longitude`](@ref) are measured from the
  same zero.

# Returns
The meridian hour angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
meridian_hour_angle(sat)
```
"""
function meridian_hour_angle(subject, body::CelestialBody)
    x̂, ŷ, _ = _equator_axes(subject, body)
    meridian = _from_fixed([1.0, 0.0, 0.0], subject, body)
    return _wrap2π(atan(dot(meridian, ŷ), dot(meridian, x̂)))
end

"""
    local_sidereal_time(subject[, body])

Sidereal time at the longitude of `subject`, the meridian hour angle plus the east longitude,
in radians in [0, 2π).

# Arguments
See [`altitude`](@ref).

# Notes
Measured from the same zero as [`meridian_hour_angle`](@ref): for the Earth, Greenwich mean
sidereal time plus the east longitude; for another body, W plus the east longitude.

# Returns
The local sidereal time in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
local_sidereal_time(sat)
```
"""
local_sidereal_time(subject, body::CelestialBody) =
    _wrap2π(meridian_hour_angle(subject, body) + longitude(subject, body))

"""
    mean_geographic_longitude(subject[, body])

Mean geographic longitude of `subject`, Ω + ω + M less the meridian hour angle, in radians in
[0, 2π).

# Arguments
See [`altitude`](@ref).

# Notes
The mean longitude of the orbit, measured from the body's prime meridian rather than an inertial
axis; on a geostationary orbit it is the station longitude. Ω, ω and M are measured in the
body's equator from the same zero as [`meridian_hour_angle`](@ref), so the difference is taken in
one plane. Defined only for a closed orbit.

# Returns
The mean geographic longitude in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([42164.0, 0.0, 0.0, 0.0, 3.0747, 0.0]))
mean_geographic_longitude(sat)
```
"""
function mean_geographic_longitude(subject, body::CelestialBody)
    r, v, μ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    x̂, ŷ, ẑ = _equator_axes(subject, body)
    to_equator(w) = [dot(w, x̂), dot(w, ŷ), dot(w, ẑ)]
    kep = KeplerianState(CartesianState([to_equator(r); to_equator(v)]), μ)
    _require_closed(kep.ecc, "Mean geographic longitude")
    M = true_to_mean_anomaly(kep.ta, kep.ecc)
    return _wrap2π(kep.raan + kep.aop + M - meridian_hour_angle(subject, body))
end

"""
    nodal_precession_rate(subject[, body]; j2)

Secular rate of the right ascension of the ascending node of `subject` from the body's J2,
−(3/2) n J2 (R/p)² cos i, in rad/s.

# Arguments
See [`altitude`](@ref).
- `j2`: the body's unnormalized J2. Defaults to the Earth's, 1.08262668e-3 (EGM96); another
  body needs it given.

# Notes
The inclination is measured from the body's equator, through its pole at the subject's epoch.
Negative for a prograde orbit, positive for a retrograde one; a sun-synchronous orbit has
1.991e-7 rad/s. Defined only for a closed orbit.

# Returns
The nodal precession rate in rad/s.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
nodal_precession_rate(sat)
```
"""
function nodal_precession_rate(subject, body::CelestialBody; j2 = nothing)
    J2 = j2 === nothing ? _default_j2(body) : j2
    r, v, μ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    ε = specific_energy(r, v, μ)
    e = norm(eccentricity_vector(r, v, μ))
    _require_closed(e, "Nodal precession rate")
    a = -μ / (2ε)
    h = cross(r, v)
    p = dot(h, h) / μ
    pole = _from_fixed([0.0, 0.0, 1.0], subject, body)
    cosi = dot(h, pole) / (norm(h) * norm(pole))
    return -1.5 * sqrt(μ / a^3) * J2 * (body.equatorial_radius / p)^2 * cosi
end

function _default_j2(body::CelestialBody)
    body.naifid == 399 && return 1.08262668e-3
    throw(ArgumentError(
        "no default J2 for $(body.name); pass it as nodal_precession_rate(sat, body; j2 = value)."))
end

# --- Sun and Moon geometry (#11, section 4) ------------------------------------

_sun_from(body::CelestialBody, subject)  = translate(body, sun, _jd_tdb(subject))
_moon_from(body::CelestialBody, subject) = translate(body, moon, _jd_tdb(subject))

"Position of `subject` from the Earth, ICRF axes."
_from_earth(subject) = _rvμ(subject, _inertial_frame(earth), NamedTuple())[1]

"""
    beta_angle(subject[, body])

Beta angle of `subject`: the elevation of the Sun direction above its orbit plane about `body`,
in radians in [−π/2, π/2].

# Arguments
See [`altitude`](@ref).

# Notes
Positive when the Sun is on the side of the orbit normal (r × v). The Sun direction is taken
from the body's centre.

# Returns
The beta angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
beta_angle(sat)
```
"""
function beta_angle(subject, body::CelestialBody)
    body.naifid == 10 && throw(ArgumentError("the beta angle needs an orbit about a body other than the Sun."))
    r, v, _ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    h = cross(r, v)
    s = _sun_from(body, subject)
    return asin(clamp(dot(h, s) / (norm(h) * norm(s)), -1, 1))
end

"""
    lunar_beta_angle(subject[, body])

Elevation of the Moon direction above the orbit plane of `subject` about `body`, in radians in
[−π/2, π/2].

# Arguments
See [`altitude`](@ref).

# Notes
The beta angle with the Moon in place of the Sun. Positive when the Moon is on the side of the
orbit normal. The Moon direction is taken from the body's centre.

# Returns
The lunar beta angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
lunar_beta_angle(sat)
```
"""
function lunar_beta_angle(subject, body::CelestialBody)
    body.naifid == 301 && throw(ArgumentError("the lunar beta angle needs an orbit about a body other than the Moon."))
    r, v, _ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    h = cross(r, v)
    m = _moon_from(body, subject)
    return asin(clamp(dot(h, m) / (norm(h) * norm(m)), -1, 1))
end

"""
    orbit_normal_sun_angle(subject[, body])

Angle between the orbit normal of `subject` about `body` and the line from `subject` to the
Sun, in radians in [0, π].

# Arguments
See [`altitude`](@ref).

# Notes
Measured from the spacecraft rather than the body's centre, which is what differs from the
complement of [`beta_angle`](@ref).

# Returns
The angle in radians.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
orbit_normal_sun_angle(sat)
```
"""
function orbit_normal_sun_angle(subject, body::CelestialBody)
    r, v, _ = _rvμ(subject, _inertial_frame(body), NamedTuple())
    return _angle_between(cross(r, v), _sun_from(body, subject) .- r)
end

# The six angles among the Sun, the Earth, the Moon and the spacecraft. Each is named for its
# three points with the vertex in the middle: sun_earth_spacecraft_angle is the angle at the
# Earth between the Sun and the spacecraft.
for (fn, vertex, a, b, doc_a, doc_v, doc_b) in (
        (:sun_earth_spacecraft_angle,  :earth, :sun,  :sc,    "the Sun",  "the Earth",      "the spacecraft"),
        (:sun_spacecraft_earth_angle,  :sc,    :sun,  :earth, "the Sun",  "the spacecraft", "the Earth"),
        (:sun_spacecraft_moon_angle,   :sc,    :sun,  :moon,  "the Sun",  "the spacecraft", "the Moon"),
        (:moon_earth_sun_angle,        :earth, :moon, :sun,   "the Moon", "the Earth",      "the Sun"),
        (:moon_earth_spacecraft_angle, :earth, :moon, :sc,    "the Moon", "the Earth",      "the spacecraft"),
        (:moon_spacecraft_earth_angle, :sc,    :moon, :earth, "the Moon", "the spacecraft", "the Earth"))
    doc = """
        $fn(subject)

    Angle at $doc_v between $doc_a and $doc_b, in radians in [0, π], at the epoch of `subject`.

    # Arguments
    - `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
      subject interface. Its position is used wherever it is centred.

    # Returns
    The angle in radians.

    # Examples
    ```julia
    sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
    $fn(sat)
    ```
    """
    @eval begin
        function $fn(subject)
            pts = (earth = zeros(3), sc = _from_earth(subject),
                   sun = _sun_from(earth, subject), moon = _moon_from(earth, subject))
            o = pts.$vertex
            return _angle_between(pts.$a .- o, pts.$b .- o)
        end
        @doc $doc $fn
    end
end

# --- local time ----------------------------------------------------------------

"""
    true_local_time(subject[, body])

Apparent solar time at the point below `subject` on `body`, in hours in [0, 24).

# Arguments
See [`altitude`](@ref).

# Notes
Twelve hours plus the hour angle between the subject's meridian and the Sun's, both measured in
the body-fixed frame, so it is 12 h with the Sun overhead.

# Returns
The true local time in hours.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
true_local_time(sat)
```
"""
function true_local_time(subject, body::CelestialBody)
    s = _to_fixed(_sun_from(body, subject), subject, body)
    Δλ = longitude(subject, body) - atan(s[2], s[1])
    return mod(12 + Δλ * 12 / π, 24)
end

"""
    mean_local_time(subject)

Mean solar time at the longitude of `subject` on the Earth, UTC plus the east longitude at 15°
per hour, in hours in [0, 24).

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.

# Notes
Defined for the Earth, where it is the time a sun-synchronous orbit's local time of the node is
quoted in. It differs from [`true_local_time`](@ref) by the equation of time, up to about 16
minutes. UTC stands in for UT1, which it follows to within 0.9 s.

# Returns
The mean local time in hours.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
mean_local_time(sat)
```
"""
function mean_local_time(subject)
    jd = epoch_of(subject).utc.jd
    utc_hours = mod(jd - 0.5, 1.0) * 24
    return mod(utc_hours + longitude(subject, earth) * 12 / π, 24)
end

# --- shadow --------------------------------------------------------------------

"""
Fraction of the Sun's disk hidden by an occulting sphere, seen from the origin of `s` and `o`:
`s` and `o` point from the observer to the Sun and to the occulter. Circular disks; the area of
overlap of two circles of angular radius `a` (Sun) and `b` (occulter) at separation `c`.
"""
function _shadow_fraction(s, o, R_sun, R_occ)
    norm(o) < norm(s) || return zero(eltype(s))          # the occulter is behind the Sun
    a = asin(R_sun / norm(s))
    b = asin(clamp(R_occ / norm(o), -1, 1))
    c = _angle_between(s, o)
    c >= a + b && return zero(c)                         # full sun
    c <= b - a && return one(c)                          # umbra
    c <= a - b && return (b / a)^2                       # annular, the occulter inside the disk
    x = (c^2 + a^2 - b^2) / (2c)
    y = sqrt(max(a^2 - x^2, zero(x)))
    A = a^2 * acos(clamp(x / a, -1, 1)) + b^2 * acos(clamp((c - x) / b, -1, 1)) - c * y
    return A / (π * a^2)
end

"""
    earth_shadow_fraction(subject)

Fraction of the Sun's disk hidden by the Earth, seen from `subject`: 0 in full sun, 1 in umbra,
in between in penumbra.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.

# Notes
The Sun and the Earth are spheres of their equatorial radii; the Earth's atmosphere and
oblateness are not modeled.

# Returns
The hidden fraction, from 0 to 1.

# Examples
```julia
sat = Spacecraft(state = CartesianState([-7000.0, 0.0, 0.0, 0.0, -7.5, 0.0]))
earth_shadow_fraction(sat)
```
"""
function earth_shadow_fraction(subject)
    r = _from_earth(subject)
    return _shadow_fraction(_sun_from(earth, subject) .- r, -r, sun.equatorial_radius,
                            earth.equatorial_radius)
end

"""
    moon_shadow_fraction(subject)

Fraction of the Sun's disk hidden by the Moon, seen from `subject`: 0 in full sun, 1 in umbra,
in between in penumbra.

# Arguments
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.

# Notes
The Sun and the Moon are spheres of their equatorial radii.

# Returns
The hidden fraction, from 0 to 1.

# Examples
```julia
sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))
moon_shadow_fraction(sat)
```
"""
function moon_shadow_fraction(subject)
    r = _from_earth(subject)
    return _shadow_fraction(_sun_from(earth, subject) .- r, _moon_from(earth, subject) .- r,
                            sun.equatorial_radius, moon.equatorial_radius)
end

# --- the shorter methods, and traits -------------------------------------------

const _BODY_PARAMETERS = (
    :altitude, :geodetic_latitude, :longitude, :velocity_azimuth, :apoapsis_altitude,
    :periapsis_altitude, :meridian_hour_angle, :local_sidereal_time, :mean_geographic_longitude,
    :beta_angle, :lunar_beta_angle, :orbit_normal_sun_angle, :true_local_time)

for fn in _BODY_PARAMETERS
    @eval $fn(subject) = $fn(subject, _central_body(subject))
end
nodal_precession_rate(subject; j2 = nothing) =
    nodal_precession_rate(subject, _central_body(subject); j2 = j2)

for (f, l, u) in (
        (:altitude,                    "Altitude",                       "km"),
        (:geodetic_latitude,           "Geodetic latitude",              "rad"),
        (:longitude,                   "Longitude",                      "rad"),
        (:velocity_azimuth,            "Velocity azimuth",               "rad"),
        (:apoapsis_altitude,           "Apoapsis altitude",              "km"),
        (:periapsis_altitude,          "Periapsis altitude",             "km"),
        (:meridian_hour_angle,         "Meridian hour angle",            "rad"),
        (:local_sidereal_time,         "Local sidereal time",            "rad"),
        (:mean_geographic_longitude,   "Mean geographic longitude",      "rad"),
        (:nodal_precession_rate,       "Nodal precession rate",          "rad/s"),
        (:beta_angle,                  "Beta angle",                     "rad"),
        (:lunar_beta_angle,            "Lunar beta angle",               "rad"),
        (:orbit_normal_sun_angle,      "Orbit-normal-to-Sun angle",      "rad"),
        (:sun_earth_spacecraft_angle,  "Sun–Earth–spacecraft angle",     "rad"),
        (:sun_spacecraft_earth_angle,  "Sun–spacecraft–Earth angle",     "rad"),
        (:sun_spacecraft_moon_angle,   "Sun–spacecraft–Moon angle",      "rad"),
        (:moon_earth_sun_angle,        "Moon–Earth–Sun angle",           "rad"),
        (:moon_earth_spacecraft_angle, "Moon–Earth–spacecraft angle",    "rad"),
        (:moon_spacecraft_earth_angle, "Moon–spacecraft–Earth angle",    "rad"),
        (:true_local_time,             "True local time",                "h"),
        (:mean_local_time,             "Mean local time",                "h"),
        (:earth_shadow_fraction,       "Earth shadow fraction",          ""),
        (:moon_shadow_fraction,        "Moon shadow fraction",           ""))
    @eval EpicycleBase.label(::typeof($f)) = $l
    isempty(u) || @eval EpicycleBase.unit(::typeof($f)) = $u
    @eval subject_type(::typeof($f)) = Union{Spacecraft, Coordinate}
end

for (f, c) in ((:longitude, 2π), (:velocity_azimuth, 2π), (:meridian_hour_angle, 2π),
               (:local_sidereal_time, 2π), (:mean_geographic_longitude, 2π),
               (:true_local_time, 24.0), (:mean_local_time, 24.0))
    @eval EpicycleBase.is_cyclic(::typeof($f)) = true
    @eval EpicycleBase.cycle(::typeof($f))     = $c
end
