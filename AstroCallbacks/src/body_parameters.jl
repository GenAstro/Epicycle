# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbit parameters relative to the Sun and Moon (#11, section 4).
#
# A parameter tied to a body takes the body, not a coordinate system: `beta_angle(sat, moon)`.
# Without one it uses the central body of the subject's own frame. Sun and Moon positions come
# from the loaded planetary ephemeris at the subject's epoch.
#
# Every parameter is generic in the state's element type, so ForwardDiff passes through. The Sun
# and Moon positions depend on the epoch only.
# =============================================================================

using AstroFrames: CoordinateSystem, ICRF, epoch_of
using AstroUniverse: CelestialBody, translate, earth, moon, sun

# --- bodies and frames ---------------------------------------------------------

function _central_body(subject)
    body = frame_of(subject).origin
    body isa CelestialBody || throw(ArgumentError(
        "the subject's frame is centred on $(body), not a celestial body; pass the body as the " *
        "second argument."))
    return body
end

_inertial_frame(body::CelestialBody) = CoordinateSystem(body, ICRF())

_jd_tdb(subject) = epoch_of(subject).tdb.jd

_angle_between(a, b) = atan(norm(cross(a, b)), dot(a, b))

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
- `subject`: a `Spacecraft`, an `AstroFrames.Coordinate`, or anything else implementing the
  subject interface.
- `body`: the `CelestialBody` whose orbit plane is meant. Without it, the central body of the
  subject's own frame.

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
See [`beta_angle`](@ref).

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
See [`beta_angle`](@ref).

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

# --- the shorter methods, and traits -------------------------------------------

const _BODY_PARAMETERS = (:beta_angle, :lunar_beta_angle, :orbit_normal_sun_angle)

for fn in _BODY_PARAMETERS
    @eval $fn(subject) = $fn(subject, _central_body(subject))
end

for (f, l, u) in (
        (:beta_angle,                  "Beta angle",                     "rad"),
        (:lunar_beta_angle,            "Lunar beta angle",               "rad"),
        (:orbit_normal_sun_angle,      "Orbit-normal-to-Sun angle",      "rad"),
        (:sun_earth_spacecraft_angle,  "Sun–Earth–spacecraft angle",     "rad"),
        (:sun_spacecraft_earth_angle,  "Sun–spacecraft–Earth angle",     "rad"),
        (:sun_spacecraft_moon_angle,   "Sun–spacecraft–Moon angle",      "rad"),
        (:moon_earth_sun_angle,        "Moon–Earth–Sun angle",           "rad"),
        (:moon_earth_spacecraft_angle, "Moon–Earth–spacecraft angle",    "rad"),
        (:moon_spacecraft_earth_angle, "Moon–spacecraft–Earth angle",    "rad"))
    @eval EpicycleBase.label(::typeof($f)) = $l
    isempty(u) || @eval EpicycleBase.unit(::typeof($f)) = $u
    @eval subject_type(::typeof($f)) = Union{Spacecraft, Coordinate}
end
