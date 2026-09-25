# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A ground station is a participant in a measurement: it has a place on a body,
# it moves as that body rotates, and it can see a spacecraft or not. That makes
# it a model, and it owns quantities the way a spacecraft does.

"""
    AbstractGeodeticReference

`AbstractGeodeticReference` is the common type for shapes against which
geodetic latitude, longitude, and altitude are measured. [`Ellipsoid`](@ref)
uses the body's equatorial radius and flattening.

# Example
```jldoctest
Ellipsoid() isa AbstractGeodeticReference

# output
true
```
"""
abstract type AbstractGeodeticReference end

"""
    Ellipsoid()

An oblate ellipsoid of revolution, taken from the body's equatorial radius and
flattening.

# Example
```jldoctest
Ellipsoid() isa Ellipsoid

# output
true
```
"""
struct Ellipsoid <: AbstractGeodeticReference end

"""
    GroundStation(; name, body, latitude, longitude, altitude,
                  reference = Ellipsoid(), min_elevation = -90.0)

A fixed site on a rotating body, from which a spacecraft is tracked.

# Arguments
- `name`: How the station is identified, including in a tracking data file.
- `body`: The body it sits on.
- `latitude`, `longitude`: Geodetic, in degrees. Latitude in `[-90, 90]`.
- `altitude`: Height above the reference shape, in km.
- `reference`: The shape latitude is measured against. See
  [`AbstractGeodeticReference`](@ref).
- `min_elevation`: How far above the horizon a spacecraft must be to be seen,
  in degrees. In `[-90, 90)`.

# Fields
- `name::String`: Station identifier.
- `body::CelestialBody`: Body on which the station is fixed.
- `reference::AbstractGeodeticReference`: Shape used for geodetic coordinates.
- `latitude::Float64`: Geodetic latitude in degrees.
- `longitude::Float64`: Geodetic longitude in degrees.
- `altitude::Float64`: Height above the reference shape in km.
- `min_elevation_deg::Float64`: Visibility cutoff in degrees.

# Returns
A `GroundStation`. [`get_state`](@ref) evaluates its position and velocity at
an epoch, and [`is_visible`](@ref) applies its elevation cutoff.

# Example
```julia
goldstone = GroundStation(name = "DSS-14", body = earth,
                          latitude = 35.4267, longitude = -116.89,
                          altitude = 1.0, min_elevation = 10.0)
```
"""
struct GroundStation{B<:CelestialBody, R<:AbstractGeodeticReference}
    name::String
    body::B
    reference::R
    latitude::Float64
    longitude::Float64
    altitude::Float64
    min_elevation_deg::Float64
end

function GroundStation(; name::AbstractString,
                         body::CelestialBody,
                         reference::AbstractGeodeticReference = Ellipsoid(),
                         latitude::Real,
                         longitude::Real,
                         altitude::Real,
                         min_elevation::Real = -90.0)
    me = Float64(min_elevation)
    -90.0 <= me < 90.0 || throw(ArgumentError(
        "GroundStation: min_elevation must lie in [-90, 90) degrees; got $me"))
    -90.0 <= latitude <= 90.0 || throw(ArgumentError(
        "GroundStation: latitude must lie in [-90, 90] degrees; got $latitude"))
    return GroundStation(String(name), body, reference,
                         Float64(latitude), Float64(longitude),
                         Float64(altitude), me)
end

"""
    geodetic_to_body_fixed(body, lat_deg, lon_deg, alt_km, ::Ellipsoid)

Geodetic coordinates to a body-fixed position vector, in km.

# Returns
The closed-form conversion, using the body's equatorial radius and flattening.
"""
function geodetic_to_body_fixed(body::CelestialBody,
                                lat_deg::Real, lon_deg::Real, alt_km::Real,
                                ::Ellipsoid)
    a  = float(body.equatorial_radius)
    f  = float(body.flattening)
    e2 = 2f - f^2

    φ = deg2rad(lat_deg)
    λ = deg2rad(lon_deg)
    sφ, cφ = sincos(φ)
    sλ, cλ = sincos(λ)

    N = a / sqrt(1 - e2 * sφ^2)
    h = float(alt_km)

    return SVector{3,Float64}((N + h) * cφ * cλ,
                              (N + h) * cφ * sλ,
                              (N * (1 - e2) + h) * sφ)
end

"""
    get_state(gs::GroundStation, t; axes = GCRF())

Where the station is, and how fast it is moving, at an epoch.

# Arguments
- `gs`: The station.
- `t`: The epoch.
- `axes`: The axes to express the answer in. Body-fixed by construction, so an
  inertial answer carries the velocity the body's rotation gives it.

# Returns
A tuple of position and velocity, in km and km/s.

# Notes
The station is fixed in body-fixed axes, so its velocity there is zero and
everything an inertial observer sees comes out of the frame conversion.

Which precession-nutation theory applies, and so which Earth orientation data,
is decided by the pair of axes rather than passed by the caller.

# Example
```julia
using AstroModels, AstroUniverse, AstroEpochs, AstroFrames
goldstone = GroundStation(name = "DSS-14", body = earth,
                          latitude = 35.4267, longitude = -116.89, altitude = 1.0)
epoch = Time("2020-06-01T00:00:00", UTC(), ISOT())
r, v = get_state(goldstone, epoch)                 # GCRF
r, v = get_state(goldstone, epoch; axes = ITRF())  # body-fixed, v is zero
```
"""
function get_state(gs::GroundStation, t; axes::AbstractAxes = GCRF())
    r_bf = geodetic_to_body_fixed(gs.body, gs.latitude, gs.longitude,
                                  gs.altitude, gs.reference)
    v_bf = SVector{3,Float64}(0.0, 0.0, 0.0)

    fixed = Coordinate(vcat(r_bf, v_bf), CoordinateSystem(gs.body, ITRF()), t)
    out   = Coordinate(fixed, CoordinateSystem(gs.body, axes))
    x     = to_vector(CartesianState(out))
    return SVector{3,Float64}(x[1], x[2], x[3]), SVector{3,Float64}(x[4], x[5], x[6])
end

"""
    is_visible(gs::GroundStation, r, t) -> Bool

Whether a position is above the station's minimum elevation at an epoch.

# Arguments
- `gs`: The station.
- `r`: The position to check, in km, in GCRF, the axes `get_state` returns by default.
- `t`: The epoch.

# Returns
`true` when the elevation angle is at or above `min_elevation`.

# Notes
Elevation is measured from the local horizontal, the plane normal to the
geodetic vertical at the station.

# Example
```jldoctest
using AstroUniverse: earth
using AstroEpochs: Time, UTC, ISOT
station = GroundStation(name = "Equator", body = earth,
                        latitude = 0.0, longitude = 0.0, altitude = 0.0)
epoch = Time("2020-06-01T00:00:00", UTC(), ISOT())
r_station, _ = get_state(station, epoch)
is_visible(station, 1.1 .* r_station, epoch)

# output
true
```
"""
function is_visible(gs::GroundStation, r::AbstractVector, t)
    # Elevation is an angle, so it is the same in any axes; body-fixed is where the station and
    # its vertical are constant.
    inertial = Coordinate(vcat(SVector{3,Float64}(r[1], r[2], r[3]), zero(SVector{3,Float64})),
                          CoordinateSystem(gs.body, GCRF()), t)
    r_bf  = to_vector(CartesianState(Coordinate(inertial, CoordinateSystem(gs.body, ITRF()))))
    r_gs  = geodetic_to_body_fixed(gs.body, gs.latitude, gs.longitude, gs.altitude, gs.reference)
    ρ     = SVector{3,Float64}(r_bf[1], r_bf[2], r_bf[3]) - r_gs
    nρ    = norm(ρ)
    nρ > 0.0 || return false
    sinel = clamp(dot(ρ, _local_vertical(gs.latitude, gs.longitude, gs.reference)) / nρ, -1.0, 1.0)
    return rad2deg(asin(sinel)) >= gs.min_elevation_deg
end

# The outward normal to the reference shape at the station, in body-fixed axes. On an
# ellipsoid that is the geodetic vertical, which is what geodetic latitude is measured from.
function _local_vertical(lat_deg::Real, lon_deg::Real, ::Ellipsoid)
    sφ, cφ = sincosd(lat_deg)
    sλ, cλ = sincosd(lon_deg)
    return SVector{3,Float64}(cφ * cλ, cφ * sλ, sφ)
end
