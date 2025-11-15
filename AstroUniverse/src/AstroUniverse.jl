# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

"""
Module containing models such as celestial bodies, ephemerides, and related utilities.
"""
module AstroUniverse

using SPICE
using Scratch
using Downloads

using EpicycleBase

import Base: show

export CelestialBody, translate
export sun, mercury, venus, earth, moon, mars, jupiter
export saturn, uranus, neptune, pluto

export get_gravparam, set_gravparam!

"""
    function __init__()

Load SPICE kernels from managed scratch space on module initialization.
"""
function __init__()
    # Get managed cache directory for SPICE kernels
    kernel_cache = @get_scratch!("spice_kernels")
    
    # Download and load kernels
    ensure_kernel(kernel_cache, "naif0012.tls", 
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls")
    ensure_kernel(kernel_cache, "de440.bsp", 
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp")
end

"""
    ensure_kernel(cache_dir, filename, url)

Download SPICE kernel if it doesn't exist and load it.
"""
function ensure_kernel(cache_dir, filename, url)
    kernel_path = joinpath(cache_dir, filename)
    if !isfile(kernel_path)
        @info "Downloading SPICE kernel: $filename"
        Downloads.download(url, kernel_path)
    end
    furnsh(kernel_path)
end

""" 
    const EARTH_DEFAULTS

Default physical parameters for Earth CelestialBody.
"""
const EARTH_DEFAULTS = (
    name = "Earth",
    mu = 398600.4418,
    equatorial_radius = 6378.137,
    flattening = 0.00335281,
    naifid = 399,
)

"""
    CelestialBody(name::AbstractString, mu::Real, equatorial_radius::Real, flattening::Real, naifid::Integer)

Represents a celestial body with physical parameters.

Fields (units):
- name::String
- mu::T — gravitational parameter 
- equatorial_radius::T — equatorial radius 
- flattening::T — geometric flattening 
- naifid::Int — NAIF body ID

# Notes:
- Numeric fields (mu, equatorial_radius, flattening) are promoted to a common element type `T`
  to ensure type stability (e.g., passing a BigFloat will promote the others to BigFloat).
- Units default to km and seconds for built-in celestial bodies. 
  If changing units, be consistent throughout the simulation.

# Examples
```julia
using AstroUniverse
moon_like = CelestialBody(name="MyMoon", 
                                 mu=4902.8, 
                                 equatorial_radius=1737.4, 
                                 flattening=0.0,
                                 naifid=301);
show(moon_like)

# output
CelestialBody:
  name               = MyMoon
  μ                  = 4902.8
  Equatorial Radius  = 1737.4
  Flattening         = 0.0
  NAIF ID            = 301
```
"""
mutable struct CelestialBody{T<:Real} <: AbstractPoint
    name::String
    mu::T
    equatorial_radius::T
    flattening::T
    naifid::Int

    function CelestialBody{T}(
        name::String,
        mu::T,
        equatorial_radius::T,
        flattening::T,
        naifid::Int,
    ) where {T<:Real}
        if !isfinite(mu) || mu <= 0
            throw(ArgumentError("CelestialBody: μ must be finite and > 0; got $(mu)."))
        end
        if !isfinite(equatorial_radius) || equatorial_radius <= 0
            throw(ArgumentError("CelestialBody: equatorial_radius must be finite and > 0; got $(equatorial_radius)."))
        end
        if !isfinite(flattening) || flattening < 0 || flattening >= 1
            throw(ArgumentError("CelestialBody: flattening must be finite and in [0, 1); got $(flattening)."))
        end
        return new{T}(name, mu, equatorial_radius, flattening, naifid)
    end
end

"""
    CelestialBody(name::AbstractString, mu::Real, equatorial_radius::Real, flattening::Real, naifid::Integer)

Positional outer constructor that promotes numeric fields to a common type.
"""
function CelestialBody(
    name::AbstractString,
    mu,
    equatorial_radius,
    flattening,
    naifid::Integer,
)
    T = promote_type(typeof(mu), typeof(equatorial_radius), typeof(flattening))
    return CelestialBody{T}(String(name), T(mu), T(equatorial_radius), T(flattening), Int(naifid))
end

"""
    CelestialBody(; name="unnamed", mu=earth.mu,
                    equatorial_radius=earth.equatorial_radius,
                    flattening=earth.flattening, naifid=earth.naifid)

Keyword outer constructor that defaults all fields to Earth's values.
Numeric fields are promoted to a common element type.
"""
function CelestialBody(;
    name::AbstractString = "unnamed",
    mu::Real = EARTH_DEFAULTS.mu,
    equatorial_radius::Real = EARTH_DEFAULTS.equatorial_radius,
    flattening::Real = EARTH_DEFAULTS.flattening,
    naifid::Integer = EARTH_DEFAULTS.naifid,
)
    T = promote_type(typeof(mu), typeof(equatorial_radius), typeof(flattening))
    return CelestialBody{T}(String(name), T(mu), T(equatorial_radius), T(flattening), Int(naifid))
end

"""
    function show(io::IO, ::MIME"text/plain", body::CelestialBody)

Show method for text/plain output.
"""
function show(io::IO, ::MIME"text/plain", body::CelestialBody)
    println(io, "CelestialBody:")
    println(io, "  name               = ", body.name)
    println(io, "  μ                  = ", body.mu)
    println(io, "  Equatorial Radius  = ", body.equatorial_radius)
    println(io, "  Flattening         = ", body.flattening)
    println(io, "  NAIF ID            = ", body.naifid)
end

"""
    Base.show(io::IO, body::CelestialBody)

Delegate show to MIME"text/plain" output.
"""
function show(io::IO, body::CelestialBody)
    show(io, MIME"text/plain"(), body)
end

"""
Sun (NAIF ID 10) CelestialBody model.
"""
sun = CelestialBody("Sun", 1.32712440018e11, 696342.0, 0.0, 10)

"""
Mercury (NAIF ID 199) CelestialBody model.
"""
mercury = CelestialBody("Mercury", 22032.0, 2439.7, 0.0, 199)

"""
Venus (NAIF ID 299) CelestialBody model.
"""
venus = CelestialBody("Venus", 324858.592, 6051.8, 0.0, 299)

"""
Earth (NAIF ID 399) CelestialBody model.
"""
earth = CelestialBody(
    EARTH_DEFAULTS.name,
    EARTH_DEFAULTS.mu,
    EARTH_DEFAULTS.equatorial_radius,
    EARTH_DEFAULTS.flattening,
    EARTH_DEFAULTS.naifid,
)

"""
Moon (NAIF ID 301) CelestialBody model.
"""
moon = CelestialBody("Moon", 4902.8, 1737.4, 0.0, 301)

"""
Mars (NAIF ID 499) CelestialBody model.
"""
mars = CelestialBody("Mars", 42828.375214, 3396.2, 0.005, 499)

"""
Jupiter (NAIF ID 599) CelestialBody model.
"""
jupiter = CelestialBody("Jupiter", 126686534.0, 71492.0, 0.06487, 599)

"""
Saturn (NAIF ID 699) CelestialBody model.
"""
saturn = CelestialBody("Saturn", 37931187.0, 60268.0, 0.09796, 699)

"""
Uranus (NAIF ID 799) CelestialBody model.
"""
uranus = CelestialBody("Uranus", 5793959.0, 25559.0, 0.0229, 799)

"""
Neptune (NAIF ID 899) CelestialBody model.
"""
neptune = CelestialBody("Neptune", 6836529.0, 24764.0, 0.0171, 899)

"""
Pluto (NAIF ID 999) CelestialBody model.
""" 
pluto = CelestialBody("Pluto", 870.3, 1188.3, 0.0, 999)

"""
    function translate(from::CelestialBody, to::CelestialBody, jd_tdb::Real)

Compute the ICRF position vector from one body to another at a given TDB Julian date.

Arguments
- from: Observing/origin body.
- to: Target body.
- jd_tdb: Julian date in the TDB time scale.

# Notes:
- Requires SPICE kernels to be loaded with SPICE.furnsh before calling.
- Uses the J2000/ICRF frame; distances are kilometers.

# Returns
- 3-element position vector [x, y, z] in kilometers, from `from` to `to`, in ICRF (J2000).

# Examples
```julia
# vector from Earth to Moon
using AstroUniverse
r_em = translate(earth, moon, 2458018.0)
println(r_em)

# output (may differ slightly due to SPICE kernel versions):
3-element Vector{Real}:
 -375694.5992365016
  -96115.68241892057
  -12226.882894748915
```
"""
function translate(from::CelestialBody, to::CelestialBody, jd_tdb::Real)
    et = (jd_tdb - 2451545.0) * 86400.0
    pos, _lt = spkpos(string(to.naifid), et, "J2000", "NONE", string(from.naifid))
    return pos
end

"""
    get_gravparam(body)

Return the body’s gravitational parameter μ.
"""
@inline get_gravparam(body::CelestialBody) = body.mu

"""
    set_gravparam!(body, μ)

Set the body’s gravitational parameter μ.
"""
function set_gravparam!(body::CelestialBody, newmu::Real)
    # Validate μ (constructor invariant mirrored here)
    if !isfinite(newmu) || newmu <= 0
        throw(ArgumentError("CelestialBody: μ must be finite and > 0; got $(newmu)."))
    end
    # Preserve numeric/AD type of the field
    setfield!(body, :mu, oftype(getfield(body, :mu), newmu))
    return body
end

end