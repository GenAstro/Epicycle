# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

__precompile__()

"""
Module containing models such as celestial bodies, ephemerides, and related utilities.
"""
module AstroUniverse

using SPICE
using StaticArrays
using Scratch
using Downloads

import SHA

using EpicycleBase

import Base: show

export CelestialBody, translate, translate_state
export sun, mercury, venus, earth, moon, mars, jupiter
export saturn, uranus, neptune, pluto

export get_gravparam, set_gravparam!
export download_spice_kernel, load_spice_kernel, unload_spice_kernel, unload_all_spice_kernels
export get_spice_directory, list_cached_spice_kernels, list_downloaded_spice_kernels
export KernelSource, KERNEL_SOURCES, kernel_source

export Mu

export IAU2015Orientation, iau2015_orientation

# Body orientation models - how a body is oriented in space, set on the body.
export AbstractOrientationModel, IauPolynomialOrientation, SpiceOrientation
export body_axes_rotation, pole_axes_rotation
export set_orientation!, orientation_model, has_orientation_model
export orientation_parameters, set_orientation_parameters, set_orientation_parameters!

import EpicycleBase: AbstractParamTag, get_field, set_field!

"""
    KernelSource(filename, url, sha256, purpose)

Where a SPICE kernel comes from, and what it should contain when it arrives.

The checksum is verified once, against the bytes just downloaded, before the
file is put in place. It is not rechecked when a kernel loads: what it guards
against is a bad transfer, and rehashing a hundred megabytes on every
`using AstroUniverse` would cost far more than it catches.

# Fields
- `filename::String`: Local kernel filename.
- `url::String`: Download location.
- `sha256::String`: Expected SHA-256 checksum.
- `purpose::String`: Data supplied by the kernel.

# Example
```jldoctest
julia> using AstroUniverse

julia> kernel_source("naif0012.tls") isa KernelSource
true
```
"""
struct KernelSource
    filename::String
    url::String
    sha256::String
    purpose::String
end

"""
    KERNEL_SOURCES

Every kernel Epicycle knows how to fetch by name.

All but one are NAIF files served unchanged by JPL. Their file names carry
their own version -- a new leap second brings a `naif0013.tls` rather than an
edit to `naif0012.tls` -- so a pinned checksum stays valid.

The exception is `epicycle_de440_1950-2100.bsp`, which Gen Astro builds and
hosts. It is DE440 for the whole solar system over 1950 to 2100, and it carries
both the barycenter and the center of every planet. That matters because every
`CelestialBody` in Epicycle names a body center: JPL's own `de440.bsp` stops at
Mercury and Venus, so against it `translate(earth, mars, t)` raises rather than
returning a position.

# Example
```jldoctest
julia> using AstroUniverse

julia> first(KERNEL_SOURCES).filename
"naif0012.tls"
```
"""
const KERNEL_SOURCES = (
    KernelSource(
        "naif0012.tls",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
        "678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b",
        "Leap seconds. Required for any conversion between UTC and ephemeris time."),
    KernelSource(
        "epicycle_de440_1950-2100.bsp",
        "https://github.com/GenAstro/Epicycle/releases/download/v0.4.0/epicycle_de440_1950-2100.bsp",
        "a201d7facce002f97b65ede7684f827103681741f908da7c49ada66cd4f38506",
        "Positions of the Sun, planets and Moon, 1950 to 2100. Both the barycenter " *
        "and the center of every planet, so either NAIF id resolves."),
    KernelSource(
        "moon_pa_de440_200625.bpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc",
        "60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f",
        "Lunar orientation. Needed for the MOON_PA and MOON_ME frames; load the " *
        "matching frame kernel moon_de440_250416.tf with it."),
    KernelSource(
        "moon_de440_250416.tf",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf",
        "a47c71e9c9f33796bdafb2c9d69a7ee447b6016ecad80f71cd6f3e479f9cf768",
        "Definitions of the lunar frames. The older moon_080317.tf is not " *
        "compatible with DE440 orientation and fails naming MOON_PA_DE421."),
    KernelSource(
        "pck00011.tpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc",
        "3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1",
        "Body radii, and the IAU spin axis and prime meridian of each planet."),
)

"""
    kernel_source(filename) -> Union{KernelSource,Nothing}

The manifest entry for a kernel, or `nothing` if Epicycle does not know it.

# Returns
A [`KernelSource`](@ref) when `filename` is registered, otherwise `nothing`.

# Example
```jldoctest
julia> using AstroUniverse

julia> kernel_source("naif0012.tls").purpose isa String
true
```
"""
function kernel_source(filename::AbstractString)
    for source in KERNEL_SOURCES
        source.filename == filename && return source
    end
    return nothing
end

"""
    DEFAULT_KERNELS

The kernels fetched and loaded when AstroUniverse starts.

Leap seconds and solar system positions, because nothing works without them,
and lunar orientation, because `MoonPA` and `MoonME` are frames Epicycle
offers and each has to work when it is named, without a separate download.

Together these are about 106 MB, downloaded once and cached across sessions.
"""
const DEFAULT_KERNELS = ("naif0012.tls",
                         "epicycle_de440_1950-2100.bsp",
                         "moon_pa_de440_200625.bpc",
                         "moon_de440_250416.tf")

"""
    ensure_kernel_download(cache_dir, filename, url; sha256 = nothing)

Download a SPICE kernel unless it is already cached. Does not load the kernel.

Pass `sha256` to have the downloaded bytes checked before they are cached.
"""
function ensure_kernel_download(cache_dir, filename, url; sha256 = nothing)
    kernel_path = joinpath(cache_dir, filename)
    isfile(kernel_path) && return kernel_path

    @info "Downloading SPICE kernel: $filename"

    # Download under a temporary name and move it into place only once it is
    # whole. An interrupted transfer leaves a partial file, and `isfile` cannot
    # tell a partial kernel from a good one -- so a partial must never be
    # allowed to occupy the real name, where it would be trusted forever.
    partial_path = kernel_path * ".part"
    try
        Downloads.download(url, partial_path)

        if sha256 !== nothing
            found = bytes2hex(open(SHA.sha256, partial_path))
            if found != lowercase(sha256)
                error("""
                      Checksum mismatch for $(filename); nothing was cached.
                          expected  $(lowercase(sha256))
                          received  $(found)
                      The download completed but did not deliver the expected file.
                      Source: $(url)
                      """)
            end
        end

        mv(partial_path, kernel_path; force = true)
    finally
        isfile(partial_path) && rm(partial_path; force = true)
    end

    return kernel_path
end

"""  
    download_spice_kernel(filename::AbstractString, url::AbstractString)

Download a SPICE kernel from the given URL to the cache directory.

The kernel is cached in the AstroUniverse managed scratch directory, so it will only be
downloaded once and persist across Julia sessions. This function does NOT load the kernel
into SPICE - use [`load_spice_kernel`](@ref) after downloading.

# Returns
`nothing` after the kernel is present in the cache.

# Arguments
- `filename`: Name of the kernel file (e.g., "de441.bsp")
- `url`: Full URL to download the kernel from

# Examples
```julia
using AstroUniverse

# Download Mars satellite ephemeris (Phobos, Deimos)
download_spice_kernel("mar099.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar099.bsp")

# Then load it
load_spice_kernel("mar099.bsp")
```

All kernels available at: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/

See also: [`load_spice_kernel`](@ref), [`get_spice_directory`](@ref)
"""
function download_spice_kernel(filename::AbstractString, url::AbstractString)
    cache_dir = @get_scratch!("spice_kernels")
    ensure_kernel_download(cache_dir, String(filename), String(url))
    return nothing
end

"""
    download_spice_kernel(filename::AbstractString)

Download one of the kernels Epicycle knows about, by name, and check it.

The download location and expected checksum come from [`KERNEL_SOURCES`](@ref), so
no URL is required and a bad transfer is caught rather than cached. Use
the two-argument form for a kernel that is not in the manifest.

The kernel is cached and not loaded. Call [`load_spice_kernel`](@ref) to use it.

# Returns
`nothing` after the kernel is present in the cache.

# Examples

```julia
# The lunar frames need orientation data and the frames that go with it.
download_spice_kernel("moon_pa_de440_200625.bpc")
download_spice_kernel("moon_de440_250416.tf")
load_spice_kernel("moon_pa_de440_200625.bpc")
load_spice_kernel("moon_de440_250416.tf")
```

See also: [`load_spice_kernel`](@ref), [`kernel_source`](@ref)
"""
function download_spice_kernel(filename::AbstractString)
    source = kernel_source(filename)
    if source === nothing
        known = join(("  " * k.filename for k in KERNEL_SOURCES), "
")
        error("""
              Epicycle has no manifest entry for $(filename).
              Kernels it knows by name:
              $(known)
              For anything else, give the URL too: download_spice_kernel(name, url)
              """)
    end
    cache_dir = @get_scratch!("spice_kernels")
    ensure_kernel_download(cache_dir, source.filename, source.url; sha256 = source.sha256)
    return nothing
end

"""  
    load_spice_kernel(filename::AbstractString)

Load a SPICE kernel from the cache directory into the SPICE system.

The kernel file must already exist in the cache (use [`download_spice_kernel`](@ref) first
if needed). This calls SPICE.furnsh() to register the kernel.

# Returns
`nothing` after the kernel is loaded or found already loaded.

# Arguments
- `filename`: Name of the kernel file in the cache (e.g., "de441.bsp")

# Examples
```julia
using AstroUniverse

# Download then load
download_spice_kernel("de440.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp")
load_spice_kernel("de440.bsp")

# Or just load if already downloaded
load_spice_kernel("naif0012.tls")  # Default kernel
```

See also: [`download_spice_kernel`](@ref), [`unload_spice_kernel`](@ref), [`unload_all_spice_kernels`](@ref)
"""
function load_spice_kernel(filename::AbstractString)
    cache_dir = @get_scratch!("spice_kernels")
    kernel_path = joinpath(cache_dir, String(filename))
    
    if !isfile(kernel_path)
        error("Kernel file not found in cache: $filename. Use download_spice_kernel() first.")
    end
    
    # Check if already loaded (SPICE uses 1-based indexing)
    count = ktotal("ALL")
    for i in 1:count
        result = kdata(i, "ALL")
        if result !== nothing
            file, filtyp, source, handle = result
            if basename(file) == filename
                @warn "Kernel '$filename' is already loaded. Skipping duplicate load."
                return nothing
            end
        end
    end
    
    furnsh(kernel_path)
    return nothing
end

"""  
    unload_spice_kernel(filename::AbstractString)

Unload a specific SPICE kernel from the SPICE system.

This calls SPICE.unload() to remove the kernel from memory. The file remains in the cache.

# Returns
`nothing` after every loaded instance is removed.

# Arguments
- `filename`: Name of the kernel file to unload (e.g., "de440.bsp")

# Examples
```julia
using AstroUniverse

# Swap planetary ephemeris versions
unload_spice_kernel("de440.bsp")
load_spice_kernel("de440.bsp")
```

See also: [`load_spice_kernel`](@ref), [`unload_all_spice_kernels`](@ref)
"""
function unload_spice_kernel(filename::AbstractString)
    cache_dir = @get_scratch!("spice_kernels")
    kernel_path = joinpath(cache_dir, String(filename))
    
    if !isfile(kernel_path)
        error("Kernel file not found in cache: $filename. Cannot unload.")
    end
    
    # Check if kernel is actually loaded (SPICE uses 1-based indexing)
    is_loaded = false
    count = ktotal("ALL")
    for i in 1:count
        result = kdata(i, "ALL")
        if result !== nothing
            file, filtyp, source, handle = result
            if basename(file) == filename
                is_loaded = true
                break
            end
        end
    end
    
    if !is_loaded
        @warn "Kernel '$filename' is not currently loaded. Nothing to unload."
        return nothing
    end
    
    # Unload all instances of this kernel
    unload_count = 0
    while true
        # Check if still loaded (SPICE uses 1-based indexing)
        still_loaded = false
        count = ktotal("ALL")
        for i in 1:count
            result = kdata(i, "ALL")
            if result !== nothing
                file, filtyp, source, handle = result
                if basename(file) == filename
                    still_loaded = true
                    break
                end
            end
        end
        
        if !still_loaded
            break
        end
        
        unload(kernel_path)
        unload_count += 1
    end
    
    if unload_count > 1
        @warn "Kernel '$filename' was loaded $unload_count times. All instances have been unloaded."
    end
    
    return nothing
end

"""  
    unload_all_spice_kernels()

Unload all SPICE kernels from memory.

This calls SPICE.kclear() to remove all loaded kernels. Useful for creating custom kernel
configurations. The files remain in the storage directory.

# Returns
`nothing` after the SPICE kernel pool is cleared.

# Examples
```julia
using AstroUniverse

# Create custom configuration
unload_all_spice_kernels()
list_cached_spice_kernels()  # Shows none loaded
```

See also: [`load_spice_kernel`](@ref), [`unload_spice_kernel`](@ref)
"""
function unload_all_spice_kernels()
    kclear()
    return nothing
end

"""  
    get_spice_directory()

Return the path to the AstroUniverse SPICE kernel cache directory.

This directory persists across Julia sessions and is managed by Scratch.jl.
Kernel files may also be placed here manually.

# Returns
The absolute path to the managed kernel directory.

# Examples
```julia
using AstroUniverse

dir = get_spice_directory()
println("SPICE kernels cached at: ", dir)

# Manually copy a kernel file to the cache
# cp("my_kernel.bsp", joinpath(dir, "my_kernel.bsp"))
# Then load it
# load_spice_kernel("my_kernel.bsp")
```

See also: [`download_spice_kernel`](@ref), [`list_downloaded_spice_kernels`](@ref)
"""
get_spice_directory() = @get_scratch!("spice_kernels")

"""
    list_downloaded_spice_kernels()

List all SPICE kernel files (.bsp and .tls) downloaded to the persistent storage directory.

Displays kernel filenames with their file sizes in a formatted table. These are files available
for loading, not necessarily currently loaded in SPICE memory. Use [`list_cached_spice_kernels`](@ref)
to see which kernels are currently loaded.

# Returns
`nothing` after printing the downloaded kernels.

# Examples
```julia
using AstroUniverse

list_downloaded_spice_kernels()

# output
Downloaded SPICE Kernels:
  naif0012.tls              (5.3 KB)
  de440.bsp                 (114.0 MB)
  mar099.bsp                (2.1 MB)

Storage location: /path/to/scratch/spice_kernels
```

See also: [`download_spice_kernel`](@ref), [`list_cached_spice_kernels`](@ref), [`get_spice_directory`](@ref)
"""
function list_downloaded_spice_kernels()
    cache_dir = get_spice_directory()
    
    # Find all .bsp and .tls files
    kernel_files = filter(readdir(cache_dir)) do f
        endswith(lowercase(f), ".bsp") || endswith(lowercase(f), ".tls")
    end
    
    if isempty(kernel_files)
        println("No SPICE kernels found in storage directory.")
        println("Storage location: ", cache_dir)
        return
    end
    
    println("Downloaded SPICE Kernels:")
    
    # Sort files: .tls first, then .bsp alphabetically
    sort!(kernel_files, by = f -> (!endswith(lowercase(f), ".tls"), lowercase(f)))
    
    for file in kernel_files
        filepath = joinpath(cache_dir, file)
        size_bytes = filesize(filepath)
        
        # Format file size nicely
        if size_bytes < 1024
            size_str = string(size_bytes, " B")
        elseif size_bytes < 1024^2
            size_str = string(round(size_bytes / 1024, digits=1), " KB")
        elseif size_bytes < 1024^3
            size_str = string(round(size_bytes / 1024^2, digits=1), " MB")
        else
            size_str = string(round(size_bytes / 1024^3, digits=1), " GB")
        end
        
        println("  ", rpad(file, 25), " (", size_str, ")")
    end
    
    println("\nStorage location: ", cache_dir)
end

"""
    list_cached_spice_kernels()

List all SPICE kernels currently loaded in memory.

Displays the filenames of kernels that have been furnished to SPICE and are actively being used.
This shows what's actually in the SPICE kernel pool, not what's downloaded to disk.

# Returns
`nothing` after printing the loaded kernels.

# Examples
```julia
using AstroUniverse

list_cached_spice_kernels()

# output
Cached SPICE Kernels (2 loaded):
  naif0012.tls
  de440.bsp
```

See also: [`load_spice_kernel`](@ref), [`unload_spice_kernel`](@ref), [`list_downloaded_spice_kernels`](@ref)
"""
function list_cached_spice_kernels()
    count = ktotal("ALL")
    
    if count == 0
        println("No SPICE kernels currently loaded.")
        return
    end
    
    # Collect kernel names (SPICE uses 1-based indexing like Fortran)
    kernel_names = String[]
    for i in 1:count
        result = kdata(i, "ALL")
        if result !== nothing
            file, filtyp, source, handle = result
            push!(kernel_names, basename(file))
        end
    end
    
    println("Cached SPICE Kernels ($(length(kernel_names)) loaded):")
    
    # Display kernel names
    for name in kernel_names
        println("  ", name)
    end
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
- texture_file::String — path to texture image file for visualization (optional)

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
                                 naifid=301,
                                 texture_file="path/to/moon_texture.jpg");
show(moon_like)

# output
CelestialBody:
  name               = MyMoon
  μ                  = 4902.8
  Equatorial Radius  = 1737.4
  Flattening         = 0.0
  NAIF ID            = 301
  Texture File       = path/to/moon_texture.jpg
```
"""
mutable struct CelestialBody{T<:Real} <: AbstractPoint
    name::String
    mu::T
    equatorial_radius::T
    flattening::T
    naifid::Int
    texture_file::String

    function CelestialBody{T}(
        name::String,
        mu::T,
        equatorial_radius::T,
        flattening::T,
        naifid::Int,
        texture_file::String,
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
        return new{T}(name, mu, equatorial_radius, flattening, naifid, texture_file)
    end
end

"""
   CelestialBody(name::AbstractString, mu::Real, equatorial_radius::Real, flattening::Real, naifid::Integer, texture_file::AbstractString)

Positional outer constructor that promotes numeric fields to a common type.
"""
function CelestialBody(
    name::AbstractString,
    mu,
    equatorial_radius,
    flattening,
    naifid::Integer,
    texture_file::AbstractString="",
)
    T = promote_type(typeof(mu), typeof(equatorial_radius), typeof(flattening))
    return CelestialBody{T}(String(name), T(mu), T(equatorial_radius), T(flattening), Int(naifid), String(texture_file))
end

"""
    CelestialBody(; name="unnamed", mu=earth.mu,
                    equatorial_radius=earth.equatorial_radius,
                    flattening=earth.flattening, naifid=earth.naifid,
                    texture_file="")

Keyword outer constructor that defaults all fields to Earth's values.
Numeric fields are promoted to a common element type.
"""
function CelestialBody(;
    name::AbstractString = "unnamed",
    mu::Real = EARTH_DEFAULTS.mu,
    equatorial_radius::Real = EARTH_DEFAULTS.equatorial_radius,
    flattening::Real = EARTH_DEFAULTS.flattening,
    naifid::Integer = EARTH_DEFAULTS.naifid,
    texture_file::AbstractString = "",
)
    T = promote_type(typeof(mu), typeof(equatorial_radius), typeof(flattening))
    return CelestialBody{T}(String(name), T(mu), T(equatorial_radius), T(flattening),
           Int(naifid), String(texture_file))
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
    texture_display = isempty(body.texture_file) ? "(none)" : body.texture_file
    println(io, "  Texture File       = ", texture_display)
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

# Example
```jldoctest
julia> using AstroUniverse; sun.naifid
10
```
"""
sun = CelestialBody("Sun", 1.32712440018e11, 696342.0, 0.0, 10, 
      joinpath(dirname(@__DIR__), "data", "SunTexture.jpg"))

"""
Mercury (NAIF ID 199) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; mercury.naifid
199
```
"""
mercury = CelestialBody("Mercury", 22032.0, 2439.7, 0.0, 199, 
          joinpath(dirname(@__DIR__), "data", "MercuryTexture.jpg"))

"""
Venus (NAIF ID 299) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; venus.naifid
299
```
"""
venus = CelestialBody("Venus", 324858.592, 6051.8, 0.0, 299, 
        joinpath(dirname(@__DIR__), "data", "VenusTexture.jpg"))

"""
Earth (NAIF ID 399) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; earth.naifid
399
```
"""
earth = CelestialBody(
    EARTH_DEFAULTS.name,
    EARTH_DEFAULTS.mu,
    EARTH_DEFAULTS.equatorial_radius,
    EARTH_DEFAULTS.flattening,
    EARTH_DEFAULTS.naifid,
    joinpath(dirname(@__DIR__), "data", "EarthTexture.jpg"),
)

"""
Moon (NAIF ID 301) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; moon.naifid
301
```
"""
moon = CelestialBody("Moon", 4902.8, 1737.4, 0.0, 301, 
       joinpath(dirname(@__DIR__), "data", "MoonTexture.jpg"))

"""
Mars (NAIF ID 499) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; mars.naifid
499
```
"""
mars = CelestialBody("Mars", 42828.375214, 3396.2, 0.005, 499, 
       joinpath(dirname(@__DIR__), "data", "MarsTexture.jpg"))

"""
Jupiter (NAIF ID 599) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; jupiter.naifid
599
```
"""
jupiter = CelestialBody("Jupiter", 126686534.0, 71492.0, 0.06487, 599, 
          joinpath(dirname(@__DIR__), "data", "JupiterTexture.jpg"))

"""
Saturn (NAIF ID 699) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; saturn.naifid
699
```
"""
saturn = CelestialBody("Saturn", 37931187.0, 60268.0, 0.09796, 699, 
         joinpath(dirname(@__DIR__), "data", "SaturnTexture.jpg"))

"""
Uranus (NAIF ID 799) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; uranus.naifid
799
```
"""
uranus = CelestialBody("Uranus", 5793959.0, 25559.0, 0.0229, 799, 
         joinpath(dirname(@__DIR__), "data", "UranusTexture.jpg"))

"""
Neptune (NAIF ID 899) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; neptune.naifid
899
```
"""
neptune = CelestialBody("Neptune", 6836529.0, 24764.0, 0.0171, 899, 
          joinpath(dirname(@__DIR__), "data", "NeptuneTexture.jpg"))

"""
Pluto (NAIF ID 999) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; pluto.naifid
999
```
""" 
pluto = CelestialBody("Pluto", 870.3, 1188.3, 0.0, 999, "")

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
    translate_state(from::CelestialBody, to::CelestialBody, jd_tdb::Real)

Position **and velocity** of `to` relative to `from` at TDB Julian date
`jd_tdb`, in ICRF.

# Arguments
- `from`: Observing/origin body.
- `to`: Target body.
- `jd_tdb`: Julian date in the TDB time scale.

# Returns
- 6-element state `[x, y, z, vx, vy, vz]` in km and km/s, from `from` to `to`.

# Notes
- The name differs from [`translate`](@ref) because the *return shape*
  differs: one is a position, the other a state. Overloading a single name to
  return three or six elements depending on the call would make the shape
  depend on something the reader has to infer.
- SPICE labels this frame `"J2000"`, which is **ICRF-aligned** and is not the
  FK5 mean equator and equinox of J2000. The two differ by the ~23 mas frame
  bias. Callers should treat the result as ICRF.
- Requires SPICE kernels to be loaded.
- **Not differentiable.** SPICE is a C library, so an origin shift built on
  this sits outside any AD path.

# Examples
```julia
using AstroUniverse
s_em = translate_state(earth, moon, 2458018.0)
```
"""
function translate_state(from::CelestialBody, to::CelestialBody, jd_tdb::Real)
    et = (jd_tdb - 2451545.0) * 86400.0
    state, _lt = spkezr(string(to.naifid), et, "J2000", "NONE", string(from.naifid))
    return state
end

"""
    get_gravparam(body)

Return the body’s gravitational parameter μ.

# Returns
The gravitational parameter in km^3/s^2.

# Example
```jldoctest
julia> using AstroUniverse; get_gravparam(earth) == earth.mu
true
```
"""
@inline get_gravparam(body::CelestialBody) = body.mu

"""
    set_gravparam!(body, μ)

Set the body’s gravitational parameter μ.

# Returns
The modified body.

# Example
```jldoctest
julia> using AstroUniverse

julia> body = CelestialBody("Test", 1.0, 1.0, 0.0, 1000001);

julia> set_gravparam!(body, 2.0) === body
true
```
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

"""
    function __init__()

Load SPICE kernels from managed scratch space on module initialization.
"""
# Defined at module scope rather than inside `__init__`. A function written inside `__init__`
# gets a fresh type every time that body is compiled, so it cannot come from the precompile
# cache and is rebuilt at every startup. `--trace-compile` reported
# `precompile(Tuple{typeof(AstroUniverse.__init__)}) # recompile` on every run, and
# `@time_imports` charged this module 102 ms at 100% recompilation.
"""Whether SPICE already has a kernel of this basename loaded."""
function is_kernel_loaded(filename)
    for i in 1:ktotal("ALL")
        result = kdata(i, "ALL")
        result === nothing && continue
        basename(first(result)) == filename && return true
    end
    return false
end

function __init__()
    # Get managed cache directory for SPICE kernels
    kernel_cache = @get_scratch!("spice_kernels")

    # Download default kernels if needed. Each returns immediately when the file is already
    # cached, so a normal start does no network I/O.
    for filename in DEFAULT_KERNELS
        source = kernel_source(filename)
        ensure_kernel_download(kernel_cache, source.filename, source.url;
                               sha256 = source.sha256)
    end

    # Load default kernels only if not already loaded
    for filename in DEFAULT_KERNELS
        path = joinpath(kernel_cache, filename)
        (isfile(path) && !is_kernel_loaded(filename)) || continue
        try
            furnsh(path)
            @debug "Loaded $filename"
        catch e
            @warn "Failed to load $filename: $e"
        end
    end
end

# ---------------------------------------------------------------------------
# Tag / Variable System — Mu
# ---------------------------------------------------------------------------

"""
    Mu <: AbstractParamTag

Tag identifying the gravitational parameter μ = GM on a `CelestialBody`.

# Example
```jldoctest
julia> using AstroUniverse; Mu() isa Mu
true
```
"""
struct Mu <: AbstractParamTag end

"""Return the gravitational parameter of `body` (km³/s²)."""
get_field(body::CelestialBody, ::Mu) = body.mu

"""Set the gravitational parameter of `body` to `v` (km³/s²)."""
function set_field!(body::CelestialBody, ::Mu, v::Real)
    body.mu = convert(typeof(body.mu), v)
    return nothing
end

# ---------------------------------------------------------------------------
# Earth Orientation Parameters (global service)
# ---------------------------------------------------------------------------

include("eop.jl")

# ---------------------------------------------------------------------------
# IAU 2015 planetary orientation (pole + prime meridian + rates)
# ---------------------------------------------------------------------------

include("iau_2015_orientation.jl")
include("orientation_models.jl")

end
