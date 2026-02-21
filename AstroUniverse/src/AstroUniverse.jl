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
export download_spice_kernel, load_spice_kernel, unload_spice_kernel, unload_all_spice_kernels
export get_spice_directory, list_cached_spice_kernels, list_downloaded_spice_kernels

"""
    ensure_kernel_download(cache_dir, filename, url)

Download SPICE kernel if it doesn't exist. Does not furnish/load the kernel.
"""
function ensure_kernel_download(cache_dir, filename, url)
    kernel_path = joinpath(cache_dir, filename)
    if !isfile(kernel_path)
        @info "Downloading SPICE kernel: $filename"
        Downloads.download(url, kernel_path)
    end
    return kernel_path
end

"""  
    download_spice_kernel(filename::AbstractString, url::AbstractString)

Download a SPICE kernel from the given URL to the cache directory.

The kernel is cached in the AstroUniverse managed scratch directory, so it will only be
downloaded once and persist across Julia sessions. This function does NOT load the kernel
into SPICE - use [`load_spice_kernel`](@ref) after downloading.

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
    load_spice_kernel(filename::AbstractString)

Load a SPICE kernel from the cache directory into the SPICE system.

The kernel file must already exist in the cache (use [`download_spice_kernel`](@ref) first
if needed). This calls SPICE.furnsh() to register the kernel.

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
You can manually place kernel files here to avoid downloading them.

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
"""
sun = CelestialBody("Sun", 1.32712440018e11, 696342.0, 0.0, 10, 
      joinpath(dirname(@__DIR__), "data", "SunTexture.jpg"))

"""
Mercury (NAIF ID 199) CelestialBody model.
"""
mercury = CelestialBody("Mercury", 22032.0, 2439.7, 0.0, 199, 
          joinpath(dirname(@__DIR__), "data", "MercuryTexture.jpg"))

"""
Venus (NAIF ID 299) CelestialBody model.
"""
venus = CelestialBody("Venus", 324858.592, 6051.8, 0.0, 299, 
        joinpath(dirname(@__DIR__), "data", "VenusTexture.jpg"))

"""
Earth (NAIF ID 399) CelestialBody model.
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
"""
moon = CelestialBody("Moon", 4902.8, 1737.4, 0.0, 301, 
       joinpath(dirname(@__DIR__), "data", "MoonTexture.jpg"))

"""
Mars (NAIF ID 499) CelestialBody model.
"""
mars = CelestialBody("Mars", 42828.375214, 3396.2, 0.005, 499, 
       joinpath(dirname(@__DIR__), "data", "MarsTexture.jpg"))

"""
Jupiter (NAIF ID 599) CelestialBody model.
"""
jupiter = CelestialBody("Jupiter", 126686534.0, 71492.0, 0.06487, 599, 
          joinpath(dirname(@__DIR__), "data", "JupiterTexture.jpg"))

"""
Saturn (NAIF ID 699) CelestialBody model.
"""
saturn = CelestialBody("Saturn", 37931187.0, 60268.0, 0.09796, 699, 
         joinpath(dirname(@__DIR__), "data", "SaturnTexture.jpg"))

"""
Uranus (NAIF ID 799) CelestialBody model.
"""
uranus = CelestialBody("Uranus", 5793959.0, 25559.0, 0.0229, 799, 
         joinpath(dirname(@__DIR__), "data", "UranusTexture.jpg"))

"""
Neptune (NAIF ID 899) CelestialBody model.
"""
neptune = CelestialBody("Neptune", 6836529.0, 24764.0, 0.0171, 899, 
          joinpath(dirname(@__DIR__), "data", "NeptuneTexture.jpg"))

"""
Pluto (NAIF ID 999) CelestialBody model.
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

"""
    function __init__()

Load SPICE kernels from managed scratch space on module initialization.
"""
function __init__()
    # Get managed cache directory for SPICE kernels
    kernel_cache = @get_scratch!("spice_kernels")
    
    # Download default kernels if needed
    ensure_kernel_download(kernel_cache, "naif0012.tls", 
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls")
    ensure_kernel_download(kernel_cache, "de440.bsp", 
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp")
    
    # Helper function to check if a kernel is already loaded
    function is_kernel_loaded(filename)
        count = ktotal("ALL")
        for i in 1:count
            result = kdata(i, "ALL")
            if result !== nothing
                file, filtyp, source, handle = result
                if basename(file) == filename
                    return true
                end
            end
        end
        return false
    end
    
    # Load default kernels only if not already loaded
    bsp_path = joinpath(kernel_cache, "de440.bsp")
    if isfile(bsp_path) && !is_kernel_loaded("de440.bsp")
        try
            furnsh(bsp_path)
            @debug "Loaded de440.bsp"
        catch e
            @warn "Failed to load de440.bsp: $e"
        end
    end
    
    tls_path = joinpath(kernel_cache, "naif0012.tls")
    if isfile(tls_path) && !is_kernel_loaded("naif0012.tls")
        try
            furnsh(tls_path)
            @debug "Loaded naif0012.tls"
        catch e
            @warn "Failed to load naif0012.tls: $e"
        end
    end
end

end