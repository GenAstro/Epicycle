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
using AstroEpochs: AstroEpochs, Time

import Base: show

export CelestialBody, translate, translate_state
export sun, mercury, venus, earth, moon, mars, jupiter
export saturn, uranus, neptune, pluto

export get_gravparam, set_gravparam!
export download_spice_kernel, load_spice_kernel, unload_spice_kernel, unload_all_spice_kernels
export get_spice_directory, loaded_spice_kernels, cached_spice_kernels
export KernelSource, KERNEL_SOURCES, kernel_source

export Mu

export IAU2015Orientation, iau2015_orientation

# Body orientation models - how a body is oriented in space, set on the body.
export AbstractOrientationModel, IAU2015, IAU1991, LunarPA, LunarME, SpiceOrientation
export body_axes_rotation, pole_axes_rotation
export set_orientation!, orientation_model, has_orientation_model
export orientation_parameters, set_orientation_parameters, set_orientation_parameters!

import EpicycleBase: AbstractParamTag, get_field, set_field!

"""
    KernelSource(filename, url, sha256, size, purpose)

Where a SPICE kernel comes from, and what it should contain when it arrives.

The checksum is verified once, against the bytes just downloaded, before the
file is put in place. Rehashing a hundred megabytes on every `using AstroUniverse`
would cost far more than it catches, so a cached file has its size checked instead,
which finds a truncated copy for almost no cost.

# Fields
- `filename::String`: Local kernel filename.
- `url::String`: Download location.
- `sha256::String`: Expected SHA-256 checksum.
- `size::Int`: Expected file size in bytes.
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
    size::Int
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
        5_257,
        "Leap seconds. Required for any conversion between UTC and ephemeris time."),
    KernelSource(
        "epicycle_de440_1950-2100.bsp",
        "https://github.com/GenAstro/Epicycle/releases/download/v0.4.0/epicycle_de440_1950-2100.bsp",
        "a201d7facce002f97b65ede7684f827103681741f908da7c49ada66cd4f38506",
        98_342_912,
        "Positions of the Sun, planets and Moon, 1950 to 2100. Both the barycenter " *
        "and the center of every planet, so either NAIF id resolves."),
    KernelSource(
        "moon_pa_de440_200625.bpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc",
        "60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f",
        12_863_488,
        "Lunar orientation. Needed for the MOON_PA and MOON_ME frames; load the " *
        "matching frame kernel moon_de440_250416.tf with it."),
    KernelSource(
        "moon_de440_250416.tf",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_250416.tf",
        "a47c71e9c9f33796bdafb2c9d69a7ee447b6016ecad80f71cd6f3e479f9cf768",
        19_478,
        "Definitions of the lunar frames. The older moon_080317.tf is not " *
        "compatible with DE440 orientation and fails naming MOON_PA_DE421."),
    KernelSource(
        "pck00011.tpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc",
        "3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1",
        131_226,
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

The kernels Epicycle needs for ephemeris and lunar orientation.

Leap seconds and solar system positions, because nothing works without them,
and lunar orientation, because `MoonPA` and `MoonME` are frames Epicycle
offers and each has to work when it is named, without a separate download.

Together these are about 111 MB. They are downloaded at the first call that needs
them and cached across sessions, and those already cached load when AstroUniverse
starts, so `using AstroUniverse` never reaches the network.
"""
const DEFAULT_KERNELS = ("naif0012.tls",
                         "epicycle_de440_1950-2100.bsp",
                         "moon_pa_de440_200625.bpc",
                         "moon_de440_250416.tf")

"""
    ensure_kernel_download(cache_dir, filename, url; sha256 = nothing, size = nothing)

Download a SPICE kernel unless it is already cached. Does not load the kernel.

Pass `sha256` to have the downloaded bytes checked before they are cached, and
`size` to have a cached file of any other size replaced.
"""
function ensure_kernel_download(cache_dir, filename, url; sha256 = nothing, size = nothing)
    kernel_path = joinpath(cache_dir, filename)
    good(path) = isfile(path) && (size === nothing || filesize(path) == size)
    good(kernel_path) && return kernel_path

    if isfile(kernel_path)
        @warn "The cached SPICE kernel $(filename) is $(filesize(kernel_path)) bytes and should " *
              "be $(size). Downloading it again."
    else
        @info "Downloading SPICE kernel: $filename"
    end

    # Download under a temporary name and move it into place only once it is
    # whole, so a partial file never occupies the real name. Each attempt has its own
    # name, so processes fetching the same kernel at once never share a file.
    partial_path = tempname(cache_dir; cleanup = false)
    try
        Downloads.download(url, partial_path; timeout = _KERNEL_DOWNLOAD_TIMEOUT)

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

        try
            _replace_file(partial_path, kernel_path)
        catch
            # Windows refuses to replace a file another process has open. If that process
            # got a good copy here first, this one is finished too.
            good(kernel_path) || rethrow()
        end
    finally
        rm(partial_path; force = true)
    end

    return kernel_path
end

const _KERNEL_DOWNLOAD_TIMEOUT = 3600.0     # seconds

"""  
    download_spice_kernel(filename::AbstractString, url::AbstractString)

Download a SPICE kernel from the given URL to the cache directory.

The kernel is cached in the AstroUniverse managed scratch directory, so it will only be
downloaded once and persist across Julia sessions. A file of that name already in the
cache is kept; delete it to fetch a new one. This function does NOT load the kernel
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
    # A kernel the manifest knows, from the place the manifest knows, is checked as one.
    source = kernel_source(filename)
    known = source !== nothing && source.url == url
    ensure_kernel_download(cache_dir, String(filename), String(url);
                           sha256 = known ? source.sha256 : nothing,
                           size = known ? source.size : nothing)
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
    ensure_kernel_download(cache_dir, source.filename, source.url;
                           sha256 = source.sha256, size = source.size)
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
    
    @lock _SPICE_LOCK begin
        # A different file that shares the name is a different kernel, so only this file counts.
        if any(path -> samefile(path, kernel_path), _loaded_paths(filename))
            @warn "Kernel '$filename' is already loaded. Skipping duplicate load."
            return nothing
        end
        furnsh(kernel_path)
    end
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
    @lock _SPICE_LOCK begin
        # Unloaded by the path SPICE reports, which is where the file was loaded from and
        # need not be the cache.
        paths = _loaded_paths(filename)
        if isempty(paths)
            isfile(joinpath(get_spice_directory(), String(filename))) ||
                error("Kernel file not found in cache: $filename. Cannot unload.")
            @warn "Kernel '$filename' is not currently loaded. Nothing to unload."
            return nothing
        end
        foreach(unload, paths)
        length(paths) > 1 &&
            @warn "Kernel '$filename' was loaded $(length(paths)) times. All instances have been unloaded."
    end
    return nothing
end

# The path of every loaded kernel with this file name, once for each time it was loaded.
function _loaded_paths(filename)
    paths = String[]
    for i in 1:ktotal("ALL")
        result = kdata(i, "ALL")
        result === nothing && continue
        basename(first(result)) == filename && push!(paths, first(result))
    end
    return paths
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
loaded_spice_kernels()  # Empty
```

After this call the default kernels are no longer fetched or loaded for you, so the
configuration you load is the one in force.

See also: [`load_spice_kernel`](@ref), [`unload_spice_kernel`](@ref)
"""
function unload_all_spice_kernels()
    @lock _SPICE_LOCK begin
        kclear()
        _KERNELS_READY[] = true
    end
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

See also: [`download_spice_kernel`](@ref), [`cached_spice_kernels`](@ref)
"""
get_spice_directory() = @get_scratch!("spice_kernels")

"""
    cached_spice_kernels() -> Vector{String}

The SPICE kernel files in the cache directory, sorted by name. These are the files available
to load; [`loaded_spice_kernels`](@ref) gives the ones in use.

# Returns
The file names, of every kernel type.

# Examples
```julia
using AstroUniverse

cached_spice_kernels()
```

See also: [`download_spice_kernel`](@ref), [`loaded_spice_kernels`](@ref), [`get_spice_directory`](@ref)
"""
cached_spice_kernels() = sort!(filter(_is_kernel_file, readdir(get_spice_directory())); by = lowercase)

_is_kernel_file(name) = lowercase(last(splitext(name))) in
    (".bsp", ".bpc", ".tls", ".tpc", ".tf", ".bc", ".tsc", ".ti", ".tm", ".dsk")

"""
    loaded_spice_kernels() -> Vector{String}

The SPICE kernels loaded in this session, in the order they were loaded. A kernel loaded
twice appears twice. [`cached_spice_kernels`](@ref) gives the files on disk.

# Returns
The file names.

# Examples
```julia
using AstroUniverse

"naif0012.tls" in loaded_spice_kernels()
```

See also: [`load_spice_kernel`](@ref), [`unload_spice_kernel`](@ref), [`cached_spice_kernels`](@ref)
"""
function loaded_spice_kernels()
    @lock _SPICE_LOCK begin
        names = String[]
        for i in 1:ktotal("ALL")
            result = kdata(i, "ALL")
            result === nothing || push!(names, basename(first(result)))
        end
        return names
    end
end

"""
    CelestialBody(name, mu, equatorial_radius, flattening, naifid, texture_file = "")
    CelestialBody(; name, mu, equatorial_radius, naifid, flattening = 0.0, texture_file = "")

A celestial body and its physical constants.

# Fields
- `name::String`
- `mu::T`: gravitational parameter, km³/s².
- `equatorial_radius::T`: km.
- `flattening::T`: geometric flattening, dimensionless, in [0, 1).
- `naifid::Int`: NAIF body ID, which selects the body's ephemeris and orientation.
- `texture_file::String`: path to a texture image for visualization, or empty.

# Notes
- The numeric fields share one element type `T`, at least `Float64`. A `BigFloat` or a dual
  number in any of them promotes the others.
- Every value is checked when it is set, in the constructor or afterward.

# Examples
```julia
using AstroUniverse
phobos = CelestialBody(name = "Phobos", mu = 7.087e-4, equatorial_radius = 13.0, naifid = 401)

# output
CelestialBody("Phobos", naifid = 401)
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
        _check_body_value(Val(:mu), mu)
        _check_body_value(Val(:equatorial_radius), equatorial_radius)
        _check_body_value(Val(:flattening), flattening)
        return new{T}(name, mu, equatorial_radius, flattening, naifid, texture_file)
    end
end

# The one place a body's values are checked: the constructors and every way of setting a value
# after construction come through here.
function _check_body_value(::Val{:mu}, v)
    isfinite(v) && v > 0 ||
        throw(ArgumentError("CelestialBody: μ must be finite and > 0; got $(v)."))
    return v
end
function _check_body_value(::Val{:equatorial_radius}, v)
    isfinite(v) && v > 0 ||
        throw(ArgumentError("CelestialBody: equatorial_radius must be finite and > 0; got $(v)."))
    return v
end
function _check_body_value(::Val{:flattening}, v)
    isfinite(v) && 0 <= v < 1 ||
        throw(ArgumentError("CelestialBody: flattening must be finite and in [0, 1); got $(v)."))
    return v
end

function Base.setproperty!(body::CelestialBody, name::Symbol, v)
    x = convert(fieldtype(typeof(body), name), v)
    name in (:mu, :equatorial_radius, :flattening) && _check_body_value(Val(name), x)
    return setfield!(body, name, x)
end

function CelestialBody(
    name::AbstractString,
    mu::Real,
    equatorial_radius::Real,
    flattening::Real,
    naifid::Integer,
    texture_file::AbstractString = "",
)
    T = promote_type(Float64, typeof(mu), typeof(equatorial_radius), typeof(flattening))
    return CelestialBody{T}(String(name), T(mu), T(equatorial_radius), T(flattening), Int(naifid),
                            String(texture_file))
end

function CelestialBody(;
    name::AbstractString,
    mu::Real,
    equatorial_radius::Real,
    naifid::Integer,
    flattening::Real = 0.0,
    texture_file::AbstractString = "",
)
    return CelestialBody(name, mu, equatorial_radius, flattening, naifid, texture_file)
end

function show(io::IO, ::MIME"text/plain", body::CelestialBody)
    println(io, "CelestialBody:")
    println(io, "  name               = ", body.name)
    println(io, "  μ                  = ", body.mu, " km³/s²")
    println(io, "  Equatorial Radius  = ", body.equatorial_radius, " km")
    println(io, "  Flattening         = ", body.flattening)
    println(io, "  NAIF ID            = ", body.naifid)
    print(io,   "  Texture File       = ", isempty(body.texture_file) ? "(none)" : body.texture_file)
end

show(io::IO, body::CelestialBody) = print(io, "CelestialBody(", repr(body.name), ", naifid = ", body.naifid, ")")

# ---------------------------------------------------------------------------
# The built-in bodies.
#
# Set up for Earth orbits, the common case. Earth carries WGS84: its μ is what the Earth
# gravity models and geodetic work use. The others are third bodies seen from far away, where a
# planet and its moons act as one mass, so each carries its system μ from NAIF's gm_de440.tpc
# (the barycenter IDs 4 to 9, with the planet's own ID where it has no moons). Radii and
# flattening are NAIF's pck00011.tpc, with flattening computed from its radii.
#
# Work centered on another planet, or in barycentric coordinates, changes μ with
# `set_gravparam!`; the documentation says what to change.
# ---------------------------------------------------------------------------

"""
Sun (NAIF ID 10) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; sun.naifid
10
```
"""
const sun = CelestialBody("Sun", 1.3271244004127942e11, 695700.0, 0.0, 10)

"""
Mercury (NAIF ID 199) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; mercury.naifid
199
```
"""
const mercury = CelestialBody("Mercury", 2.2031868551400003e4, 2440.53, 1 - 2438.26 / 2440.53, 199)

"""
Venus (NAIF ID 299) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; venus.naifid
299
```
"""
const venus = CelestialBody("Venus", 3.24858592e5, 6051.8, 0.0, 299)

"""
Earth (NAIF ID 399) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; earth.naifid
399
```
"""
const earth = CelestialBody("Earth", 398600.4418, 6378.137, 1 / 298.257223563, 399)

"""
Moon (NAIF ID 301) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; moon.naifid
301
```
"""
const moon = CelestialBody("Moon", 4.9028001184575496e3, 1737.4, 0.0, 301)

"""
Mars (NAIF ID 499) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; mars.naifid
499
```
"""
const mars = CelestialBody("Mars", 4.2828375815756102e4, 3396.19, 1 - 3376.20 / 3396.19, 499)

"""
Jupiter (NAIF ID 599) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; jupiter.naifid
599
```
"""
const jupiter = CelestialBody("Jupiter", 1.2671276409999998e8, 71492.0, 1 - 66854.0 / 71492.0, 599)

"""
Saturn (NAIF ID 699) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; saturn.naifid
699
```
"""
const saturn = CelestialBody("Saturn", 3.7940584841799997e7, 60268.0, 1 - 54364.0 / 60268.0, 699)

"""
Uranus (NAIF ID 799) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; uranus.naifid
799
```
"""
const uranus = CelestialBody("Uranus", 5.7945563999999985e6, 25559.0, 1 - 24973.0 / 25559.0, 799)

"""
Neptune (NAIF ID 899) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; neptune.naifid
899
```
"""
const neptune = CelestialBody("Neptune", 6.8365271005803989e6, 24764.0, 1 - 24341.0 / 24764.0, 899)

"""
Pluto (NAIF ID 999) CelestialBody model.

# Example
```jldoctest
julia> using AstroUniverse; pluto.naifid
999
```
""" 
const pluto = CelestialBody("Pluto", 9.755e2, 1188.3, 0.0, 999)

# Texture images ship in the package, and where the package is can differ between the machine
# that precompiled it and the one running it, so the paths are set when it loads.
const _TEXTURES = (sun => "SunTexture.jpg", mercury => "MercuryTexture.jpg",
                   venus => "VenusTexture.jpg", earth => "EarthTexture.jpg",
                   moon => "MoonTexture.jpg", mars => "MarsTexture.jpg",
                   jupiter => "JupiterTexture.jpg", saturn => "SaturnTexture.jpg",
                   uranus => "UranusTexture.jpg", neptune => "NeptuneTexture.jpg")

# ---------------------------------------------------------------------------
# SPICE access.
#
# CSPICE is not thread-safe: two threads inside it abort the Julia process. Every call into it
# holds this lock, so an ephemeris call from threads waits its turn.
# ---------------------------------------------------------------------------

const _SPICE_LOCK = ReentrantLock()

# The default kernels are in the pool, or the user has taken over its configuration.
const _KERNELS_READY = Threads.Atomic{Bool}(false)

"""
    spice_lock() -> ReentrantLock

The lock held around every call into CSPICE, which is not thread-safe. A package that calls
SPICE directly holds it too.
"""
spice_lock() = _SPICE_LOCK

"""
    require_spice_kernels()

Make sure the default kernels are in the SPICE pool, downloading any that are not cached. A
package that reads the pool through SPICE directly calls this first.

Throws when a kernel cannot be fetched, naming each file, where it comes from and where to
put it by hand.
"""
@inline require_spice_kernels() = (_KERNELS_READY[] || _fetch_default_kernels(); nothing)

@noinline function _fetch_default_kernels()
    @lock _SPICE_LOCK begin
        _KERNELS_READY[] && return nothing
        cache_dir = get_spice_directory()
        failures = String[]
        for filename in DEFAULT_KERNELS
            source = kernel_source(filename)
            try
                ensure_kernel_download(cache_dir, source.filename, source.url;
                                       sha256 = source.sha256, size = source.size)
            catch err
                err isa InterruptException && rethrow()
                push!(failures, "  $(filename)\n    from $(source.url)\n    $(_one_line(err))")
            end
        end
        isempty(failures) || error(
            "Ephemeris and lunar orientation need SPICE kernels that are not in the cache and " *
            "could not be downloaded:\n$(join(failures, "\n"))\nTo install them by hand, " *
            "download each file from the address shown and put it in\n  $(cache_dir)\nthen call " *
            "again; there is no need to restart Julia.")
        _KERNELS_READY[] = _load_cached_kernels(cache_dir)
    end
    return nothing
end

# Load each default kernel whose file is in the cache at the right size. Reports whether all of
# them are loaded. Called with the lock held, or from `__init__`.
function _load_cached_kernels(cache_dir)
    all_loaded = true
    for filename in DEFAULT_KERNELS
        path = joinpath(cache_dir, filename)
        if !isfile(path) || filesize(path) != kernel_source(filename).size
            all_loaded = false
        elseif !is_kernel_loaded(filename)
            furnsh(path)
        end
    end
    return all_loaded
end

# SPICE takes seconds of TDB past J2000. The date is carried in two parts through the
# subtraction, so the result keeps the precision of the `Time`.
function _ephemeris_seconds(t::Time{<:Union{AbstractFloat, Integer, Rational}})
    tdb = t.tdb
    et = ((tdb.jd1 - 2451545.0) + tdb.jd2) * 86400.0
    isfinite(et) || throw(ArgumentError("The ephemeris cannot be read at $(t): it is not finite."))
    return Float64(et)
end

_ephemeris_seconds(t::Time{T}) where {T} = throw(ArgumentError(
    "The ephemeris is read from SPICE, which is not differentiable, so it cannot take a time of " *
    "element type $(T). Differentiate with respect to the state instead, or use an orientation " *
    "model written in Julia."))

_bare_date_error(name) = ArgumentError(
    "$(name) takes an `AstroEpochs.Time`, which knows its time scale; a bare Julian date does " *
    "not. For a TDB Julian date, use `Time(jd, TDB(), JD())`.")

"""
    translate(from::CelestialBody, to::CelestialBody, t::Time) -> Vector{Float64}

Position of `to` relative to `from` at the epoch `t`, in ICRF.

# Arguments
- `from`: Observing/origin body.
- `to`: Target body.
- `t`: Epoch, in any time scale. The ephemeris is read in TDB, using both parts of the date.

# Returns
The 3-element position `[x, y, z]`, in km.

# Notes
SPICE labels this frame `"J2000"`; it is ICRF-aligned and is not the FK5 mean equator and
equinox of J2000. Safe to call from threads, one at a time, since SPICE is not thread-safe.
Converting `t` from UTC or TT costs several times the lookup, so a loop over many bodies at one
epoch is faster with `t.tdb`.

# Examples
```jldoctest
julia> using AstroUniverse, AstroEpochs

julia> t = Time(2458018.0, TDB(), JD());

julia> round.(translate(earth, moon, t); digits = 3)
3-element Vector{Float64}:
 -375694.599
  -96115.682
  -12226.883
```
"""
function translate(from::CelestialBody, to::CelestialBody, t::Time)
    et = _ephemeris_seconds(t)
    require_spice_kernels()
    pos, _lt = @lock _SPICE_LOCK spkpos(string(to.naifid), et, "J2000", "NONE", string(from.naifid))
    return pos
end

translate(::CelestialBody, ::CelestialBody, ::Real) = throw(_bare_date_error("translate"))

"""
    translate_state(from::CelestialBody, to::CelestialBody, t::Time) -> Vector{Float64}

Position **and velocity** of `to` relative to `from` at the epoch `t`, in ICRF.

# Arguments
- `from`: Observing/origin body.
- `to`: Target body.
- `t`: Epoch, in any time scale. The ephemeris is read in TDB, using both parts of the date.

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
- **Not differentiable.** SPICE is a C library, so an origin shift built on
  this sits outside any AD path.

# Examples
```julia
using AstroUniverse, AstroEpochs
s_em = translate_state(earth, moon, Time(2458018.0, TDB(), JD()))
```
"""
function translate_state(from::CelestialBody, to::CelestialBody, t::Time)
    et = _ephemeris_seconds(t)
    require_spice_kernels()
    state, _lt = @lock _SPICE_LOCK spkezr(string(to.naifid), et, "J2000", "NONE", string(from.naifid))
    return state
end

translate_state(::CelestialBody, ::CelestialBody, ::Real) = throw(_bare_date_error("translate_state"))

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
    body.mu = newmu
    return body
end

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
    for (body, file) in _TEXTURES
        body.texture_file = joinpath(pkgdir(@__MODULE__), "data", file)
    end

    # Loading the package never reaches the network, and a problem with the cache never stops it
    # from loading: the constants and frame theories need no kernel, and the first call that
    # does need one fetches what is missing and reports what went wrong.
    try
        _KERNELS_READY[] = _load_cached_kernels(@get_scratch!("spice_kernels"))
    catch err
        @debug "SPICE kernels were not loaded at startup" exception = err
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
    body.mu = v
    return nothing
end

# ---------------------------------------------------------------------------
# Body orientation: IAU 2015 planetary orientation (pole + prime meridian +
# rates) and the orientation models. Before eop.jl, because its frame theories
# are Earth's orientation models.
# ---------------------------------------------------------------------------

include("iau_2015_orientation.jl")
include("orientation_models.jl")

# ---------------------------------------------------------------------------
# Earth Orientation Parameters (global service)
# ---------------------------------------------------------------------------

include("eop.jl")

end
