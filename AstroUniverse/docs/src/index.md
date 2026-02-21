```@meta
CurrentModule = AstroUniverse
```

# AstroUniverse

The AstroUniverse module provides models for celestial bodies, their physical properties, and related utilities for astrodynamics applications. It includes predefined celestial body objects with standard gravitational parameters and other physical constants commonly used in orbital mechanics.

The module automatically downloads and manages SPICE kernels (NASA's ephemeris data) to provide accurate celestial body positions and orientations using Scratch.jl. 

## Quick Start

The example below shows how to access predefined celestial bodies and their properties and how to add a celestial body:

```julia
using AstroUniverse

# Access predefined celestial bodies
earth.mu
venus.naifid

# Create a custom body
phobos = CelestialBody(
    name = "Phobos",
    naifid = 401,                    # NAIF ID for Phobos
    mu = 7.0875e-4,       # km³/s² (gravitational parameter)
    equatorial_radius = 11.1,                   # km (mean radius)
)

# Define a texture map for 3D graphics
custom_body = CelestialBody(
    texture_file = "/path/to/texture.jpg"
)
```

Note: built-in bodies include `sun`, `mercury`, `venus`, `earth`, `moon`, `mars`, `jupiter`, `saturn`, `uranus`, `neptune`, and `pluto`.

## SPICE Kernels and Ephemeris Data

AstroUniverse uses NASA's SPICE system for high-fidelity ephemeris calculations. SPICE provides accurate positions and orientations of celestial bodies across time.

### Default Kernels

The following kernels are automatically downloaded and loaded when you first use AstroUniverse:

- **naif0012.tls**: Leap second kernel for time conversions
- **de440.bsp**: Planetary ephemeris covering years 1550-2650

These kernels are stored using Scratch.jl and persist across Julia sessions, so they only download once.

### Basic Workflow

Loading additional SPICE kernels is a two-step process:

```julia
using AstroUniverse

# Step 1: Download kernel to persistent storage (only happens once)
download_spice_kernel("de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp")

# Step 2: Load kernel into SPICE system
load_spice_kernel("de440s.bsp")
```

This separation allows you to download kernels once and selectively load different combinations in different sessions.

### Common Use Cases

!!! warning "Large File Downloads"
    The examples below download ephemeris files ranging from ~13 MB to ~114 MB. These files are stored locally and only download once, but be aware of the initial download time and bandwidth usage.

**Extended Planetary Ephemeris:**
```julia
# Smaller file size, same time range as de440 (1550-2650)
download_spice_kernel("de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp")
load_spice_kernel("de440s.bsp")

# Extended time range planetary ephemeris (1550-2650, larger file)
download_spice_kernel("de430.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp")
load_spice_kernel("de430.bsp")
```

**Satellite Ephemerides:**
```julia
# Mars satellites (Phobos, Deimos)
download_spice_kernel("mar099.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar099.bsp")
load_spice_kernel("mar099.bsp")

# Jupiter satellites (Io, Europa, Ganymede, Callisto)
download_spice_kernel("jup365.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp")
load_spice_kernel("jup365.bsp")

# Saturn satellites
download_spice_kernel("sat441.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/sat441.bsp")
load_spice_kernel("sat441.bsp")
```

**Browse all available kernels:**
- Planetary: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
- Satellites: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/
- Leap seconds: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/

### Managing Kernels

**List downloaded kernels (files on disk):**
```julia
list_downloaded_spice_kernels()
# Output:
# Downloaded SPICE Kernels:
#   naif0012.tls              (5.3 KB)
#   de440.bsp                 (114.0 MB)
#   de440s.bsp                (31.0 MB)
#   jup365.bsp                (13.2 MB)
```

**List loaded kernels (in SPICE memory):**
```julia
list_cached_spice_kernels()
# Output:
# Cached SPICE Kernels (3 loaded):
#   naif0012.tls
#   de440.bsp
#   jup365.bsp
```

**Get storage directory:**
```julia
cache_dir = get_spice_directory()
println("Kernels stored at: ", cache_dir)
```

**Unload specific kernel:**
```julia
# Swap ephemeris versions
unload_spice_kernel("de440.bsp")
load_spice_kernel("de440s.bsp")
```

**Clear all loaded kernels:**
```julia
# Start fresh with custom configuration
unload_all_spice_kernels()
load_spice_kernel("naif0012.tls")
load_spice_kernel("de440s.bsp")
load_spice_kernel("jup365.bsp")
```

### Custom Configurations

For specialized analyses, you can create custom kernel configurations:

```julia
using AstroUniverse

# Clear default kernels
unload_all_spice_kernels()

# Load only what you need
load_spice_kernel("naif0012.tls")      # Required for time conversions
load_spice_kernel("de440s.bsp")        # Smaller planetary ephemeris
load_spice_kernel("mar099.bsp")        # Mars satellites only

# Your analysis code here...
```

### Advanced Usage

**Manual file placement:**
```julia
# Get storage directory
storage_dir = get_spice_directory()

# Copy a local kernel file to storage
cp("my_custom_kernel.bsp", joinpath(storage_dir, "my_custom_kernel.bsp"))

# Load it
load_spice_kernel("my_custom_kernel.bsp")
```

**Delete downloaded kernels:**

Kernels persist across sessions. To remove them, use standard filesystem operations:
```julia
storage_dir = get_spice_directory()
rm(joinpath(storage_dir, "old_kernel.bsp"))
```

The storage directory is managed by Scratch.jl and will be automatically cleaned if the package is removed.

## Texture Maps

Relatively small texture files for the Sun and planets are distributed with AstroUniverse in the AstroUniverse/data folder. Thanks to https://www.solarsystemscope.com/ for texture maps which are licensed using the Creative Commons 4.0 BY license.

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroUniverse]
Order = [:type, :function, :macro, :constant]
```
