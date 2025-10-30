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


```

## CelestialBody Structure

The `CelestialBody` struct is the core data structure for representing astronomical objects in AstroUniverse. It stores the essential physical and orbital properties needed for astrodynamics calculations.

### Fields

- **`name::String`** - Human-readable name of the celestial body (e.g., "Earth", "Mars")
- **`naifid::Int`** - NAIF ID number used by NASA SPICE for ephemeris lookups
- **`mu::Float64`** - Standard gravitational parameter (GM) in km³/s²
- **`equatorial_radius::Float64`** - Equatorial radius in kilometers
- **`parent::Union{CelestialBody, Nothing}`** - Parent body for moons/satellites (Nothing for planets)

### Usage Examples

```julia
# Access predefined properties
println("Earth's GM: $(earth.mu) km³/s²")
println("Moon's parent: $(moon.parent.name)")

# Create asteroid Ceres
ceres = CelestialBody(
    name = "Ceres",
    naifid = 2000001,
    mu = 62.6284,           # km³/s²
    equatorial_radius = 469.7,  # km
    parent = nothing        # Dwarf planet, no parent
)
```

### Predefined Bodies

AstroUniverse includes all major solar system bodies:
- **Planets**: `sun`, `mercury`, `venus`, `earth`, `mars`, `jupiter`, `saturn`, `uranus`, `neptune`
- **Moons**: `moon` (Earth's), plus major moons of other planets
- **Dwarf Planets**: `pluto`

All predefined bodies include accurate physical parameters from NASA/JPL sources and are ready for immediate use in calculations.

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroUniverse]
Order = [:type, :function, :macro, :constant]
```
