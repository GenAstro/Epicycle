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

## Table of Contents

```@index
```

## API Reference

```@autodocs
Modules = [AstroUniverse]
Order = [:type, :function, :macro, :constant]
```
