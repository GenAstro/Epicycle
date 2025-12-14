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
