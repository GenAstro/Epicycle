```@meta
CurrentModule = AstroUniverse
```

# AstroUniverse

The AstroUniverse module provides models for celestial bodies, their physical properties, and related utilities for astrodynamics applications. It includes predefined celestial body objects with standard gravitational parameters and other physical constants commonly used in orbital mechanics.

The module automatically downloads and manages SPICE kernels (NASA's ephemeris data) to provide accurate celestial body positions and orientations. All major solar system bodies are included with their standard gravitational parameters, allowing for immediate use in trajectory analysis without manual setup.

AstroUniverse integrates seamlessly with other Epicycle modules, providing the foundational celestial body models used throughout the astrodynamics toolkit.

## Quick Start

The example below shows how to access predefined celestial bodies and their properties:

```julia
using AstroUniverse

# Access predefined celestial bodies
earth_body = earth
mars_body = mars

# Get gravitational parameters
μ_earth = get_gravparam(earth)  # km³/s²
μ_mars = get_gravparam(mars)

# All major solar system bodies are available
bodies = [sun, mercury, venus, earth, moon, mars, jupiter, saturn, uranus, neptune, pluto]
```

## Features

- **Predefined celestial bodies**: All major solar system bodies with standard parameters
- **SPICE kernel management**: Automatic download and loading of NASA ephemeris data  
- **Gravitational parameters**: Standard GM values for all bodies
- **Type-safe interfaces**: Consistent API across all celestial body objects
- **Scratch space management**: Kernels cached locally to avoid repeated downloads
