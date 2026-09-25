# Epicycle

## Overview

Epicycle is a Julia package ecosystem for astrodynamics and space mission design, built with a modular architecture that spans mission analysis workflows from preliminary design through trajectory optimization and orbit determination. 

## Package Architecture

The Epicycle ecosystem implements a layered architecture where each package provides focused functionality while maintaining clean interfaces for composition in workflows. Users can access the complete system through `using Epicycle` in Julia, or leverage individual components independently for specialized applications.

| Package | Purpose | 
|:--------|:--------|
| `Epicycle` | Integration layer providing unified interface and common workflows |
| `EpicycleBase` | Core abstractions, type hierarchy, and fundamental constants |
| `AstroRoutines`| Low level astrodynamics functions and algorithms | 
| `AstroStates` | Spacecraft state representations and state transformations | 
| `AstroEpochs` | Time systems, epoch handling, and temporal conversions | 
| `AstroUniverse` | Celestial body models, ephemeris data, and gravitational parameters | 
| `AstroFrames` | Reference frames and coordinate system transformations | 
| `AstroModels` | Spacecraft and physical models  | 
| `AstroManeuvers` | Maneuver models and algorithms |
| `AstroCallbacks` | Utilities for constraints, objectives, and events |
| `AstroProp` | Numerical integration and trajectory propagation methods | 
| `AstroSolve` | Optimization algorithms and constraint solving capabilities |
| `EpicycleIO` | Plotting, reporting, and 3D trajectory visualization |
| `EpicycleGraphics` | 3D graphics and data visualization |

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version Epicycle is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("Epicycle")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

## Documentation

Full documentation is available at: [Epicycle Documentation](https://genastro.github.io/Epicycle/Epicycle/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

Packages in the Epicycle ecosystem are not all under the same license. Each package carries its own
license file, which is the authoritative statement of its terms, so read the one in the package you
intend to use.

The `Epicycle` package itself is under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude and ChatGPT are used in the development of Epicycle.
