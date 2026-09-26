# Welcome to Epicycle

[![CI](https://github.com/GenAstro/Epicycle/workflows/CI/badge.svg)](https://github.com/GenAstro/Epicycle/actions)
[![codecov](https://codecov.io/gh/GenAstro/Epicycle/branch/main/graph/badge.svg?token=FNHOVC5O5N)](https://codecov.io/gh/GenAstro/Epicycle)
[![Documentation (stable)](https://img.shields.io/badge/docs-stable-blue.svg)](https://genastro.github.io/Epicycle/Epicycle/stable/)
[![Documentation (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://genastro.github.io/Epicycle/Epicycle/dev/)

Epicycle is an application for space systems with a nod to the giants before us and to how far we've come. 

## Documentation

### The Epicycle Application

- **[Epicycle Application](https://genastro.github.io/Epicycle/Epicycle/dev/)** - Complete application guide

### Component Documentation

- **[EpicycleBase Documentation](https://genastro.github.io/Epicycle/EpicycleBase/dev/)** - Foundation types and abstractions
- **[AstroRoutines Documentation](https://genastro.github.io/Epicycle/AstroRoutines/dev/)** - Low-level astrodynamics routines, usable without the rest of Epicycle
- **[AstroStates Documentation](https://genastro.github.io/Epicycle/AstroStates/dev/)** - State vector representations and conversions
- **[AstroEpochs Documentation](https://genastro.github.io/Epicycle/AstroEpochs/dev/)** - Time systems and epoch handling
- **[AstroUniverse Documentation](https://genastro.github.io/Epicycle/AstroUniverse/dev/)** - Celestial body models and ephemeris
- **[AstroFrames Documentation](https://genastro.github.io/Epicycle/AstroFrames/dev/)** - Coordinate systems and transformations
- **[AstroModels Documentation](https://genastro.github.io/Epicycle/AstroModels/dev/)** - Spacecraft and physical models
- **[AstroManeuvers Documentation](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)** - Maneuver models and algorithms
- **[AstroCallbacks Documentation](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)** - I/O, stopping conditions, cost, and constraints
- **[AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/)** - Trajectory propagation algorithms
- **[AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/)** - Optimization and constraint solving
- **[EpicycleIO Documentation](https://genastro.github.io/Epicycle/EpicycleIO/dev/)** - Plotting, three-dimensional views, and reports, drawn in a browser

## What's New

The beta release adds targeting, optimal control and orbit determination to the propagation core,
and brings EpicycleIO's browser-drawn plots, three-dimensional views and reports into the public
tree. It also adds force models, a new callback interface, new state types, leap second and EOP
management, and 32 runnable example scripts rendered as documentation pages.

[Current Capabilities](https://genastro.github.io/Epicycle/Epicycle/dev/#Current-Capabilities) contains 
a complete summary of Epicycle's functionality.

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

Each package is licensed under the terms in its own `LICENSE.md`. [LICENSE.md](LICENSE.md) lists
which package is under which license.

## Notes
Claude and ChatGPT are used in the development of Epicycle.
