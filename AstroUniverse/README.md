# AstroUniverse

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroUniverse Overview

The AstroUniverse module provides models for celestial bodies, their physical properties, and related utilities for astrodynamics applications. It includes predefined celestial body objects with standard gravitational parameters and other physical constants commonly used in orbital mechanics.

AstroUniverse uses SPICE.jl for ephemeris modelling. The module automatically downloads and manages SPICE kernels (NASA's ephemeris data) to provide accurate celestial body positions and orientations using Scratch.jl.

## Installation

```julia
using Pkg
Pkg.add("AstroUniverse")
```

The first `using AstroUniverse` downloads about 110 MB of SPICE kernels, most of it the DE440
ephemeris, and needs a network connection. Later sessions use the stored copy and work offline.

## Documentation

Full documentation is available at: [AstroUniverse Documentation](https://genastro.github.io/Epicycle/AstroUniverse/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroUniverse is open source under the [MIT License](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
