# AstroManeuvers

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroManeuvers Overview

The AstroManeuvers module provides utilities and functions for orbital maneuver calculations in astrodynamics applications. The initial release includes impulsive maneuver models and functions for applying maneuvers to spacecraft objects.

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroManeuvers is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add both registries once, then install as usual:

```julia
using Pkg
pkg"registry add General https://github.com/GenAstro/GenAstroRegistry.git"
pkg"add AstroManeuvers"
```

General is named in that command for two reasons: the dependencies live there, and Julia
installs it by itself only while no registry is present at all, so adding the Gen Astro
registry alone would leave it out. Installing without the Gen Astro registry resolves to 0.2.0,
the last version General carries, and reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroManeuvers Documentation](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroManeuvers is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
