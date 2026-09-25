# AstroManeuvers

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroManeuvers Overview

The AstroManeuvers module provides utilities and functions for orbital maneuver calculations in astrodynamics applications. The initial release includes impulsive maneuver models and functions for applying maneuvers to spacecraft objects.

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroManeuvers is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroManeuvers")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.2.0, the last version General carries, and
reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroManeuvers Documentation](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroManeuvers is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
