# AstroModels

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroModels Overview

The AstroModels module provides physical models for astrodynamics applications, including spacecraft representations with state, time, and mass properties. The primary model is the `Spacecraft` type, which integrates orbital states from AstroStates with temporal information from AstroEpochs.

The module supports multiple initialization patterns for spacecraft objects, allowing users to specify orbital states directly or construct them from component data. All models are designed for compatibility with automatic differentiation libraries.

## What's New

- **`SphericalDrag` spacecraft geometry** — isotropic (spherical-body) drag geometry
  (`SphericalDrag(; c_d, drag_area)`), attached to `Spacecraft.drag`. Read by
  `AstroProp.AtmosphericDrag` at each integration step.
- **`SphericalSRP` spacecraft geometry** — isotropic (spherical-body) SRP geometry
  (`SphericalSRP(; c_r, srp_area)`), attached to `Spacecraft.srp`. Read by
  `AstroProp.SolarRadiationPressure` at each integration step.
- **`Spacecraft.save_history`** — new field to opt out of trajectory-segment storage during
  propagation. Default `true` (unchanged behavior). When `false`, `sc.state` and `sc.time`
  still update at the end of `propagate!`; only the segment push into `sc.history` is
  skipped. The returned `ODESolution` still carries the full trajectory. Useful for
  benchmarking or when only the final state is needed.

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroModels is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add both registries once, then install as usual:

```julia
using Pkg
pkg"registry add General https://github.com/GenAstro/GenAstroRegistry.git"
pkg"add AstroModels"
```

General is named in that command for two reasons: the dependencies live there, and Julia
installs it by itself only while no registry is present at all, so adding the Gen Astro
registry alone would leave it out. Installing without the Gen Astro registry resolves to 0.2.0,
the last version General carries, and reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroModels Documentation](https://genastro.github.io/Epicycle/AstroModels/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroModels is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
