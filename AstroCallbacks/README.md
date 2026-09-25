# AstroCallbacks

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroCallbacks Overview

AstroCallbacks provides a common interface for reading and writing spacecraft, maneuver, and celestial-body quantities. AstroProp uses these quantities for stopping conditions, AstroSolve uses them for variables and constraints, and reporting tools use them for columns and plots.

A quantity is a function of a subject, such as `semi_major_axis(sat)`, and frame-dependent quantities accept a coordinate system after the subject, as in `inclination(sat, EarthMJ2000Ec)`. The older Calc structs (`OrbitCalc`, `BodyCalc`, `ManeuverCalc`) are still accepted and will be deprecated in a future release. The documentation covers both.

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version AstroCallbacks is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroCallbacks")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroCallbacks Documentation](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroCallbacks is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
