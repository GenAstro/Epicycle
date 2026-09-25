# AstroProp

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroProp Overview

AstroProp provides force models, orbital propagators, and stopping conditions for modelling spacecraft motion. AstroProp provides interfaces to the extensive numerical integration libraries in Julia's OrdinaryDiffEq.jl. The package can integate multiple spacecraft as a coupled dynamic system.

The module supports various force models including point mass gravity, and flexible stopping conditions based on orbital elements, time, or custom calculations. AstroProp is tested against the General Mission Analysis Tool (GMAT).

## What's New

- **Zonal Earth gravity (J₂–J₅)** — `HarmonicGravity(earth; model = Zonal())`. 
- **Exponential Earth atmosphere** — `Exponential()` density model for `AtmosphericDrag`.
- **Spherical atmospheric drag** — `AtmosphericDrag` force, reading `SphericalDrag` geometry
  from the spacecraft. Uses the drag-geometry type introduced in AstroModels; see the
  [AstroModels What's New](../AstroModels/README.md#whats-new).
- **Spherical solar radiation pressure** — `SolarRadiationPressure` force with a dual-cone
  eclipse shadow, reading `SphericalSRP` geometry from the spacecraft. Uses the SRP-geometry
  type introduced in AstroModels; see the
  [AstroModels What's New](../AstroModels/README.md#whats-new).
- **Faster and more flexible stopping conditions** — `StopAt` gains `detection` (`:discrete`
  default, `:continuous` escape hatch) and `rootfind_tol` kwargs. The default `:discrete`
  path polls the condition once per accepted step on `integrator.u` and, when a sign change
  appears, bisects on the Vern9 dense-output interpolant to locate the exact root — same
  root precision as before, ~2× faster on typical LEO drag workloads. Drop-in: existing
  `StopAt(subject, var, target; direction=…)` calls need no change.

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version AstroProp is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroProp")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroProp is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
