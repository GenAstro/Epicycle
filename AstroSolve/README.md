# AstroSolve

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroSolve Overview

AstroSolve provides a framework for solving astrodynamics design problems through constrained optimization and event-driven architecture. AstroSolve builds an event sequence as a Directed Acyclic Graph (DAG), which will soon be fully differentiable using AD. The DAG model for optimization is inspired by the architecture used in NASA's Copernicus software.

AstroSolve's architecture is based on key user-facing components: Models (spacecraft, maneuvers, propagators), SolverVariables (optimization variables), Constraints (equality and inequality conditions), Events (discrete mission phases), and Sequence (DAG structure defining event dependencies).

## Installation

Versions through 0.4.0 are in Julia's General registry. From the next version AstroSolve is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroSolve")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.4.0, the last version General carries, and
reports nothing about the newer ones.

## Documentation

Full documentation is available at: [AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroSolve is source available under the [Gen Astro Source Available License 1.0](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
