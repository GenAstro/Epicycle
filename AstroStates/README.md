# AstroStates

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroStates Overview

The AstroStates module provides models, structs, utilities, and conversions for orbital state representations. A state representation is a set of quantities that uniquely define an orbit. Supported forms include Cartesian, Keplerian, Modified Equinoctial, and others.

The module offers multiple interfaces for transforming and storing states. Low‑level conversion functions (e.g., `cart_to_kep.jl`) can be used directly. A type system automatically provides concrete structs for each representation (e.g., `CartesianState`) and converts between all supported permutations. The `OrbitState` utility preserves type stability when the representation may change by storing the numeric state and a type tag in separate fields. The library supports automatic differentiation with ForwardDiff.jl and Zygote.jl.

AstroStates is tested against output from the General Mission Analysis Tool (GMAT) R2022a.

## What's New

- **Brouwer mean-element state types** — `BrouwerMeanShortState` and `BrouwerMeanLongState`.
- **`MeanSMA` calc** (in `AstroCallbacks`) — target or report the Brouwer long-period mean
  semi-major axis inside stopping conditions and constraints.

## Installation

```julia
using Pkg
Pkg.add("AstroStates")
```

## Documentation

Full documentation is available at: [AstroStates Documentation](https://genastro.github.io/Epicycle/AstroStates/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

AstroStates is open source under the [MIT License](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
