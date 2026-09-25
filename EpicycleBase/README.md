# EpicycleBase

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## EpicycleBase Overview

EpicycleBase provides core abstract types shared across the Epicycle ecosystem. These fundamental types form the public type hierarchy used by higher-level packages including AstroStates, AstroEpochs, AstroFrames, AstroProp, and AstroSolve.

The module defines base abstractions for variables (states, controls, time, parameters), functions, calculation variables, orbit state types, and geometric points. This shared foundation enables consistent interfaces and type hierarchies across all Epicycle packages.

## Installation

```julia
using Pkg
Pkg.add("EpicycleBase")
```

## Documentation

Full documentation is available at: [EpicycleBase Documentation](https://genastro.github.io/Epicycle/EpicycleBase/dev/)

## Contributing

The terms for contributions are in [LICENSE.md](LICENSE.md).

## License

EpicycleBase is open source under the [MIT License](LICENSE.md).

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
