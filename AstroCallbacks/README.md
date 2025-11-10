# AstroCallbacks

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroCallbacks Overview

AstroCallbacks provides fundamental calculations used throughout the Epicycle ecosystem. These calculations can be used in targeting and optimization variables and constraints, propagator stopping conditions, I/O operations, and eventually solver-for and consider parameters for orbit determination.

The module serves as a foundation for higher-level astrodynamics operations, providing core mathematical functions and utilities needed across multiple Epicycle packages. AstroCallbacks contains the "Calcs" - a collection of structs that provide a unified interface for setting and getting quantities throughout Epicycle.

The Calc framework provides type-stable access to commonly needed astrodynamics quantities while maintaining automatic differentiation compatibility. Calc types include OrbitCalc (orbital elements), BodyCalc (celestial body properties), and ManeuverCalc (Δv and thrust properties).

## Installation

```julia
using Pkg
Pkg.add("AstroCallbacks")
```

## Documentation

Full documentation is available at: [AstroCallbacks Documentation](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)

## Contributing to Epicycle 

Contributing is easy.

1. Fork the project
2. Create a new feature branch
3. Make your changes
4. Submit a pull request

We use the Linux Kernel's Developer's Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## License

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL V3.0
2) Evaluation and Education use Only
3) Commercial License

See LICENSE.txt for terms each license option.  For commercial licensing, 
email licensing at genastro.org.

## Notes
Claude Sonnet and ChatGPT are used in the development of Epicycle.
