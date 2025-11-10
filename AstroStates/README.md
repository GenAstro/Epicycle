# AstroStates

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroStates Overview

The AstroStates module provides models, structs, utilities, and conversions for orbital state representations. A state representation is a set of quantities that uniquely define an orbit. Supported forms include Cartesian, Keplerian, Modified Equinoctial, and others.

The module offers multiple interfaces for transforming and storing states. Low‑level conversion functions (e.g., `cart_to_kep.jl`) can be used directly. A type system automatically provides concrete structs for each representation (e.g., `CartesianState`) and converts between all supported permutations. The `OrbitState` utility preserves type stability when the representation may change by storing the numeric state and a type tag in separate fields. The library supports automatic differentiation with ForwardDiff.jl and Zygote.jl.

AstroStates is tested against output from the General Mission Analysis Tool (GMAT) R2022a.

## Installation

```julia
using Pkg
Pkg.add("AstroStates")
```

## Documentation

Full documentation is available at: [AstroStates Documentation](https://genastro.github.io/Epicycle/AstroStates/dev/)

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
