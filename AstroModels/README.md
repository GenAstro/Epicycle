# AstroModels

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroModels Overview

The AstroModels module provides physical models for astrodynamics applications, including spacecraft representations with state, time, and mass properties. The primary model is the `Spacecraft` type, which integrates orbital states from AstroStates with temporal information from AstroEpochs.

The module supports multiple initialization patterns for spacecraft objects, allowing users to specify orbital states directly or construct them from component data. All models are designed for compatibility with automatic differentiation libraries.

## Installation

```julia
using Pkg
Pkg.add("AstroModels")
```

## Documentation

Full documentation is available at: [AstroModels Documentation](https://genastro.github.io/Epicycle/AstroModels/dev/)

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
