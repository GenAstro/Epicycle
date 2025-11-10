# AstroUniverse

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroUniverse Overview

The AstroUniverse module provides models for celestial bodies, their physical properties, and related utilities for astrodynamics applications. It includes predefined celestial body objects with standard gravitational parameters and other physical constants commonly used in orbital mechanics.

AstroUniverse uses SPICE.jl for ephemeris modelling. The module automatically downloads and manages SPICE kernels (NASA's ephemeris data) to provide accurate celestial body positions and orientations using Scratch.jl.

## Installation

```julia
using Pkg
Pkg.add("AstroUniverse")
```

## Documentation

Full documentation is available at: [AstroUniverse Documentation](https://genastro.github.io/Epicycle/AstroUniverse/dev/)

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
