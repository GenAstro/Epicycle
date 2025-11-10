# AstroProp

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroProp Overview

AstroProp provides force models, orbital propagators, and stopping conditions for modelling spacecraft motion. AstroProp provides interfaces to the extensive numerical integration libraries in Julia's OrdinaryDiffEq.jl. The package can integate multiple spacecraft as a coupled dynamic system.

The module supports various force models including point mass gravity, and flexible stopping conditions based on orbital elements, time, or custom calculations. AstroProp is tested against the General Mission Analysis Tool (GMAT).

## Installation

```julia
using Pkg
Pkg.add("AstroProp")
```

## Documentation

Full documentation is available at: [AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/)

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
