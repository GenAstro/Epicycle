# AstroSolve

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroSolve Overview

AstroSolve provides a framework for solving astrodynamics design problems through constrained optimization and event-driven architecture. AstroSolve builds an event sequence as a Directed Acyclic Graph (DAG), which will soon be fully differentiable using AD. The DAG model for optimization is inspired by the architecture used in NASA's Copernicus software.

AstroSolve's architecture is based on key user-facing components: Models (spacecraft, maneuvers, propagators), SolverVariables (optimization variables), Constraints (equality and inequality conditions), Events (discrete mission phases), and Sequence (DAG structure defining event dependencies).

## Installation

```julia
using Pkg
Pkg.add("AstroSolve")
```

## Documentation

Full documentation is available at: [AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/)

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
