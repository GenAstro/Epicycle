# AstroFrames

## Epicycle Overview

Epicycle is an application and package ecosystem for space mission design and navigation. It contains packages that handle astrodynamics models and algorithms that integrate seamlessly to allow users to setup and solve hard problems, fast.

## AstroFrames Overview

The AstroFrames package provides coordinate systems and reference frames for astrodynamics applications. AstroFrames provides types for defining coordinate systems with customizable origins and axes and conversions between different coordinates.

A coordinate system consists of an origin and axis system, allowing you to create frames like Earth-centered ICRF or spacecraft-relative VNB coordinates.

The initial release of AstroFrames only has types and supports ICRF, VNB, and Inertial axes. The long-term plan is to interface with Julia Space's SatelliteToolboxTransformations.jl for reference frame conversions. The first version of AstroFrames contains a minimal API to be consistent with Epicycle and SatelliteToolboxTransformations.jl. 


## Installation

```julia
using Pkg
Pkg.add("AstroFrames")
```

## Documentation

Full documentation is available at: [AstroFrames Documentation](https://genastro.github.io/Epicycle/AstroFrames/dev/)

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
