# Epicycle

## Overview

Epicycle is a Julia package ecosystem for astrodynamics and space mission design, built with a modular architecture that spans mission analysis workflows from preliminary design through trajectory optimization. The current implementation - which is the initial release - focuses on establishing a robust architecture with key components—coordinate systems, time standards, spacecraft state representations, basic propagation, targeting and optimization—while providing extensible interfaces for systematic expansion toward operational mission design and navigation capabilities.

The ecosystem consists of eleven specialized packages organized in architectural layers, from core abstractions (EpicycleBase) through utilities (AstroEpochs, AstroStates, and AstroUniverse) to  integrated workflows (AstroManeuvers, AstroProp, and AstroSolve). The structure enables users to use low level utilities independently of the full system, or compose mission-specific analyses using interfaces designed to solve complex design problems, fast. The architecture in the first release is well tested and documented, and provides a clear pathway for extending capabilities toward high-fidelity force modeling, advanced optimization algorithms, and operational navigation applications as the framework matures.

## Package Architecture

The Epicycle ecosystem implements a layered architecture where each package provides focused functionality while maintaining clean interfaces for composition in workflows. Users can access the complete system through `using Epicycle` in Julia, or leverage individual components independently for specialized applications.

| Package | Purpose | 
|:--------|:--------|
| `Epicycle` | Integration layer providing unified interface and common workflows |
| `EpicycleBase` | Core abstractions, type hierarchy, and fundamental constants |
| `AstroStates` | Spacecraft state representations and state transformations | 
| `AstroEpochs` | Time systems, epoch handling, and temporal conversions | 
| `AstroUniverse` | Celestial body models, ephemeris data, and gravitational parameters | 
| `AstroFrames` | Reference frames and coordinate system transformations | 
| `AstroModels` | Spacecraft and physical models  | 
| `AstroManeuvers` | Maneuver models and algorithms |
| `AstroCallbacks` | Utilities for constraints, objectives, and events |
| `AstroProp` | Numerical integration and trajectory propagation methods | 
| `AstroSolve` | Optimization algorithms and constraint solving capabilities |

## Installation

```julia
using Pkg
Pkg.add("Epicycle")
```

## Documentation

Full documentation is available at: [Epicycle Documentation](https://genastro.github.io/Epicycle/Epicycle/dev/)

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
