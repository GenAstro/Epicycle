# Welcome to Epicycle

[![CI](https://github.com/GenAstro/Epicycle/workflows/CI/badge.svg)](https://github.com/GenAstro/Epicycle/actions)
[![codecov](https://codecov.io/gh/GenAstro/Epicycle/branch/main/graph/badge.svg?token=FNHOVC5O5N)](https://codecov.io/gh/GenAstro/Epicycle)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://genastro.github.io/Epicycle/Epicycle/dev/)

Epicycle is a Julia package ecosystem for astrodynamics and space mission design, built with a modular architecture that spans mission analysis workflows from preliminary design through trajectory optimization. The current implementation - which is the initial release - focuses on establishing a robust architecture with key components—coordinate systems, time standards, spacecraft state representations, basic propagation, targeting and optimization—while providing extensible interfaces for systematic expansion toward operational mission design and navigation capabilities.

The ecosystem consists of eleven specialized packages organized in architectural layers, from core abstractions (AstroBase) through utilities (AstroEpochs, AstroStates, and AstroUniverse) to  integrated workflows (AstroManeuvers, AstroProp, and AstroSolve). The structure enables users to use low level utilities independently of the full system, or compose mission-specific analyses using interfaces designed to solve complex design problems, fast. The architecture in the first release is well tested and documented, and provides a clear pathway for extending capabilities toward high-fidelity force modeling, advanced optimization algorithms, and operational navigation applications as the framework matures.

## Documentation

### The Epicycle Application

- **[Epicycle Application](https://genastro.github.io/Epicycle/Epicycle/dev/)** - Complete application guide

### Component Documentation

- **[AstroBase Documentation](https://genastro.github.io/Epicycle/AstroBase/dev/)** - Foundation types and abstractions
- **[AstroStates Documentation](https://genastro.github.io/Epicycle/AstroStates/dev/)** - State vector representations and conversions
- **[AstroEpochs Documentation](https://genastro.github.io/Epicycle/AstroEpochs/dev/)** - Time systems and epoch handling
- **[AstroUniverse Documentation](https://genastro.github.io/Epicycle/AstroUniverse/dev/)** - Celestial body models and ephemeris
- **[AstroFrames Documentation](https://genastro.github.io/Epicycle/AstroFrames/dev/)** - Coordinate systems and transformations
- **[AstroModels Documentation](https://genastro.github.io/Epicycle/AstroModels/dev/)** - Spacecraft and physical models
- **[AstroManeuvers Documentation](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)** - Maneuver models and algorithms
- **[AstroCallbacks Documentation](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)** - I/O, stopping conditions, cost, and constraints
- **[AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/)** - Trajectory propagation algorithms
- **[AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/)** - Optimization and constraint solving

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