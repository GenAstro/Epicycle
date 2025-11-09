# Epicycle: An Application for Space Mission Design and Navigation

Epicycle is a Julia package ecosystem for astrodynamics and space mission design, built with a modular architecture that spans mission analysis workflows from preliminary design through trajectory optimization. The current implementation - which is the initial release - focuses on establishing robust architecture with key components—coordinate systems, time standards, spacecraft state representations, basic propagation, targeting and optimization—while providing extensible interfaces for systematic expansion toward operational mission design and navigation capabilities.

The ecosystem consists of eleven specialized packages organized in architectural layers, from core abstractions through domain-specific functionality to integrated workflows. This structure enables users to compose mission-specific analyses from reusable components while maintaining a clear pathway for extending capabilities toward high-fidelity force modeling, advanced optimization algorithms, and operational navigation applications as the framework matures.

## Architecture and Components

### Development Environment

Epicycle is designed to work seamlessly with modern development tools:

- **Visual Studio Code** with Julia Language Server provides comprehensive IDE features
- **Code Intelligence** through language-based completion and GitHub Copilot integration
- **Interactive Development** with integrated debugger and REPL workflow
- **Workspace Management** for organizing multi-package projects

### Package Architecture

The Epicycle ecosystem implements a layered architecture where each package provides focused functionality while maintaining clean interfaces for composition in workflows. Users can access the complete system through `using Epicycle` in Julia, or leverage individual components independently for specialized applications. The packages are organized from the integration layer down through foundational components: 

| Package | Purpose | 
|:--------|:--------|
| `Epicycle` | Integration layer providing unified interface and common workflows |
| `AstroBase` | Core abstractions, type hierarchy, and fundamental constants |
| `AstroStates` | Spacecraft state representations and state transformations | 
| `AstroEpochs` | Time systems, epoch handling, and temporal conversions | 
| `AstroUniverse` | Celestial body models, ephemeris data, and gravitational parameters | 
| `AstroFrames` | Reference frames and coordinate system transformations | 
| `AstroModels` | Spacecraft and physical models  | 
| `AstroManeuvers` | Maneuver models and algorithms |
| `AstroCallbacks` | Utilities for constraints, objectives, and events |
| `AstroProp` | Numerical integration and trajectory propagation methods | 
| `AstroSolve` | Optimization algorithms and constraint solving capabilities | 

## Current Status

The initial release of Epicycle provides a core architecture that integrates fundamental algorithms and models into an intuitive application for solving engineering problems. The model packages establish essential functionality while the architecture is designed for systematic expansion. The implementation maintains loose coupling between packages so lower-level components can be used independently of the full Epicycle application.

The system includes comprehensive test suites and documentation across all packages, with 95% code coverage and continuous integration on GitHub.  The software has been tested and validated on macOS, Linux, and Windows environments, with the architectural foundation established and core functionality implemented, making Epicycle ready for mission analysis workflows while maintaining an extensible design for advanced capabilities. 

Astrodynamics computations are tested against the General Mission Analysis Tool (GMAT) R2022a. Time conversion calculations are tested against Astropy.Time.

## Why New Software

Julia is a modern, high-performance language designed for technical computing. It combines the ease of use found in MATLAB and Python with the performance of C/C++.

Most aerospace tools require custom scripting interfaces or domain-specific languages. Julia serves as both the implementation language and the user interface, providing direct access to the full computational ecosystem. The language's design emphasizes scientific computing and automatic differentiation, both essential for aerospace optimization and navigation applications.

- **High-Performance Numerical Analysis** - Julia is designed for high-performance numerical analysis, making it suitable for complex scientific computations.
- **Efficient Linear Algebra** - Julia excels in linear algebra with efficient matrix operations and optimized algorithms.
- **Differential Equations** - Julia provides advanced features for solving differential equations, making it suitable for complex scientific and engineering problems.
- **SciML Machine Learning** - Julia seamlessly integrates with SciML for machine learning, enhancing the capabilities for scientific machine learning applications.
- **Optimization Tools** - Julia interfaces seamlessly with optimization tools like SNOPT and IPOPT, facilitating the handling of complex tasks in technical computing.

## Acknowledgments

Epicycle builds upon the foundational work of many contributors to the aerospace and scientific computing communities:

**Astrodynamics Standards**
- NASA GMAT Development Team for orbital mechanics specifications and validation test cases
- David Vallado, "Fundamentals of Astrodynamics and Applications, 4th Edition" (2013), Microcosm Press, for mathematical formulations and algorithmic references
- The Astropy Project for rigorous time system standards and implementations

**Julia Scientific Computing Ecosystem**
- SciML Organization for OrdinaryDiffEq.jl used in AstroProp
- Julia Astro community for SPICE.jl used in AstroUniverse
- BYU FLOW Lab for SNOW.jl used in AstroSolve
- Julia Space Mission Design for the TEMPO.jl library used in AstroEpochs 
- Wächter & Biegler for the IPOPT nonlinear programming solver

**Open Source Foundations**
- Julia Computing and contributors to the Julia language
- The Documenter.jl team for documentation generation
- GitHub Actions and the CI/CD community for automated testing infrastructure
- Visual Studio Code, used to develop Epicycle and the recommended user interface

We gratefully acknowledge these projects and their maintainers, whose work makes Epicycle possible.

### Core Contributors

- Steve Hughes (steven.hughes at genastro.org), architect and lead developer.

## License

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL V3.0
2) Evaluation and Education use Only
3) Commercial License

See LICENSE.txt for the terms of each license option. For licensing questions contact licensing [at] genastro.com

## Contributing

To protect both contributors and our company, we use the Linux Kernel's Developer's 
Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## Support

For support, including technical support and services to apply Epicycle to your application, contact support [at] genastro.com

## What is an Epicycle?

Humankind has been studying planetary motion for millennia. An epicycle is a geometric theory developed by Ptolemy to explain why planets appear to reverse direction and perform small loops in their celestial paths. While this model represented a significant advancement over earlier theories, it was ultimately incorrect—and it would be nearly 1500 years before Kepler developed a more accurate framework for understanding orbital mechanics.

We've come remarkably far in our understanding, yet fundamental questions remain. Either our theories of relativity, quantum mechanics, or both may be incomplete—reminding us that scientific discovery is an ongoing journey.

The Epicycle software is a tribute to the brilliant minds who came before us, celebrating how far we've advanced while embracing the excitement of continuing to push the boundaries of knowledge and make new discoveries. 