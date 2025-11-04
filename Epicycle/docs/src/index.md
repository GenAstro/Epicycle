# Epicycle.jl 

Welcome to Epicycle, an application for space mission analysis, trajectory optimization, and navigation. 

## Design Philosophy

Epicycle uses a modular architecture. You can import individual packages like `AstroStates` for specific functionality or `using Epicycle` to load the complete system.

The framework is tested against GMAT validation cases and runs on Windows, macOS, and Linux. It uses Julia for performance and supports automatic differentiation with ForwardDiff.jl and Zygote.jl.

## Core Capabilities

Epicycle provides an integrated application for astrodynamics analysis with focus on breadth and extensible interfaces. The system covers orbital state representations, coordinate transformations, trajectory propagation, and optimization.

The system handles Cartesian, Keplerian, and Modified Equinoctial orbital elements. Trajectory propagation uses Julia's differential equation solvers. Optimization connects SNOW algorithms with IPOPT.

Current implementation emphasizes application integration over individual model depth. Interfaces are designed for systematic expansion of the model library.


Julia, VS Code, Epicycle

### Package Architecture

| Package | Purpose | 
|:--------|:--------|
| `AstroBase` | Core abstractions and type hierarchy |
| `AstroStates` | Orbital state representations and conversions | 
| `AstroEpochs` | Time systems and epoch handling | 
| `AstroUniverse` | Celestial body models and ephemeris | 
| `AstroCoords` | Coordinate systems and transformations | 
| `AstroModels` | Spacecraft and physicsl models | 
| `AstroMan`| Maneuver models and algorithms |
| `AstroFun` | Quantities used in I/O, stopping conditions, cost, constraints
| `AstroProp` | Trajectory propagation algorithms | 
| `AstroSolve` | Optimization and constraint solving | 
| `Epicycle` | The application.  `using Epicycle` loads everything. 

## Acknowledgments

Epicycle builds upon the foundational work of many contributors to the aerospace and scientific computing communities:

**Astrodynamics Standards**
- NASA GMAT Development Team for orbital mechanics algorithms and validation test cases
- David Vallado for "Fundamentals of Astrodynamics and Applications" formulations
- The Astropy Project for rigorous time system standards and implementations

**Julia Scientific Computing Ecosystem**
- SciML Organization for OrdinaryDiffEq.jl and the broader differential equations ecosystem
- Julia Astro community for SPICE.jl and astronomical coordinate systems
- BYU FLOW Lab for SNOW.jl optimization framework
- Wächter & Biegler for the IPOPT nonlinear programming solver

**Open Source Foundations**
- Julia Computing and contributors for the Julia language
- The Documenter.jl team for documentation generation
- GitHub Actions and the CI/CD community for automated testing infrastructure
TODO: Visual Studio Code

We gratefully acknowledge these projects and their maintainers, whose work makes Epicycle possible.

### Support


### License

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL V3.0
2) Evaluation and Education use Only
3) Commercial License

### Contributing

To protect both contributors and our company, we use the Linux Kernel's Developer's 
Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

### 📚 [Using Epicycle](welcome.md)
Learn how to install, configure, and use the Epicycle ecosystem for your mission analysis needs.

### 🚀 [Tutorials](unit_examples.md) 
Hands-on examples and complete mission workflows using real aerospace scenarios.

---
