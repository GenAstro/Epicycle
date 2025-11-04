# Epicycle.jl 

Welcome to Epicycle, an application for space mission analysis, trajectory optimization, and navigation. 

## Design Philosophy

Epicycle embraces a modular architecture that adapts to your workflow. Whether you need a single capability like orbital state conversions (`using AstroStates`) or the complete mission design environment (`using Epicycle`), the system scales from focused analysis to comprehensive mission planning.

The framework prioritizes production readiness through extensive validation against GMAT test cases and comprehensive testing across Windows, macOS, and Linux platforms. Every algorithm undergoes rigorous verification to ensure reliability in mission-critical applications.

Performance remains central to the design, leveraging Julia's computational speed while maintaining compatibility with automatic differentiation frameworks like ForwardDiff.jl and Zygote.jl. This enables gradient-based optimization and sensitivity analysis throughout the mission design process.

## Core Capabilities

Epicycle delivers a complete, integrated application for astrodynamics analysis that prioritizes breadth and extensibility. The system provides working implementations across the full mission design workflow - from orbital state representations and coordinate transformations to trajectory propagation and optimization - with interfaces designed for expanding the model library as capabilities mature.

The application architecture handles essential orbital state representations including Cartesian, Keplerian, and Modified Equinoctial elements, enabling seamless transitions between formulations. Trajectory propagation integrates with Julia's differential equation ecosystem, while the optimization framework connects SNOW-based algorithms with IPOPT for nonlinear programming applications.

Rather than focusing on depth in individual models, Epicycle emphasizes the integration of components into a cohesive application. The extensible interface design supports systematic expansion of the model library, starting with enhanced physical models and growing toward more sophisticated capabilities as the framework matures.

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
| `AstroFun` | Quantities used in I/O, stopping conditions, cost, constraints.
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

### 📚 [Using Epicycle](welcome.md)
Learn how to install, configure, and use the Epicycle ecosystem for your mission analysis needs.

### 🚀 [Tutorials](unit_examples.md) 
Hands-on examples and complete mission workflows using real aerospace scenarios.

---
