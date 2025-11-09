# Overview

Epicycle is a comprehensive Julia package ecosystem for astrodynamics and space mission design. Built with a modular architecture, it provides a complete suite of tools for orbital mechanics, spacecraft trajectory analysis, and mission planning workflows.

The ecosystem consists of 11 specialized packages that work together to provide everything from fundamental astronomical calculations to advanced trajectory optimization:

- **Foundation Layer**: AstroBase (frames & constants), AstroStates (spacecraft states), AstroEpochs (time systems)
- **Dynamics Layer**: AstroDynamics (propagation), AstroForces (force modeling), AstroSimulation (mission simulation)  
- **Analysis Layer**: AstroTrajectories (orbital analysis), AstroManeuvers (spacecraft maneuvers), AstroEvents (mission events)
- **Mission Layer**: AstroMissions (mission design), AstroOptimization (trajectory optimization)
- **Integration Layer**: Epicycle (unified interface and workflows)

## Why New Software?

The astrodynamics field has historically relied on legacy Fortran codebases and proprietary tools that create barriers to innovation and collaboration. Existing open-source solutions often suffer from fragmentation, limited extensibility, or performance constraints.

### Modern Architecture Needs

Contemporary space missions require:
- **Scalable Performance**: Handle everything from CubeSat missions to interplanetary trajectories
- **Modular Design**: Compose mission-specific analysis workflows from reusable components  
- **Open Collaboration**: Enable researchers and engineers to build upon shared foundations
- **Rapid Prototyping**: Quickly test new algorithms and mission concepts
- **Integration Flexibility**: Work seamlessly with modern data science and optimization tools

### Julia Ecosystem Advantages

Julia provides unique benefits for astrodynamics applications:
- **Performance**: Near-C/Fortran speeds with high-level expressiveness
- **Composability**: Multiple dispatch enables seamless package integration
- **Scientific Computing**: Native differential equations, optimization, and automatic differentiation
- **Growing Ecosystem**: Active scientific computing community with modern tooling

Epicycle bridges the gap between specialized astrodynamics knowledge and modern computational capabilities.

## Design Philosophy

Epicycle is built around three core principles:

### Modularity & Composability
Each package focuses on a specific domain while maintaining clean interfaces. Users can combine components to build custom workflows without unnecessary complexity.

### Performance & Scalability  
From preliminary mission design to high-fidelity simulation, Epicycle scales efficiently across problem sizes while maintaining numerical accuracy.

### Extensibility & Interoperability
Open architecture allows researchers to extend capabilities and integrate with existing tools and datasets seamlessly.

## Software Architecture

### Package Hierarchy

| Package | Purpose | 
|:--------|:--------|
| `AstroBase` | Core abstractions and type hierarchy |
| `AstroStates` | Orbital state representations and conversions | 
| `AstroEpochs` | Time systems and epoch handling | 
| `AstroUniverse` | Celestial body models and ephemeris | 
| `AstroFrames` | Coordinate systems and transformations | 
| `AstroModels` | Spacecraft and physical models | 
| `AstroManeuvers`| Maneuver models and algorithms |
| `AstroCallbacks` | Quantities used in I/O, stopping conditions, cost, constraints |
| `AstroProp` | Trajectory propagation algorithms | 
| `AstroSolve` | Optimization and constraint solving | 
| `Epicycle` | The application.  `using Epicycle` loads everything. |

### Dependency Structure

The packages are organized in layers to ensure clean separation of concerns:

1. **Foundation**: Core types and constants that all other packages depend on
2. **Domain Specific**: Specialized functionality for states, time, coordinates, etc.
3. **Analysis Tools**: Higher-level capabilities built on foundation components
4. **Integration**: User-facing interfaces and complete workflows

## Key Features

### Comprehensive Force Modeling
- Gravitational perturbations (J2-J6, third-body, solid tides)
- Atmospheric drag with exponential and NRLMSISE-00 models
- Solar radiation pressure with cylindrical and spherical Earth shadow models
- Relativistic effects for high-precision applications

### Advanced Propagation Methods
- High-order Runge-Kutta integrators with adaptive stepping
- Specialized methods for different orbit regimes
- Event detection and handling during propagation
- Parallel processing for large trajectory sets

### Mission Design Tools
- Lambert problem solvers for transfer trajectory design
- Maneuver planning and optimization
- Launch window analysis
- Ground track and coverage analysis

### Modern Development Practices
- Comprehensive test coverage across all packages
- Continuous integration and automated testing
- Clear documentation with worked examples
- Type-stable implementations for optimal performance

## Licensing Model

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL v3.0
2) Evaluation and Education use Only
3) Commercial License

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

We gratefully acknowledge these projects and their maintainers, whose work makes Epicycle possible.

## Contributing

To protect both contributors and our company, we use the Linux Kernel's Developer's 
Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## Getting Started

New users should begin with the [Getting Started](getting_started.md) guide, which walks through installation and basic usage patterns. The [Tutorials](unit_examples.md) section provides hands-on examples for common astrodynamics tasks.

For developers interested in extending Epicycle, each package maintains its own documentation with implementation details and contribution guidelines.

## Software Architecture

### Package Hierarchy

| Package | Purpose | 
|:--------|:--------|
| `AstroBase` | Core abstractions and type hierarchy |
| `AstroStates` | Orbital state representations and conversions | 
| `AstroEpochs` | Time systems and epoch handling | 
| `AstroUniverse` | Celestial body models and ephemeris | 
| `AstroFrames` | Coordinate systems and transformations | 
| `AstroModels` | Spacecraft and physical models | 
| `AstroManeuvers`| Maneuver models and algorithms |
| `AstroCallbacks` | Quantities used in I/O, stopping conditions, cost, constraints |
| `AstroProp` | Trajectory propagation algorithms | 
| `AstroSolve` | Optimization and constraint solving | 
| `Epicycle` | The application.  `using Epicycle` loads everything. |

### Dependency Structure

The packages are organized in layers to ensure clean separation of concerns:

1. **Foundation**: Core types and constants that all other packages depend on
2. **Domain Specific**: Specialized functionality for states, time, coordinates, etc.
3. **Analysis Tools**: Higher-level capabilities built on foundation components
4. **Integration**: User-facing interfaces and complete workflows

## Key Features

### Comprehensive Force Modeling
- Gravitational perturbations (J2-J6, third-body, solid tides)
- Atmospheric drag with exponential and NRLMSISE-00 models
- Solar radiation pressure with cylindrical and spherical Earth shadow models
- Relativistic effects for high-precision applications

### Advanced Propagation Methods
- High-order Runge-Kutta integrators with adaptive stepping
- Specialized methods for different orbit regimes
- Event detection and handling during propagation
- Parallel processing for large trajectory sets

### Mission Design Tools
- Lambert problem solvers for transfer trajectory design
- Maneuver planning and optimization
- Launch window analysis
- Ground track and coverage analysis

### Modern Development Practices
- Comprehensive test coverage across all packages
- Continuous integration and automated testing
- Clear documentation with worked examples
- Type-stable implementations for optimal performance

## Licensing Model

We believe in the power of open source to foster innovation and community-driven 
development and also recognize the need for a sustainable business model and a model
that can handle export-controlled aerospace content. 

For these reasons, Epicycle is offered under a tri-licensing model. The license allows
users to choose between the following three options:

1) LGPL v3.0
2) Evaluation and Education use Only
3) Commercial License

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

We gratefully acknowledge these projects and their maintainers, whose work makes Epicycle possible.

## Contributing

To protect both contributors and our company, we use the Linux Kernel's Developer's 
Certificate of Origin (DCO) as detailed in CONTRIBUTING.txt.

## Getting Started

New users should begin with the [Getting Started](getting_started.md) guide, which walks through installation and basic usage patterns. The [Tutorials](unit_examples.md) section provides hands-on examples for common astrodynamics tasks.

For developers interested in extending Epicycle, each package maintains its own documentation with implementation details and contribution guidelines.