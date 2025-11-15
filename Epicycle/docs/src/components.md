# Components

The Epicycle ecosystem consists of eleven specialized packages organized in architectural layers. Each package provides focused functionality while maintaining clean interfaces for composition into mission analysis workflows.

## Foundation Packages

Core abstractions and fundamental representations that underpin all mission analysis capabilities.

### [EpicycleBase](https://genastro.github.io/Epicycle/EpicycleBase/dev/)
Core abstractions, type hierarchy, and fundamental constants. Provides the foundational types and interfaces used throughout the ecosystem.

### [AstroStates](https://genastro.github.io/Epicycle/AstroStates/dev/)
Spacecraft state representations and coordinate transformations. Handles Cartesian, Keplerian, and Modified Equinoctial element representations with conversion capabilities.

### [AstroEpochs](https://genastro.github.io/Epicycle/AstroEpochs/dev/)
Time systems, epoch handling, and temporal conversions. Provides rigorous time standards including UTC, TAI, TT, and TDB with high-precision conversions.

## Physical Models

Environmental and spacecraft modeling capabilities for realistic mission simulation.

### [AstroUniverse](https://genastro.github.io/Epicycle/AstroUniverse/dev/)
Celestial body models, ephemeris data, and gravitational parameters. Includes planetary data, gravitational constants, and ephemeris interfaces.

### [AstroFrames](https://genastro.github.io/Epicycle/AstroFrames/dev/)
Reference frames and coordinate system transformations. Handles inertial and rotating frames with precise transformation matrices.

### [AstroModels](https://genastro.github.io/Epicycle/AstroModels/dev/)
Spacecraft and physical models for environmental interactions. Includes drag, solar radiation pressure, and other perturbation models.

## Mission Design

High-level capabilities for trajectory design, optimization, and mission analysis.

### [AstroManeuvers](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)
Maneuver models and algorithms for trajectory modification. Provides impulsive and finite-burn maneuver representations.

### [AstroCallbacks](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)
Mission analysis utilities for constraints, objectives, and events. Includes stopping conditions, cost functions, and constraint definitions.

### [AstroProp](https://genastro.github.io/Epicycle/AstroProp/dev/)
Numerical integration and trajectory propagation methods. High-order Runge-Kutta integrators with adaptive stepping and event detection.

### [AstroSolve](https://genastro.github.io/Epicycle/AstroSolve/dev/)
Optimization algorithms and constraint solving capabilities. Interfaces with IPOPT and other solvers for trajectory optimization problems.