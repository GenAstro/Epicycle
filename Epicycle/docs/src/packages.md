# Packages

The Epicycle ecosystem is a set of packages organized in architectural layers. Each package provides focused functionality while maintaining clean interfaces for composition into mission analysis workflows.

| Layer | Packages |
|:------|:---------|
| Foundation: core abstractions and fundamental representations that underpin all mission analysis capabilities | [EpicycleBase](#EpicycleBase), [AstroRoutines](#AstroRoutines), [AstroStates](#AstroStates), [AstroEpochs](#AstroEpochs) |
| Physical models: environmental and spacecraft modeling capabilities for realistic mission simulation | [AstroUniverse](#AstroUniverse), [AstroFrames](#AstroFrames), [AstroModels](#AstroModels) |
| Mission design and navigation: high-level capabilities for trajectory design, optimization, and orbit determination | [AstroManeuvers](#AstroManeuvers), [AstroCallbacks](#AstroCallbacks), [AstroProp](#AstroProp), [AstroSolve](#AstroSolve) |
| Visualization and reporting: output for results produced by the packages above | [EpicycleIO](#EpicycleIO) |

Each package has its own documentation, linked from its section below.

## [EpicycleBase](https://genastro.github.io/Epicycle/EpicycleBase/dev/)

Core abstractions, type hierarchy, and fundamental constants. Provides the foundational types and interfaces used throughout the ecosystem.

## [AstroRoutines](https://genastro.github.io/Epicycle/AstroRoutines/dev/)

Low-level routines for classical astrodynamics.

## [AstroStates](https://genastro.github.io/Epicycle/AstroStates/dev/)

Spacecraft state representations such as Cartesian, Keplerian, Modified Equinoctial and many others.

## [AstroEpochs](https://genastro.github.io/Epicycle/AstroEpochs/dev/)

Time systems, epoch handling, and temporal conversions including UTC, TAI, TT, TDB and others with sub-nanosecond precision.

## [AstroUniverse](https://genastro.github.io/Epicycle/AstroUniverse/dev/)

Celestial body models, ephemeris data, and gravitational parameters.

## [AstroFrames](https://genastro.github.io/Epicycle/AstroFrames/dev/)

Reference frames and coordinate frame transformations.

## [AstroModels](https://genastro.github.io/Epicycle/AstroModels/dev/)

Physical models such as Spacecraft and GroundStations.

## [AstroManeuvers](https://genastro.github.io/Epicycle/AstroManeuvers/dev/)

Maneuver models.

## [AstroCallbacks](https://genastro.github.io/Epicycle/AstroCallbacks/dev/)

Utilities and interfaces to define quantities for constraints, objectives, events, and solve-fors.

## [AstroProp](https://genastro.github.io/Epicycle/AstroProp/dev/)

Numerical integration and trajectory propagation methods.

## [AstroSolve](https://genastro.github.io/Epicycle/AstroSolve/dev/)

Targeting, optimal control, and orbit estimation algorithms.

## [EpicycleIO](https://genastro.github.io/Epicycle/EpicycleIO/dev/)

Plotting, reporting, and 3D graphics.
