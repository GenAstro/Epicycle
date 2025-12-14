```@meta
CurrentModule = AstroModels
```

# AstroModels

The AstroModels module provides physical models for astrodynamics applications, including spacecraft representations with state, time, mass and other properties.  The module supports multiple initialization patterns for spacecraft objects, allowing users to specify orbital states directly or construct them from component data. All models are designed for compatibility with automatic differentiation libraries.

## Quick Start

Create a spacecraft with orbital state and time information:

```julia
using AstroModels, AstroStates, AstroEpochs

# Method 1: Using a CartesianState struct
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0
)

# Method 2: Direct construction with OrbitState
sc = Spacecraft(
    state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03], Cartesian()),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0
)
```
**Note**: The spacecraft uses the `OrbitState` struct internally but accepts concrete state types or an `OrbitState` struct at construction. See the Reference Guide for complete details. 

Specifying coordinate system:

```julia
sc = Spacecraft(
    coord_sys = CoordinateSystem(mars, ICRFAxes())
)
```

Adding a 3D CAD model for visualization:

```julia
sc = Spacecraft(
    cad_model = CADModel(
        file_path = "path/to/model.obj",
        scale = 100.0,
        visible = true
    )
)
``` 

# Reference Guide

```@autodocs
Modules = [AstroModels]
```
# Index

```@index
```

