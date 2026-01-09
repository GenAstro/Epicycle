```@meta
CurrentModule = AstroModels
```

# AstroModels

AstroModels provides physical models for astrodynamics applications, including spacecraft representations with orbital state, time, and physical properties. All models are designed for compatibility with automatic differentiation libraries.

## Quick Start

```julia
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse

# Create spacecraft with all key properties
sc = Spacecraft(
    # Orbital state (position and velocity)
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    
    # Epoch
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    
    # Mass (kg)
    mass = 1000.0,
    
    # Coordinate system
    coord_sys = CoordinateSystem(earth, ICRFAxes()),
    
    # 3D model for visualization
    cad_model = CADModel(
        file_path = "models/satellite.obj",
        scale = 100.0,
        visible = true
    ),
    
    # Name
    name = "MySat"
)

# History is automatically populated during propagation
# See Spacecraft > History for details
```

!!! warning "Spacecraft state data type"
    The spacecraft uses the `OrbitState` struct internally but accepts concrete state types (e.g., `CartesianState`, `KeplerianState`) at construction. See the [State](spacecraft_state.md) section for complete details.

## Learn More

- [State](spacecraft_state.md) - Orbital state representations
- [Time](spacecraft_time.md) - Epoch and time scales
- [Mass](spacecraft_mass.md) - Spacecraft mass
- [Coordinate System](spacecraft_coord_sys.md) - Reference frames
- [CAD Model](spacecraft_cad_model.md) - 3D visualization
- [History](history.md) - Trajectory data
- [Reference](reference.md) - Complete API documentation

