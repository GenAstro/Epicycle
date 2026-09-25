# Spacecraft Overview

The `Spacecraft` type represents a spacecraft with orbital state, time, physical properties, and trajectory history.

## Spacecraft Fields

- `state::OrbitState` — Orbital state (position and velocity) - see [State](spacecraft_state.md)
- `time::Time` — Epoch - see [Time](spacecraft_time.md)
- `coord_sys::CoordinateSystem` — Coordinate system (origin and axes) - see [Coordinate System](spacecraft_coord_sys.md)
- `history::SpacecraftHistory` — Trajectory history - see [History](history.md)
- `name::String` — User label
- `cad_model::CADModel` — 3D model for visualization - see [CAD Model](spacecraft_cad_model.md)
- `drag` — Drag geometry such as `SphericalDrag`, or `nothing`
- `srp` — Solar radiation pressure geometry such as `SphericalSRP`, or `nothing`
- `save_history::Bool` — Whether propagation appends a segment to `history`

Total mass in kg is set at construction and read with `total_mass(sc)`; the field itself is private. See [Mass](spacecraft_mass.md).

## Basic Construction

```julia
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse

# Create spacecraft with all key properties
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0,
    coord_sys = CoordinateSystem(earth, ICRF()),
    name = "MySat"
)
```

See the `Spacecraft` constructor documentation for default values and additional options.

## Type Promotion

Spacecraft automatically promotes numeric types for automatic differentiation:

```julia
using AstroModels, AstroStates, AstroEpochs, ForwardDiff

# Mass with dual number
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = ForwardDiff.Dual(1000.0, 1.0)
)

# All numeric fields promote to Dual
sc.state  # OrbitState with Dual elements
total_mass(sc)   # Dual number
sc.time   # Time with Dual jd1, jd2
```

The numeric type `T` is determined by: `T = promote_type(eltype(state), typeof(time.jd1), typeof(mass))`


### Deep Copy

<!-- doc-continue -->
```julia
sc_copy = deepcopy(sc)
# All mutable fields (state, time, history) are independently copied
```
