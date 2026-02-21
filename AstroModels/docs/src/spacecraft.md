# Spacecraft Overview

The `Spacecraft` type represents a spacecraft with orbital state, time, physical properties, and trajectory history.

## Spacecraft Fields

- `state::OrbitState` — Orbital state (position and velocity) - see [State](spacecraft_state.md)
- `time::Time` — Epoch - see [Time](spacecraft_time.md)
- `mass::T` — Total mass (kg) - see [Mass](spacecraft_mass.md)
- `coord_sys::CoordinateSystem` — Coordinate system (origin and axes) - see [Coordinate System](spacecraft_coord_sys.md)
- `history::SpacecraftHistory` — Trajectory history - see [History](history.md)
- `name::String` — User label
- `cad_model::CADModel` — 3D model for visualization - see [CAD Model](spacecraft_cad_model.md)

## Basic Construction

```julia
using AstroModels, AstroStates, AstroEpochs

# Create spacecraft with all key properties
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = 1000.0,
    coord_sys = CoordinateSystem(earth, ICRFAxes()),
    name = "MySat"
)
```

See the `Spacecraft` constructor documentation for default values and additional options.

## Type Promotion

Spacecraft automatically promotes numeric types for automatic differentiation:

```julia
using ForwardDiff

# Mass with dual number
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    mass = ForwardDiff.Dual(1000.0, 1.0)
)

# All numeric fields promote to Dual
sc.state  # OrbitState with Dual elements
sc.mass   # Dual number
sc.time   # Time with Dual jd1, jd2
```

The numeric type `T` is determined by: `T = promote_type(eltype(state), typeof(time.jd1), typeof(mass))`


### Deep Copy

```julia
sc_copy = deepcopy(sc)
# All mutable fields (state, time, history) are independently copied
```
