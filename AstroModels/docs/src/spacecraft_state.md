# State

The `state` field stores the spacecraft's orbital state (position and velocity) using the `OrbitState` type from AstroStates.

## Creating Spacecraft with State

```julia
using AstroModels, AstroStates

# Method 1: Using CartesianState (recommended)
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03])
)

# Method 2: Using KeplerianState
sc = Spacecraft(
    state = KeplerianState([7000.0, 0.01, 45.0, 0.0, 0.0, 0.0])
)

# Method 3: Direct OrbitState construction
sc = Spacecraft(
    state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03], Cartesian())
)
```

**Note**: The spacecraft uses `OrbitState` internally but accepts concrete state types (e.g., `CartesianState`, `KeplerianState`) at construction for convenience.

## Accessing State

```julia
# Get current state in specific representation
cart_state = get_state(sc, Cartesian())
kep_state = get_state(sc, Keplerian())

# Get position/velocity vector (6-element)
pv = to_posvel(sc)  # [x, y, z, vx, vy, vz]
```

## Modifying State

```julia
# Update position and velocity
new_pv = [7050.0, 0.0, 0.0, 0.0, 7.6, 0.0]
set_posvel!(sc, new_pv)  # Mutates spacecraft in place
```

## Type Promotion

State components automatically promote for automatic differentiation:

```julia
using ForwardDiff

# Create with Dual mass
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    mass = ForwardDiff.Dual(1000.0, 1.0)
)

# All state components promoted to Dual
sc.state  # OrbitState with Dual{Float64} elements
```
