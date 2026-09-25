# Mass

A spacecraft carries its total mass in kilograms. Read it with `total_mass`; the field itself is
private, and `sc.mass` throws. Mass changes only by applying a maneuver, so there is no setter.

## Basic Usage

```julia
using AstroModels

# Specify mass at construction
sc = Spacecraft(
    mass = 1500.0  # kg
)

# Read mass
current_mass = total_mass(sc)
```

## Type Promotion

Mass participates in automatic differentiation type promotion:

```julia
using AstroModels, ForwardDiff

# Mass with Dual number
sc = Spacecraft(
    mass = ForwardDiff.Dual(1000.0, 1.0)  # Value with derivative
)

# All numeric fields promote to Dual
total_mass(sc)   # Dual{Float64}
sc.state  # OrbitState with Dual{Float64} elements
```

The spacecraft's numeric type `T` is: `T = promote_type(eltype(state), typeof(time.jd1), typeof(mass))`.
The SRP coefficient is not part of it: `c_r` belongs to the SRP geometry attached through the
`srp` field, not to the spacecraft.
