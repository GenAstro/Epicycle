# Mass

The `mass` field stores the spacecraft's total mass in kilograms.

## Basic Usage

```julia
# Specify mass at construction
sc = Spacecraft(
    mass = 1500.0  # kg
)

# Access mass
current_mass = sc.mass
```

## Type Promotion

Mass participates in automatic differentiation type promotion:

```julia
using ForwardDiff

# Mass with Dual number
sc = Spacecraft(
    mass = ForwardDiff.Dual(1000.0, 1.0)  # Value with derivative
)

# All numeric fields promote to Dual
sc.mass   # Dual{Float64}
sc.state  # OrbitState with Dual{Float64} elements
sc.c_r    # Dual{Float64}
```

The spacecraft's numeric type `T` is: `T = promote_type(eltype(state), typeof(time.jd1), typeof(mass), typeof(c_r))`
