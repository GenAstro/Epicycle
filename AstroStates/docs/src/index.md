```@meta
CurrentModule = AstroStates
```

# AstroStates

The AstroStates module provides models, structs, utilities, and conversions for orbital state representations. A state representation is a set of quantities that uniquely define an orbit. Supported forms include Cartesian, Keplerian, Modified Equinoctial, and others.

The module offers multiple interfaces for transforming and storing states. Low‑level conversion functions (e.g., `cart_to_kep.jl`) can be used directly. A type system automatically provides concrete structs for each representation (e.g., `CartesianState`) and converts between all supported permutations. The `OrbitState` utility preserves type stability when the representation may change by storing the numeric state and a type tag in separate fields. The library supports automatic differentiation with ForwardDiff.jl and Zygote.jl.

AstroStates is tested against output from the General Mission Analysis Tool (GMAT) R2022a.

References:

- Vallado, D. A. (2013), Fundamentals of Astrodynamics and Applications, 4th ed., Microcosm Press / Springer. 
- GMAT Development Team (2022), General Mission Analysis Tool (GMAT) Mathematical Specification, Version R2022a, NASA Goddard Space Flight Center. 

## Quick Start

The example below illustrates how to create a state struct, perform conversions, inspect elements of a state, and how to view all supported types.  

```julia
using AstroStates

# Define a Cartesian state
cart = CartesianState([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5])

# Convert to Keplerian then back to Cartesian
mu = 398600.4418 
kep   = KeplerianState(cart, mu)     
cart2 = CartesianState(kep, mu)     

# Display some state elements
kep.sma
kep.raan

# Generate a vector containing the state struct data
to_vector(kep)

# See a list of all supported representations
subtypes(AbstractOrbitState)
```
The API is documented with docstrings; external references are intentionally omitted to avoid duplication. In the REPL, type `?` to enter help mode, then enter a name to view its documentation. For example, `?IncomingAsymptoteState` displays the incoming hyperbolic asymptote state, and `?cart_to_sphradec` shows the spherical RA/Dec conversion helper.

## State Overview

AstroStates provides a library of state structs to create, store, and convert orbit states.  These structs derive from AbstractOrbitState. You can create states from numeric vectors or by converting from another state struct; conversions are performed automatically via overloaded constructors. State structs print readably and expose elements as fields. To list supported representations, run `subtypes(AbstractOrbitState)`. Kinematic conversions (e.g., Cartesian, Spherical) do not require mu, while conic element conversions (e.g., Keplerian, Modified Equinoctial) do.

```julia
using AstroStates

# Create a Cartesian state from a vector
c = CartesianState([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5])

# Create a Keplerian state from individual elements.
k = KeplerianState(-98000.0, 2.6, pi/4, deg2rad(145), pi/8, 0.0 )

# Create a Keplerian state from a Cartesian State performing conversion automatically
mu = 398600.4415
k2 = KeplerianState(c, mu)

# Convert the Keplerian state to outgoing asymptote representation
h = OutGoingAsymptoteState(k, mu)

# Inspect elements of the states we just created.  Use "?" to see fields on a struct.
k.sma
c.posvel
h.c3
```

## OrbitState Container

The OrbitState struct provides a type-stable container when the representation may change during a simulation but type stability is still required. For example, Epicycle’s Spacecraft uses OrbitState to accept different input representations and to switch state types during a run. OrbitState stores (1) the state data and (2) a tag that describes the representation. The tag is a state-type marker that parallels the concrete state struct names (e.g., Keplerian, Cartesian, etc.).

```julia
using AstroStates

# Create an OrbitState struct that stores the state and state type.
os = OrbitState([-98000.0, 2.6, pi/4, deg2rad(145), pi/8, 0.0 ],Keplerian())

# Print the state and type
println(os.state)
println(os.statetype)

# Create an OrbitState struct from a concrete type struct.
c = CartesianState([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5])
os = OrbitState(c)

# See all available types
subtypes(AbstractOrbitStateType)
```
---

## Conversions Overview

The conversion functions in AstroStates are contained in individual files with function names like `cart_to_kep.jl`.  These functions can be used directly without the struct-based interfaces above when appropriate and when that is easier to integrate into other applications.  

```julia
using AstroStates

# Bypass structs and work directly with vectors, etc.  
mu = 398600.4415
k  = cart_to_kep([7000.0, 0.0, 100.0, 0.0, 7.5, 2.5], mu)

# Convert an equinoctial state to alternate equinoctial state
ae = equinoctial_to_alt_equinoctial([7758.763,-0.0047,0.09769,-0.00695,0.16227, 6.2762])
```
The conversions are written to resemble astrodynamics textbooks with the intention that the code can serve as its own math spec. Here is an example from `kep_to_cart.jl`:

``` julia
using LinearAlgebra

"""
    kep_to_cart(state::Vector{<:Real}, μ::Real; tol::Float64=1e-12)

Convert a Keplerian state vector to a Cartesian state vector.

# Arguments
- `state::Vector{<:Real}`: Keplerian elements `[a, e, i, Ω, ω, ν]`
- `μ`: Gravitational parameter
- `tol`: Tolerance for singularities like p ≈ 0 (default: 1e-12)
- `a`: semi-major axis
- `e`: eccentricity
- `i`: inclination
- `Ω`: right ascension of ascending node
- `ω`: argument of periapsis
- `ν`: true anomaly

# Returns
A 6-element vector `[x, y, z, vx, vy, vz]` representing Cartesian position and velocity.

# Example
cart = kep_to_cart([7000.0, 0.01, pi/4, 0.0, 0.0, pi/3], 398600.4418)

# Notes
- Angles must be in radians.
- Dimensional quantities must be consistent units with μ.
- Returns a vector of `NaN`s if conversion is undefined.
"""
function kep_to_cart(state::Vector{<:Real}, μ::Real; tol::Float64=1e-12)
    if length(state) != 6
        error("Input vector must have exactly six elements: a, e, i, Ω, ω, ν.")
    end

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(NaN, 6)
    end

    # Unpack the elements
    a, e, i, Ω, ω, ν = state

    # Compute semi-latus rectum: p = a * (1 - e²)
    p = a * (1.0 - e^2)

    # Check for degenerate orbit (e.g., parabolic or collapsed)
    if p < tol || abs(1-e) < tol
        @warn "Conversion Failed: Orbit is parabolic or singular."
        return fill(NaN, 6)
    end

    # Compute radial distance: r = p / (1 + e * cos(ν))
    r = p / (1.0 + e * cos(ν))

    # Position and velocity in perifocal frame 
    factor = sqrt(μ / p)
    r̄ₚ = [r * cos(ν), r * sin(ν), 0.0]
    v̄ₚ = [-factor * sin(ν), factor * (e + cos(ν)), 0.0]

    # Precompute sines and cosines for rotation matrix
    cos_Ω, sin_Ω = cos(Ω), sin(Ω)
    cos_ω, sin_ω = cos(ω), sin(ω)
    cos_i, sin_i = cos(i), sin(i)

    # Rotation matrix from perifocal to inertial
    R = [
        cos_ω * cos_Ω - sin_ω * cos_i * sin_Ω   -sin_ω * cos_Ω - cos_ω * cos_i * sin_Ω   sin_i * sin_Ω;
        cos_ω * sin_Ω + sin_ω * cos_i * cos_Ω   -sin_ω * sin_Ω + cos_ω * cos_i * cos_Ω  -sin_i * cos_Ω;
        sin_ω * sin_i                                    cos_ω * sin_i                   cos_i
    ]

    # Rotate position and velocity from perifocal to inertial frame
    pos = R * r̄ₚ
    vel = R * v̄ₚ 

    return vcat(pos, vel)
end
```
---

## Automatic Differentiation 

All functions and conversions in AstroStates are fully differentiable using Julia's automatic differentation libraries ForwardDiff and Zygote. Examples for computing Jacobians are shown below.  

Note: The time to precompile AD interfaces is substantial, but those times are only incurred on the first execution and when included in loops or functions the times are orders of magnitude faster. REPL peformance for these examples is poor for that reason. 

```julia
using ForwardDiff
using AstroStates

# Define the state vector and mu
x = [7000.0, 0.0, 100.0, 0.0, 7.5, 2.5]
mu = 398600.4418

# Define a function closure that returns a vector  
f(x) = to_vector(KeplerianState(CartesianState(x, mu), mu))

# Compute the Jacobian of Keplerian state w/r/t Cartesian State at x
J = ForwardDiff.jacobian(f, x)
```

```julia
using Zygote
using AstroStates

# State vector and mu
x = [7000.0, 0.0, 100.0, 0.0, 7.5, 2.5]
mu = 398600.4418

# Define a function closure that returns a vector  
f(x) = to_vector(ModifiedEquinoctialState(CartesianState(x, mu), mu))

# Compute the Jacobian of Modified Equinoctial elements w/r/t Cartesian
J = first(Zygote.jacobian(f, x))  
```

## State Types Reference

```@autodocs
Modules = [AstroStates]
Order = [:type]
```

## Conversions Reference

```@autodocs
Modules = [AstroStates]
Order = [:function]
```

## API Index

```@index
Pages = ["index.md"]
```

