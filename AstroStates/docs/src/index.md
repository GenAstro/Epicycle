```@meta
CurrentModule = AstroStates
```

# AstroStates

The AstroStates module provides models, structs, utilities, and conversions for orbital state representations. A state representation is a set of quantities that uniquely define an orbit. Supported forms include Cartesian, Keplerian, Modified Equinoctial, and others.

The module offers multiple interfaces for transforming and storing states. Low‑level conversion functions (e.g., `cart_to_kep.jl`) can be used directly. A type system automatically provides concrete structs for each representation (e.g., `CartesianState`) and converts between all supported permutations. The `OrbitState` utility preserves type stability when the representation may change by storing the numeric state and a type tag in separate fields. The conversions are differentiable with ForwardDiff.jl.

AstroStates is tested against output from the General Mission Analysis Tool (GMAT) R2022a.

References:

- Vallado, D. A. (2013), Fundamentals of Astrodynamics and Applications, 4th ed., Microcosm Press / Springer. 
- GMAT Development Team (2022), General Mission Analysis Tool (GMAT) Mathematical Specification, Version R2022a, NASA Goddard Space Flight Center. 

## Quick Start

The example below illustrates how to create a state struct, perform conversions, inspect elements of a state, and how to view all supported types.  

```julia
using AstroStates
using InteractiveUtils        # subtypes

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

# Create a Cartesian state from position and velocity vectors
c = CartesianState([7000.0, 0.0, 100.0] , [0.0, 7.5, 2.5])

# Create a Keplerian state from individual elements.
k = KeplerianState(-98000.0, 2.6, pi/4, deg2rad(145), pi/8, 0.0 )

# Create a Keplerian state from a Cartesian State performing conversion automatically
mu = 398600.4415
k2 = KeplerianState(c, mu)

# Convert the Keplerian state to outgoing asymptote representation
h = OutGoingAsymptoteState(k, mu)

# Inspect elements of the states we just created.  Use "?" to see fields on a struct.
c.position
c.velocity
k.sma
c.posvel
h.c3
```

## OrbitState Container

The OrbitState struct holds a state whose representation may change during a run, without the
cost that normally comes from letting a type vary. For example, Epicycle’s Spacecraft uses OrbitState to accept different input representations and to switch state types during a run. OrbitState stores (1) the state data and (2) a tag that describes the representation. The tag is a state-type marker that parallels the concrete state struct names (e.g., Keplerian, Cartesian, etc.).

```julia
using AstroStates
using InteractiveUtils        # subtypes

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
function kep_to_cart(state::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)
    if length(state) != 6
        error("Input vector must have exactly six elements: a, e, i, Ω, ω, ν.")
    end
    T = float(promote_type(eltype(state), typeof(μ)))

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(T(NaN), 6)
    end

    # Unpack the elements
    a, e, i, Ω, ω, ν = state

    # Semi-latus rectum: p = a * (1 - e²)
    p = a * (1 - e^2)

    # Degenerate orbit (parabolic or collapsed)
    if p < tol || abs(1 - e) < tol
        @warn "Conversion Failed: Orbit is parabolic or singular."
        return fill(T(NaN), 6)
    end

    # Radial distance: r = p / (1 + e cos ν). On a hyperbola the denominator reaches zero at the
    # asymptote; beyond it r would be negative and the state a mirror image of no real point.
    denom = 1 + e * cos(ν)
    if denom <= tol
        @warn "Conversion Failed: True anomaly $(ν) is at or beyond the asymptote of a hyperbola " *
              "with eccentricity $(e)."
        return fill(T(NaN), 6)
    end
    r = p / denom

    # Position and velocity in the perifocal frame
    factor = sqrt(μ / p)
    sν, cν = sincos(ν)
    r̄ₚ = SVector{3,T}(r * cν, r * sν, 0)
    v̄ₚ = SVector{3,T}(-factor * sν, factor * (e + cν), 0)

    # Rotation from perifocal to inertial
    sΩ, cΩ = sincos(Ω)
    sω, cω = sincos(ω)
    si, ci = sincos(i)
    R = SMatrix{3,3,T}(cω * cΩ - sω * ci * sΩ,  cω * sΩ + sω * ci * cΩ,  sω * si,
                       -sω * cΩ - cω * ci * sΩ, -sω * sΩ + cω * ci * cΩ, cω * si,
                       si * sΩ,                 -si * cΩ,                ci)      # column-major

    pos = R * r̄ₚ
    vel = R * v̄ₚ
    return T[pos[1], pos[2], pos[3], vel[1], vel[2], vel[3]]
end
```
---

## Automatic Differentiation 

The conversions are differentiable with ForwardDiff, including at periapsis, apoapsis and a zero node or argument of periapsis. The test suite checks this for every representation: the Jacobian of the conversion from Cartesian times the Jacobian of the conversion back is the identity. Derivatives do not exist at the states where an element is undefined, an exactly circular or exactly equatorial orbit, and there the returned derivative is not meaningful.

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

## State Types Reference

```@autodocs
Modules = [AstroStates]
Public  = true
Private = false
Order = [:type]
```

## Conversions Reference

```@autodocs
Modules = [AstroStates]
Public  = true
Private = false
Order = [:function]
```

## API Index

```@index
Pages = ["index.md"]
```

