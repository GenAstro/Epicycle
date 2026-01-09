# For basic language usage style use the idiomatic julia style guidelines here

https://docs.julialang.org/en/v1/manual/style-guide/

julia/base/rational.jl at 788b2c77c10c2160f4794a4d4b6b81a95a90940c · JuliaLang/julia

Some important ones 
make low level model code look like a math spec (see cart_to_kep.jl)
Avoid writing overly-specific types
Append ! to names of functions that modify their arguments
Avoid confusion about whether something is an instance or a type

# Function Style Guide and Example 

- https://docs.julialang.org/en/v1/manual/style-guide/

- use mathematical symbols as appropriate inside functions, but not in the interface 

- For inputs that must be differentiable, annotate them using a derived type from Real.  For example state::Vector{<:Real} so ForwardDiff (and other AD tools) can supply Vector{Dual<:Real} seamlessly. Avoid types that are too specific in specific in function interfaces unless they will never change.  The use of Real ensures some type checking, but uses dispatch for fast execution.  

``` julia 

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
cart = kep_to_cart([7000.0, 0.01, pi/4,0.0,0.0,pi/3], 398600.4418)

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

Note that the type definition state::Vector{<:Real} is what allows the use of forward mode differentiation here. 

```julia
using LinearAlgebra
using ForwardDiff
using FiniteDiff
using AstroStates  # adjust if module name differs

μ = 398600.4418  # km^3/s^2 (example Earth μ)
deg2rad(θ) = θ * (pi / 180)

# Keplerian test vector: [a, e, i, Ω, ω, ν]
# Moderate inclination, small eccentricity, all angles away from singularities
kep = [
    7000.0,                 # a (km)
    0.01,                   # e
    deg2rad(28.0),          # i
    deg2rad(40.0),          # Ω
    deg2rad(15.0),          # ω
    deg2rad(60.0)           # ν
]

# Function mapping Keplerian → Cartesian (6 → 6)
f(x) = AstroStates.kep_to_cart(x, μ)

# ForwardDiff Jacobian (6x6)
J_ad = ForwardDiff.jacobian(f, kep)
```

# Struct and API Design Principles for Julia Scientific/Optimization Code

## Summary Table

| Principle                | Description                                                                 |
|--------------------------|-----------------------------------------------------------------------------|
| Parametric types         | Use for performance, type stability, and AD support                         |
| Generic numeric fields   | Use `F<:Real` for AD compatibility                                          |
| Keyword constructors     | All fields as kwargs with defaults for user-facing structs                   |
| Docstrings               | Describe all fields, kwargs, and provide examples                           |
| Custom show methods      | Implement for clean REPL output                                             |
| Abstract/parametric fields| Use for extensibility and flexibility                                      |
| Validation               | Check field values and

## 1. Performance and Type Stability
- Use parametric types for numeric fields (e.g., `{F<:Real}`) when performance, type stability, or AD support is important.
- Ensure all related numeric fields use the same type parameter for consistency and efficiency.

## 2. Automatic Differentiation (AD) Compatibility
- Use generic numeric types (`F<:Real`) for fields that may participate in differentiation.
- Avoid hardcoding types like `Float64` for fields that may be used with dual numbers or other AD types.

## 3. User-Facing Constructors
- Provide an outer constructor that accepts all fields as keyword arguments (`kwargs`), with sensible defaults for each user facing struct. 
- This enables expressive, order-independent struct creation and makes APIs easier to use and maintain.

## 4. Documentation
- Write clear, idiomatic docstrings for all structs and constructors.
    - List and describe all fields and keyword arguments.
    - Include at least one usage example.
    - All fields must have range spec or description
- Document which fields are required and which are optional (with defaults).
- When constructors are overloaded, document the core constructor in full detail, then only document the difference in overloaded outer constructors

## 4.1 Validation
- Constructors must validate user inputs against the range spec in the doc string

## 5. REPL and Display Integration
- Implement `Base.show` methods for custom structs to provide concise, readable REPL output.
- Optionally, provide a verbose `show(io, ::MIME"text/plain", obj)` for detailed inspection.

## 6. Extensibility and Flexibility
- Use abstract types for fields that may be extended by users (e.g., `AbstractSRPModel`).
- Prefer `AbstractVector{F}` or parametric vectors for fields that may hold different numeric types.

## 7. Error Checking and Validation
- Validate individual field constraints when fields are set (e.g., check for correct vector lengths, positive-definite parameters, etc.)
- Validate field coupling and relationships at execution time (e.g., in `display()`, `propagate()`, `solve()`)
- Provide informative error messages for invalid input

## 8. Constructor Patterns for Composition Structs

Different struct types require different constructor patterns based on their usage:

### Rule 1: Required Positional Arguments
Use positional constructors when all inputs are required and there are no reasonable defaults.
- All arguments are required for valid construction
- No sensible default values exist
- Optional: May provide keyword argument interface for UI consistency

**Example:**
```julia
# Forces require specific configuration
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
integ = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9)
prop = OrbitPropagator(forces, integ)
```

### Rule 2: Keyword Arguments with Defaults
For structs with many fields that have reasonable defaults, provide defaults in the constructor to allow users to set only what they need.
- Most or all fields have sensible defaults
- Users customize only what they need
- Enables quick prototyping and clear intent

**Example:**
```julia
# Most maneuver parameters have standard defaults
man = ImpulsiveManeuver(; axes=VNB(), g0=9.81, Isp=300.0, element1=0.0, element2=0.0, element3=0.0)

# User sets only what differs from defaults
man = ImpulsiveManeuver(element1=1.3, axes=VNB())
```

### Rule 3: Defaults with Incremental Composition (add! methods)
For user-composed structs, use a constructor with reasonable defaults (possibly with kwargs) to create the base struct, then use `add!` methods for users to incrementally compose the struct.
- Struct is built incrementally over multiple steps
- All fields have reasonable defaults (like GUI-based tools: GMAT, STK, FreeFlyer)
- Complex validation between parts happens at execution time
- May have ordering or dependency requirements between composed elements

**Example:**
```julia
# Sequence starts with reasonable defaults
seq = Sequence()  
add_events!(seq, event1, Event[])
add_events!(seq, event2, [event1])

# View3D provides default coordinate system
view = View3D()  # Uses Earth GCRF by default
# Or user can override
view = View3D(coord_sys=moon_gcrf)
# Then compose incrementally
add_spacecraft!(view, sc)
display(view)  # Final validation happens here
```

**Key Principles:**
- Field-level validation in setters/add! methods (e.g., empty checks, bounds)
- Coupling validation at execution time (e.g., coord_sys matching between view and spacecraft)
- Mutable structs allow modification after construction

---

## Example: ImpulsiveManeuver Struct

```julia
"""
    ImpulsiveManeuver(; axes=:VNB, g0=9.81, Isp=300.0, element1=0.0, element2=0.0, element3=0.0)

Represents an impulsive maneuver in a specified reference frame.

# Keyword Arguments
- `axes::Symbol = :VNB`: Reference frame for the maneuver (e.g., :VNB, :LVLH, :ICRF).
- `g0::Real = 9.81`: Standard gravity [m/s²].
- `Isp::Real = 300.0`: Specific impulse [s].
- `element1::Real = 0.0`: First maneuver component (e.g., ΔV₁).
- `element2::Real = 0.0`: Second maneuver component (e.g., ΔV₂).
- `element3::Real = 0.0`: Third maneuver component (e.g., ΔV₃).

# Example
```julia

mutable struct ImpulsiveManeuver{F<:Real}
    axes::Symbol
    g0::F
    Isp::F
    element1::F
    element2::F
    element3::F

    function ImpulsiveManeuver(axes::Symbol, g0::F, Isp::F, element1::F, element2::F, element3::F) where {F<:Real}
        new{F}(axes, g0, Isp, element1, element2, element3)
    end
end

function ImpulsiveManeuver(;
    axes::Symbol = :VNB,
    g0::Real = 9.81,
    Isp::Real = 300.0,
    element1::Real = 0.0,
    element2::Real = 0.0,
    element3::Real = 0.0
)
    F = promote_type(typeof(g0), typeof(Isp), typeof(element1), typeof(element2), typeof(element3))
    new{F}(axes, F(g0), F(Isp), F(element1), F(element2), F(element3))
end

function Base.show(io::IO, m::ImpulsiveManeuver)
    println(io, "ImpulsiveManeuver(")
    println(io, "  axes = ", m.axes)
    println(io, "  g0 = ", m.g0)
    println(io, "  Isp = ", m.Isp)
    println(io, "  element1 = ", m.element1)
    println(io, "  element2 = ", m.element2)
    println(io, "  element3 = ", m.element3, ")")
end

man = ImpulsiveManeuver(element1=1.3, axes=:VNB)

## Differentiable Array Constructors in Julia

When writing code that must be compatible with automatic differentiation (AD) in Julia (e.g., ForwardDiff), always construct arrays (such as identity matrices, zeros, or ones) using the element type of your input data. This ensures that arrays will work with dual numbers or other AD types.

---

### Identity Matrix

```julia
I3 = Matrix{eltype(x)}(I, 3, 3)
z = zeros(eltype(x), n)
o = ones(eltype(x), n)

function foo(x)
    I3 = Matrix{eltype(x)}(I, 3, 3)
    z = zeros(eltype(x), 3)
    o = ones(eltype(x), 3)
    # ... use I3, z, o in AD-safe computations
end
---

For generic numeric code: write formulas with plain literals (2π, 0.5, 3/2) and rely on promotion; only use zero(T) / one(T) to seed accumulators or typed containers, and allocate results as Vector{T}. Avoid wrapping every constant in T(...)—clarity beats ceremonial generality.

Use state::Vector{<:Real} and μ::Float64 for differentiable conversion/math functions (e.g., orbital element ↔ Cartesian) to keep the signature explicit, simple, and AD‑compatible (ForwardDiff supplies Vector{Dual<:Real} which matches). Don’t generalize to AbstractVector or μ::Real until a concrete need (views, StaticArrays, μ differentiation) arises—add overloads later instead of widening the original. Preserve the element type in outputs by building a literal vector whose entries derive from the inputs. Treat physical degeneracies (parabolic/singular) with a warning and a NaN vector; treat shape/argument misuse with an error. Avoid unnecessary casts; plain literals (1.0, 2π) promote correctly with Duals and BigFloat. Only broaden the interface in response to an observed requirement, not preemptively.

# Docstring format

"""
    FUNCTION SIGNATURE

DECLARATIVE ON SENTENCE SUMMARY

"Arguments" (for a function) or "Fields" for a struct
- List of items

# Notes:
   # Do not include how, that may change, just document the contract. 
# Returns (for a function)

# Examples
```jldoctest
Code here (no julia prompts, want to be able to copy and past the whole chunk!!)

# output
REPL output goes here
```
"""

# Show methods should overload to support future JSON, HTML etc. 

function show(io::IO, ::MIME"text/plain", body::CelestialBody)
    println(io, "CelestialBody: ")
    println(io, "  name               = ", body.name)
    println(io, "  μ                  = ", body.mu)
    println(io, "  Equatorial Radius  = ", body.equatorial_radius)
    println(io, "  Flattening         = ", body.flattening)
    println(io, "  NAIF ID            = ", body.naifid)
end

# Delegate the generic show to the text/plain variant for print/println.
function show(io::IO, body::CelestialBody)
    show(io, MIME"text/plain"(), body)
end

---

# Struct Development Checklist

Use this checklist when creating new structs to ensure production-ready, well-tested code. This assumes the struct has been prototyped and the design has been reviewed.

## 1. Design Phase
- [ ] **Struct name** - PascalCase, descriptive, follows Epicycle conventions (see Naming Conventions)
- [ ] **Field names** - snake_case for multi-word, consistent with ecosystem patterns
- [ ] **Field types** - Use parametric types (`F<:Real`) for numeric fields that may be differentiated
- [ ] **Constructor pattern** - Choose appropriate pattern (positional, kwargs with defaults, or incremental composition)
- [ ] **Mutability** - Choose `struct` (immutable) or `mutable struct` based on usage pattern

## 2. Documentation Phase
- [ ] **Struct docstring** - Describes purpose and usage with signature showing constructor
- [ ] **Field documentation** - Each field documented with type and valid range/constraints
- [ ] **Examples** - At least one working example in docstring (use `jldoctest` when appropriate)
- [ ] **Constructor variants** - Document all public constructor signatures
- [ ] **Range specifications** - All numeric fields have documented valid ranges (e.g., "must be positive", "0.0 to 1.0")

## 3. Implementation Phase
- [ ] **Struct definition** - Inner constructor if needed for validation or type promotion
- [ ] **Outer constructors** - Keyword constructor with defaults for user-facing structs
- [ ] **Validation logic** - Constructor validates inputs against documented range specifications
- [ ] **Error messages** - Clear, actionable `ArgumentError` messages cite field name and constraint
- [ ] **Base.show method** - Implements clean REPL display (use `MIME"text/plain"` variant)
- [ ] **Base.deepcopy** - Implement if struct contains mutable fields or is itself mutable
- [ ] **Exports** - Add to package module exports list

## 4. Testing Phase
- [ ] **Happy path tests** - Verify default constructor and valid field combinations work
- [ ] **Validation tests** - Test each constraint with `@test_throws ArgumentError`
- [ ] **Edge cases** - Test boundary values (zero, negative, empty strings, etc.)
- [ ] **Error messages** - Verify error messages contain expected substrings
- [ ] **Type promotion** - Test with different numeric types if parametric (Float64, Int, ForwardDiff.Dual)
- [ ] **Show method** - Test REPL output format matches expected structure
- [ ] **Integration tests** - Test struct usage in broader package context (if applicable)

## 5. Integration Phase
- [ ] **Package documentation** - Add to package README or docs with system-level examples
- [ ] **Epicycle docs** - Add examples to main Epicycle documentation if appropriate
- [ ] **Update related code** - Modify dependent code that uses or composes this struct
- [ ] **Update tests** - Update integration tests in dependent packages
- [ ] **Code review** - Have design and implementation reviewed by another developer

## 6. Quality Checks
- [ ] **Naming consistency** - Compare field names with similar structs in ecosystem
- [ ] **No typos** - Run spell check on docstrings and comments
- [ ] **Type stability** - Verify `@code_warntype` shows no type instabilities (for performance-critical code)
- [ ] **Test coverage** - All code paths exercised by tests
- [ ] **Documentation builds** - Verify package docs build without warnings

## Quick Reference: Constructor Patterns

**Pattern 1 - Required Positional:**
```julia
# All fields required, no defaults
ForceModel(gravity::Gravity, drag::Drag)
```

**Pattern 2 - Keyword with Defaults:**
```julia
# Most fields have defaults, users customize as needed
ImpulsiveManeuver(; axes=VNB(), g0=9.81, Isp=300.0, element1=0.0)
```

**Pattern 3 - Defaults + Incremental Composition:**
```julia
# Start with defaults, build incrementally with add! methods
seq = Sequence()
add_events!(seq, event1, Event[])
```

See "Constructor Patterns for Composition Structs" section for detailed guidance.

# Naming Conventions

Consistent naming conventions improve code readability, maintainability, and reduce cognitive load across the Epicycle ecosystem.

## General Principles

1. **Clarity over brevity** - Use full words rather than abbreviations unless the abbreviation is universally understood in the domain (e.g., `mu` for gravitational parameter, `SRP` for solar radiation pressure)
2. **snake_case for multi-word fields** - Use underscores to separate words in field names for readability
3. **Lowercase for simple field names** - Single-word fields use lowercase (e.g., `state`, `time`, `mass`, `name`)
4. **Consistent terminology** - Use the same term for the same concept across all structs

## Field Naming Patterns

### Multi-word Fields
Use `snake_case` (underscores) to separate words:
```julia
coord_sys          # Not coordSys or coordsys
lower_bound        # Not lowerBound or lowerbound
upper_bound        # Not upperBound or upperbound
file_path          # Not filePath or filepath
model_scale        # Not modelScale
```

### Common Field Names
Standardize on these names for consistency across structs:

| Concept | Field Name | Example | Notes |
|---------|-----------|---------|-------|
| Name/identifier | `name` | `"MySpacecraft"` | String identifier for user reference |
| File path | `file_path` | `"assets/model.obj"` | Full or relative path to a file |
| Coordinate system | `coord_sys` | `CoordinateSystem(earth, ICRFAxes())` | Abbreviated but widely used |
| Scale factor | `scale` | `1.0` | Dimensionless scaling parameter |
| Visibility flag | `visible` | `true` | Boolean for display control |
| Time/epoch | `time` | `Time(...)` | Temporal point or epoch |
| State vector | `state` | `CartesianState(...)` | Orbital or physical state |
| Mass | `mass` | `1000.0` | Total mass in kg |
| Lower bound | `lower_bound` | `[-5.0, -2.0, 0.0]` | Optimization lower bounds |
| Upper bound | `upper_bound` | `[5.0, 2.0, 0.0]` | Optimization upper bounds |

### Boolean Fields
- Prefer simple descriptive names: `visible`, `active`, `enabled`
- Use `is_` prefix only when clarity requires it: `is_converged`, `is_valid`
- Avoid negative booleans: use `enabled` not `disabled`

Examples:
```julia
visible::Bool      # Clear without prefix
active::Bool       # Clear without prefix
is_converged::Bool # Prefix adds clarity for state checks
```

### Abbreviations
Acceptable abbreviations (use sparingly and consistently):
- `coord_sys` - coordinate system (widely used)
- `mu` or `μ` - gravitational parameter (standard in astrodynamics)
- `Isp` - specific impulse (standard)
- `SRP` - solar radiation pressure (standard)
- `num` prefix - for counts (e.g., `numvars`, but prefer `num_vars` for new code)

Avoid abbreviations unless words have many characters:
- `path` (not `pth`)
- `model` (not `mdl`)
- `scale` (not `scl`)
- `visible` (not `vis`)

Avoid: `file`, `filename`, `path` alone (ambiguous)

### Numeric Suffixes
When fields represent ordered elements or components:
```julia
element1, element2, element3  # Maneuver components
```

For multiple similar items, prefer vectors or collections:
```julia
tanks::Vector{Tank}           # Not tank1, tank2, tank3
thrusters::Vector{Thruster}   # Not thruster1, thruster2
```

## Struct Naming

- Use **PascalCase** for struct names: `Spacecraft`, `ImpulsiveManeuver`, `CoordinateSystem`
- Be descriptive: `CartesianState` not `CartState`, `OrbitPropagator` not `Propagator`
- Use full words: `Integrator` not `Integ` (exceptions: widely understood abbreviations)

## Function Naming

Follow Julia conventions:
- **lowercase with underscores**: `kep_to_cart`, `add_spacecraft!`, `to_posvel`
- **Append `!` for mutating functions**: `set_posvel!`, `add_events!`, `push_segment!`
- **Verb-noun order**: Use `verb_noun` pattern for clarity and consistency
  - ✓ **Correct**: `solve_trajectory`, `display_view`, `propagate_orbit`, `compute_state`
  - ✗ **Avoid**: `trajectory_solve`, `view_display`, `orbit_propagate`
  - The action (verb) should come first, followed by the object (noun)
  - Exception: conversion functions use `to_` prefix: `kep_to_cart`, `cart_to_kep`

### Legacy Naming

Some existing functions use noun_verb order and will be updated in future refactors:
- `trajectory_solve` → should be `solve_trajectory` (TODO: refactor with alias for backward compatibility)

## Internal/Private Fields

Prefix with underscore to indicate internal use:
```julia
_scene::Union{Nothing, Scene}  # Internal graphics state
_cache::Dict                    # Internal cached data
```

## Inspecting Existing Naming Patterns

Before creating new structs or fields, review existing naming patterns for consistency. Use this snippet to inspect field names:

```julia
using Epicycle

# Add structs you want to inspect
structs_to_inspect = [
    Spacecraft, ImpulsiveManeuver, Event, SolverVariable, 
    Sequence, CoordinateSystem, CelestialBody, ForceModel
]

for T in structs_to_inspect
    println("\n", T, ":")
    if !isabstracttype(T) && fieldcount(T) > 0
        for fname in fieldnames(T)
            println("  ", fname)
        end
    end
end
```

This helps you identify existing patterns and maintain consistency across the codebase.

## Migration Notes

Existing inconsistencies to address in future refactors:
- `numvars` → consider `num_vars` for consistency with `lower_bound`, `upper_bound`
- Review any abbreviated field names against the clarity principle

## Examples

### Good Field Names
```julia
struct Spacecraft
    state::OrbitState
    time::Time
    mass::Float64
    name::String
    coord_sys::CoordinateSystem
    history::Vector{Segment}
end

struct CADModel
    file_path::String
    scale::Float64
    visible::Bool
end

struct SRPModel
    area::Float64
    reflectivity::Float64
    file_path::String  # For plate model file
end
```

### Avoid
```julia
struct BadExample
    st::OrbitState          # Too abbreviated
    coordSys::CoordSys      # camelCase inconsistent
    fileName::String        # camelCase
    vis::Bool               # Too abbreviated
    lowerbound::Vector      # Missing underscore
end
```