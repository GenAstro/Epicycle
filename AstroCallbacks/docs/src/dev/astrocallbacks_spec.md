# AstroCallbacks – Spec (v0.1)

## 1. Overview and Scope (What are we building?)

**Purpose:**

AstroCallbacks provides a unified interface for defining and evaluating quantities of interest from astrodynamics domain objects (spacecraft, maneuvers, celestial bodies). Calcs enable reuse - define a quantity once, use it everywhere (optimization, plotting, propagation, reporting, event detection, estimation).

**In-scope:**
- Standard calc definition that is expressive and applies across ecosystem and to custom calcs
- Standardized access patterns: `get_calc()` and `set_calc!()`
- Type-safe variable definitions using Julia's type hierarchy
- Spacecraft state quantities (position, velocity, orbital elements)
- Time-based quantities (epochs, durations, in various scales/formats)
- Control quantities (maneuver parameters)
- Celestial body properties
- User-defined custom quantities
- Composable calculations (derive new quantities from existing ones)
- Calc attribute system (settability, continuity, differentiability, etc.)

**Out-of-scope:**
- How optimizers use calcs (that's AstroSolve)
- How plotting tools use calcs (external plotting packages)
- Report formatting and generation (reporting packages)
- Event sequencing logic (AstroSolve)
- Numerical optimization algorithms

## 2. Requirements (What we are building.)

### Functional Requirements - Calc Types

AstroCallbacks shall support calcs for:
- AstroState state elements (AstroStates.jl)
- AstroTime time elements (AstroEpochs.jl)
- Maneuver elements (AstroManeuvers.jl)
- CelestialBody properties (AstroUniverse.jl)
- Spacecraft properties (AstroModels.jl)
- ForceModel properties (AstroProp.jl)
- User-defined custom calcs

Note, this list is based on current functionality.  The design needs to be extendible as the system grows. 

### Functional Requirements - Calc Interface

AstroCallbacks shall provide interfaces for:
- Reading current values from domain objects (e.g., via a get method)
- Writing new values to mutable domain objects (e.g., via a set method, when calc is settable)

AstroCallbacks shall expose calc attributes to support integration:
- Is the calc numeric (returns numeric type vs non-numeric like String)
- Calc dimensionality (scalar or vector, with size)
- Is the calc settable (supports modification)
- Is the calc continuous (suitable for gradient-based optimization)
- Does the calc support automatic differentiation

### Integration Requirements

AstroCallbacks calcs shall support integration into:
- Propagation stopping conditions (e.g., `StopAt` using time calc)
- Optimization/estimation variables (solve-fors via `SolverVariable`)
- Optimization/estimation constraints (via `Constraint`)
- Plotting (time series evaluation)
- Reporting (formatted output, including non-numeric calcs)

AstroCallbacks shall ensure:
- Automatic differentiation (AD) compatibility for numeric mutable calcs
- Changes to calc values are visible to external code (e.g., event function closures)

### Extensibility Requirements

AstroCallbacks shall enable:
- Users to define custom variable types
- Users to define custom calc implementations (extend get and optionally set methods)
- Package developers to add new subject-variable combinations

### Design Requirements

AstroCallbacks shall ensure:
- Calcs can reference multiple domain objects (composite calculations)
- Calcs can derive new quantities from existing calcs
- Type-stable operations for compiler optimization
- Clear error messages for invalid operations (type mismatches, setting read-only calcs, out-of-range values)

## 3. Design (How we are building it.)

### Design Overview

The AstroCallbacks design centers on **calc objects** - lightweight wrappers that bind one or more domain objects (the "subject(s)") to a specific quantity of interest (the "variable"). This subject-variable binding pattern enables type-safe, extensible quantity access across the Epicycle ecosystem.

**Core Architecture:**
- **Calc Types** (`OrbitCalc`, `ManeuverCalc`, `BodyCalc`, `TimeCalc`) - categorized by subject type
- **Variable Types** (`SMA`, `Cartesian`, `DeltaVMag`, etc.) - define which quantity to access
- **Access Interface** (`get_calc()`, `set_calc!()`) - retrieve and modify values via multiple dispatch
- **Attribute System** - calcs self-describe capabilities (settability, continuity, dimensionality, AD compatibility)

**Key Design Principles:**
1. **Separation of concerns**: Calcs define *what* quantities are; external modules determine *how* to use them
2. **Type safety**: Compiler catches invalid subject-variable combinations at compile time
3. **Extensibility**: Users and packages add new calcs via Julia's multiple dispatch
4. **Composability**: Calcs can reference multiple objects and aggregate values from other calcs
5. **Mutable references**: Calcs hold references to mutable domain objects, enabling in-place updates during propagation and optimization
6. **Encapsulation**: Calcs interact with domain objects through public get/set methods, not direct field access

**Information Flow:**
```
Domain Object (Spacecraft, Maneuver, Time)
    ↓ wrapped by
Calc (OrbitCalc, ManeuverCalc, TimeCalc) + Variable Type (SMA, Cartesian, JD)
    ↓ accessed via
get_calc() / set_calc!()
    ↓ returns/accepts
Numeric or Non-numeric Value
    ↓ consumed by
External Modules (AstroSolve, plotting, reporting)
```

**Calc Patterns:**

The design supports three complementary patterns for defining calcs:

1. **Subject-Variable Binding** - Standard pattern for accessing properties of domain objects
   - Example: `OrbitCalc(spacecraft, SMA())` accesses semi-major axis
   - Use when: Accessing a well-defined quantity from a single subject type
   - Characteristics: Type-safe, extensible via dispatch, clear semantics

2. **Multi-Subject Calcs** - Access quantities involving multiple domain objects
   - Example: `RelativeMotionCalc(sat1, sat2, RelativeDistance())` 
   - Use when: Quantity depends on relationships between multiple objects
   - Characteristics: Same subject-variable binding pattern, multiple subjects

3. **Function-Based Calcs** - Arbitrary computations over domain objects
   - Example: Custom cost function, complex constraint, physics-based metric
   - Use when: Computation doesn't map cleanly to a variable type, or requires complex logic
   - Characteristics: Maximum flexibility, user-defined computation
   - Implementation: Custom calc struct holds function + object references, `get_calc()` evaluates function

All three patterns share the same public interface (`get_calc()`, `set_calc!()`, attribute functions) enabling uniform integration into optimization, plotting, and reporting workflows regardless of underlying implementation.

The following subsections detail the specific types, functions, patterns, and conventions that implement this design.

### Core Components

**Calc Types (the wrappers):**

Calc types are structs that hold references to domain objects and specify what quantity to access:

- `OrbitCalc` - accesses spacecraft state quantities (mutable)
- `ManeuverCalc` - accesses maneuver parameters (mutable)
- `BodyCalc` - accesses celestial body properties (typically read-only)
- `TimeCalc` - accesses time values with specified scale/format (mutable)
- Custom calc types - user-defined for multi-subject or composite calculations

**Variable Type Hierarchy (what quantities):**

Variable types specify *which* quantity a calc accesses. They are immutable singleton types organized in a hierarchy:

- `AbstractCalcVariable` - base type for all variable types
  - `AbstractOrbitVar` - base for orbit quantities
    - `SMA`, `Inc`, `RAAN`, `Cartesian`, etc. - concrete orbit variable types
  - `AbstractManeuverVar` - base for maneuver quantities
    - `DeltaVMag`, `DeltaVVector`, etc. - concrete maneuver variable types
  - `AbstractBodyVar` - base for body quantities
    - `GravParam`, `Radius`, etc. - concrete body variable types
  - User-defined variable types for custom quantities

**Relationship:**

Calc types combine one or more subjects with a variable type: `CalcType(subject(s)..., VariableType())`

Example: `OrbitCalc(spacecraft, SMA())` creates a calc that accesses the SMA of the spacecraft.

**Functions:**

**Core Access:**
- `get_calc(calc)` - retrieve current value from calc
- `set_calc!(calc, value)` - update mutable calc (! indicates mutation)

**Calc Attributes:**
- `is_settable(calc)` - returns true if calc supports modification
- `is_numeric(calc)` - returns true if calc returns numeric type
- `length(calc)` - returns dimensionality (1 for scalar, n for vector)

**Extensibility:**
- `make_calc(subject(s), variable)` - dispatch point for creating calcs (extensibility)

**Patterns:**
- **Subject-Variable Binding**: `CalcType(subject(s), variable_type)`
  - Example: `OrbitCalc(spacecraft, SMA())`, `ManeuverCalc(maneuver, DeltaVVector())`
- **Direct Object Reference**: `CalcType(object, ...params)`
  - Example: `TimeCalc(sat.time, TDB(), JD())`
- **Composite Calcs**: Custom struct holds multiple references, aggregates in `get_calc()`

### Usage Examples

**Use Case: Full Spacecraft State and Core API**

Demonstrates complete calc API: creation, get, set, and attribute queries.

```julia
# Create spacecraft
sat = Spacecraft(
    state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
    time = Time("2025-01-01T00:00:00", UTC(), ISOT())
)

# Create calc for full Cartesian state
state_calc = OrbitCalc(sat, IncomingAsymptote())

# Read current value
current_state = get_calc(state_calc)  
# Returns: [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0] (6-element vector [x,y,z,vx,vy,vz])

# Modify value (if settable)
new_state = [7100.0, 100.0, 50.0, 0.1, 7.6, 0.2]
set_calc!(state_calc, new_state)
# Spacecraft state updated in-place

# Query calc attributes
is_settable(state_calc)  # Returns: true (can modify spacecraft state)
is_numeric(state_calc)   # Returns: true (returns Float64 vector)
length(state_calc)       # Returns: 6 (dimensionality: [x,y,z,vx,vy,vz])

# Other element types
state_calc = OrbitCalc(sat, IncomingAsymptote())
```
---

**Use Case: Spacecraft State Element**

Access individual orbital elements (scalar or vector quantities).

```julia
# Scalar orbital element
sma_calc = OrbitCalc(sat, SMA())         # Semi-major axis (km)
inc_calc = OrbitCalc(sat, Inc())         # Inclination (rad)
```
---

**Use Case: Time in Different Scales and Formats**

Time calcs support multiple time scales (UTC, TDB, TT) and formats (numeric vs string).

```julia
# Numeric formats (settable, for optimization/plotting)
jd_calc = TimeCalc(sat.time, TDB(), JD())        # Julian Date
mjd_calc = TimeCalc(sat.time, UTC(), MJD())      # Modified Julian Date  

# String format (read-only, for reporting)
iso_calc = TimeCalc(sat.time, UTC(), ISOT())           # ISO 8601 string
# get_calc(iso_calc) returns "2025-01-01T00:00:00"
```

---

**Use Case: Custom Composite Calc**

Aggregate values from multiple existing calcs.

```julia
# Define custom calc struct (holds references to multiple objects)
struct TotalDeltaVCalc
    maneuvers::Vector{ImpulsiveManeuver}
    spacecraft::Spacecraft
end

# Implement get_calc (aggregates from other calcs)
function AstroCallbacks.get_calc(calc::TotalDeltaVCalc)
    total = 0.0
    for mnv in calc.maneuvers
        dv_calc = ManeuverCalc(mnv, calc.spacecraft, DeltaVMag())
        total += get_calc(dv_calc)
    end
    return total
end

# Usage
toi = ImpulsiveManeuver(axes=VNB(), element1=2.5, element2=0.0, element3=0.0)
mcc = ImpulsiveManeuver(axes=VNB(), element1=0.6, element2=0.5, element3=0.0) 
moi = ImpulsiveManeuver(axes=VNB(), element1=0.3, element2=0.0, element3=0.0)

total_dv_calc = TotalDeltaVCalc([toi, mcc, moi], sat)
# get_calc(total_dv_calc) returns ~3.27 km/s
```

---

**Use Case: Maneuver Parameters (Multi-Subject Pattern)**

Access and modify impulsive maneuver delta-v components.

```julia
# Create maneuver in VNB frame
toi = ImpulsiveManeuver(axes=VNB(), element1=2.5, element2=0.1, element3=-0.2)

# Scalar maneuver magnitude (read-only, computed)
mag_calc = ManeuverCalc(toi, sat, DeltaVMag())
# get_calc(mag_calc) returns ~2.51 km/s

# Vector maneuver components (settable)
vec_calc = ManeuverCalc(toi, sat, DeltaVVector())
# get_calc(vec_calc) returns [2.5, 0.1, -0.2]
# set_calc!(vec_calc, [2.6, 0.0, 0.0]) updates maneuver

# Note: ManeuverCalc demonstrates Pattern 2 (Multi-Subject Calcs)
# Takes two subjects: maneuver and spacecraft (needed for frame transformations)
```

---

**Use Case: Celestial Body Properties**

Access celestial body properties (typically read-only).

```julia
# Body gravitational parameter
mu_calc = BodyCalc(earth, GravParam())
# get_calc(mu_calc) returns 398600.4418 (km³/s²)

# Body equatorial radius  
radius_calc = BodyCalc(earth, Radius())
# get_calc(radius_calc) returns 6378.137 (km)

# Body properties are typically read-only
# set_calc!(mu_calc, new_value) would throw error
```

---

**Use Case: Custom Variable and Vector Calcs**

Wrap user-defined parameters as calcs for optimization.

```julia
# Scalar parameter
stop_radius = VariableCalc(10000.0)
get_calc(stop_radius)           # Returns: 10000.0
set_calc!(stop_radius, 9500.0)  # Updates internal value
get_calc(stop_radius)           # Returns: 9500.0

# Vector parameter
gains = VariableCalc([0.1, 0.2, 0.3])
get_calc(gains)                    # Returns: [0.1, 0.2, 0.3]
set_calc!(gains, [0.15, 0.25, 0.35])

# Use in optimization
var_radius = SolverVariable(calc=stop_radius, lower=8000.0, upper=12000.0)
var_gains = SolverVariable(calc=gains, lower=[0.0, 0.0, 0.0], upper=[1.0, 1.0, 1.0])

# Use in functions (accepts calc or literal)
propagate!(prop, sat, StopAt(sat, PosMag(), stop_radius))
# StopAt evaluates: get_calc(stop_radius) during propagation

# Pattern: Any function accepting numeric values can accept calc
# Enables parameter optimization without changing function signatures
```

---

**Use Case (DEFERRED): Spacecraft Properties**

Design validation for spacecraft-level properties (mass, area, coefficients).

```julia
# Mass (changes with maneuvers)
mass_calc = SpacecraftCalc(sat, Mass())
# get_calc(mass_calc) returns current mass (kg)
# set_calc!(mass_calc, 1500.0) updates mass

# Drag coefficient (configuration parameter)
cd_calc = SpacecraftCalc(sat, DragCoefficient())
# get_calc(cd_calc) returns Cd value
# set_calc!(cd_calc, 2.2) updates coefficient

# Pattern: Same subject-variable binding, works with existing design
# DEFERRED: Implementation not required for initial release
```

---

**Use Case (DEFERRED): ForceModel Properties**

Design validation for force model configuration parameters.

```julia
# Gravity field degree/order
degree_calc = ForceModelCalc(gravity_model, GravityDegree())
# get_calc(degree_calc) returns current degree
# set_calc!(degree_calc, 8) updates degree

# Integrator tolerance
tol_calc = IntegratorCalc(integrator, RelativeTolerance())
# get_calc(tol_calc) returns current reltol
# set_calc!(tol_calc, 1e-10) updates tolerance

# Pattern: Same subject-variable binding, extends naturally to new domains
# DEFERRED: Implementation not required for initial release
```

---

**Use Case: Relative/Elapsed Time**

Elapsed time from anchor captured at construction.

```julia
# Time interval between two epochs (explicit endpoints)
t1 = Time("2025-01-01T00:00:00", UTC(), ISOT())
t2 = Time("2025-01-15T12:00:00", UTC(), ISOT())
dt_calc = TimeInterval(t1, t2, Days())
# get_calc(dt_calc) returns 14.5 (days)

# Elapsed time from anchor (captured at construction)
# Use case: "propagate for 3 days", "stop after 2 hours", time-series plotting
elapsed_calc = RelativeTimeCalc(sat, ElapsedDays())
# Constructor deep copies sat.time as immutable anchor
# get_calc(elapsed_calc) returns days elapsed from anchor to current sat.time

# Different time units via variable type
elapsed_sec = RelativeTimeCalc(sat, ElapsedSeconds())
elapsed_hr = RelativeTimeCalc(sat, ElapsedHours())

# Use in propagation stopping conditions
propagate!(prop, sat, StopAt(sat, elapsed_calc, 3.0))  # Stop after 3 days from anchor

# Use in optimization (anchor remains fixed across iterations)
var_time = SolverVariable(calc=elapsed_calc, lower=0.0, upper=10.0)

# Use in plotting (anchor set once, tracks elapsed time)
times = []
smas = []
for step in propagation_history
    push!(times, get_calc(elapsed_calc))  # Elapsed days from anchor
    push!(smas, get_calc(OrbitCalc(sat, SMA())))
end
plot(times, smas)  # X-axis: days from mission start

# Pattern: RelativeTimeCalc follows subject-variable binding
# Variable types: ElapsedDays(), ElapsedSeconds(), ElapsedHours()
# Anchor captured at construction, remains immutable thereafter
```

---

**Use Case: Function-Based Calc**

Arbitrary computation using closure pattern (Pattern 3).

```julia
# Define custom calc with function + references
struct CustomCostCalc
    func::Function
    spacecraft::Spacecraft
    target_sma::Float64
end

# Implement get_calc (evaluates custom function)
function AstroCallbacks.get_calc(calc::CustomCostCalc)
    sma_calc = OrbitCalc(calc.spacecraft, SMA())
    current_sma = get_calc(sma_calc)
    return (current_sma - calc.target_sma)^2  # Quadratic cost
end

# Usage
cost_calc = CustomCostCalc(
    (sc, tgt) -> (get_calc(OrbitCalc(sc, SMA())) - tgt)^2,
    sat,
    42000.0  # target GEO radius
)
# get_calc(cost_calc) evaluates arbitrary function

# Pattern: Maximum flexibility, any computation, same interface
```

---

**Use Case: Attribute Queries**

External modules check calc attributes before use.

```julia
# Optimizer checks before creating SolverVariable
state_calc = OrbitCalc(sat, Cartesian())
if is_settable(state_calc) && is_numeric(state_calc)
    var = SolverVariable(calc=state_calc, lower=..., upper=...)
else
    error("Calc must be settable and numeric for optimization")
end

# Plotter queries dimensionality
sma_calc = OrbitCalc(sat, SMA())
if length(sma_calc) == 1
    plot_scalar_timeseries(sma_calc)
else
    plot_vector_timeseries(sma_calc)
end

# Reporter handles non-numeric calcs
time_calc = TimeCalc(sat.time, UTC(), ISOT())
if is_numeric(time_calc)
    println("Time: ", get_calc(time_calc))
else
    println("Time: ", get_calc(time_calc))  # String output
end
```

---

**Use Case: Error Cases**

External modules check calc attributes before use.

```julia
# Optimizer checks before creating SolverVariable
state_calc = OrbitCalc(sat, Cartesian())
if is_settable(state_calc) && is_numeric(state_calc)
    var = SolverVariable(calc=state_calc, lower=..., upper=...)
else
    error("Calc must be settable and numeric for optimization")
end

# Plotter queries dimensionality
sma_calc = OrbitCalc(sat, SMA())
if length(sma_calc) == 1
    plot_scalar_timeseries(sma_calc)
else
    plot_vector_timeseries(sma_calc)
end

# Reporter handles non-numeric calcs
time_calc = TimeCalc(sat.time, UTC(), ISOT())
if is_numeric(time_calc)
    println("Time: ", get_calc(time_calc))
else
    println("Time: ", get_calc(time_calc))  # String output
end
```

---
**Error Cases**

Invalid operations produce clear error messages.

```julia
# Error: Setting read-only calc
mu_calc = BodyCalc(earth, GravParam())
# set_calc!(mu_calc, 400000.0)
# ERROR: Cannot set_calc! on read-only BodyCalc(earth, GravParam())

# Error: Type mismatch
state_calc = OrbitCalc(sat, Cartesian())
# set_calc!(state_calc, "invalid")
# ERROR: Expected Vector{Float64}, got String

# Error: Out-of-range value
sma_calc = OrbitCalc(sat, SMA())
# set_calc!(sma_calc, -7000.0)
# ERROR: SMA must be positive, got -7000.0

# Error: Unsupported variable type
# calc = OrbitCalc(sat, UnsupportedVariable())
# ERROR: No method get_calc(::OrbitCalc{Spacecraft, UnsupportedVariable})
```

---

### Design Decisions and Rationale

**Why get/set pattern instead of property access?**
- Enables dispatch on variable type
- Supports conversions (e.g., time scale/format conversions)
- Allows validation and side effects
- Works with AD frameworks

**Why subject-variable binding pattern?**
- Type-safe: compiler catches invalid combinations
- Extensible: users can add new combinations via multiple dispatch
- Self-documenting: `OrbitCalc(sat, SMA())` is clear
- Alternative considered: string-based lookup ("sma") - rejected for type safety

**Why separate calc types per subject?**
- Clear ownership: OrbitCalc operates on spacecraft, ManeuverCalc on maneuvers
- Enables subject-specific optimizations
- Better error messages
- Alternative: single universal Calc type - rejected for clarity and dispatch

**Why mutable domain objects?**
- In-place propagation requires mutable spacecraft
- Optimization requires modifying variables without recreating objects
- Alternative: immutable with copy-on-write - rejected for performance

**Calc Attributes System**

External modules query calc capabilities before use:

- **Settability**: Can optimizer modify this? (`is_settable(calc)`)
- **Value Type**: Numeric or non-numeric? (`value_type(calc)`)
- **Continuity**: Suitable for gradient-based optimization? (`is_continuous(calc)`)
- **Dimensionality**: Scalar or vector? (`length(get_calc(calc))`)
- **Differentiability**: Supports AD? (implementation-dependent)

**Why attributes instead of fixed calc categories?**
- Flexible: same calc can be numeric in one format, string in another (TimeCalc)
- Extensible: users can add custom attributes
- Self-describing: calcs declare their own capabilities
- Enables runtime checks in external modules

### Conventions & Constraints

**Naming:**
- Calc types named by subject: `OrbitCalc`, `ManeuverCalc`, `BodyCalc`, `TimeCalc`
- Variable types named by quantity: `SMA`, `Inc`, `Cartesian`, `DeltaVMag`
- Follow Julia style: CamelCase for types, snake_case for functions

**Numeric Types:**
- Default: Float64 for all numeric calculations
- AD compatibility: all numeric calcs must work with ForwardDiff.Dual
- Vector returns: `Vector{Float64}` (not SVector - simplicity over performance here)

**Mutability:**
- Domain objects (Spacecraft, Maneuver, Time) are mutable
- Calc objects themselves are typically immutable (just hold references)
- Variable type objects (SMA(), Cartesian()) are immutable singletons

**Error Handling:**
- Throw errors for invalid inputs (wrong types, out-of-range values)
- Throw errors for unsupported operations (set_calc! on read-only calc)
- Validation at calc construction when possible
- Clear error messages: "Cannot set_calc! on read-only calc OrbitCalc(sat, SMA())"

**Units:**
- Distances: km
- Velocities: km/s  
- Angles: radians (internal), may accept degrees in constructors with explicit type
- Time: Julian Date (Float64) for numeric formats
- Mass: kg
- Gravitational parameter: km³/s²

## 4. Package Interactions (Where it fits)

**Package Dependencies:**
- AstroStates - spacecraft state representations (CartesianState, etc.)
- AstroManeuvers - maneuver types (ImpulsiveManeuver, etc.)
- AstroUniverse - celestial bodies
- AstroEpochs - Time type and time scale/format handling

**Used by:**
- AstroSolve - optimization variables and constraints
- Plotting packages - time series of orbit elements, state components
- Reporting tools - formatted output of mission data
- Custom analysis scripts - user-defined calculations

**Example Interactions:**

**AstroSolve Integration:**
```julia
# AstroSolve checks calc attributes before creating SolverVariable
state_calc = OrbitCalc(sat, Cartesian())

# Hypothetical: AstroSolve validates before use
if is_settable(state_calc) && is_continuous(state_calc)
    var = SolverVariable(calc=state_calc, lower=..., upper=...)
else
    error("Calc must be settable and continuous for optimization")
end
```

**Plotting Integration:**
```julia
# Plotting package queries calc over time
sma_calc = OrbitCalc(sat, SMA())
sma_values = []
for t in time_range
    set_time!(sat, t)
    push!(sma_values, get_calc(sma_calc))
end
plot(time_range, sma_values)
```

**Reporting Integration:**
```julia
# Reporting accepts any calc, including non-numeric
time_calc = TimeCalc(sat.time, UTC(), ISOT())
sma_calc = OrbitCalc(sat, SMA())

println("Time: ", get_calc(time_calc))      # "2025-01-01T00:00:00"
println("SMA: ", get_calc(sma_calc), " km") # "7000.0 km"
```

## 5. Testing (How we verify it)

### Test categories:

**Input Validation**
- Invalid types rejected (e.g., passing non-spacecraft to OrbitCalc)
- Out-of-range values caught in set_calc! (e.g., negative SMA)
- Inconsistent combinations detected (e.g., unsupported variable types)
- Clear error messages provided

**Correctness**
- get_calc returns correct values from domain objects
- set_calc! properly updates underlying objects
- Mutability contracts honored (read-only calcs throw on set_calc!)
- Type stability maintained
- Vector calcs return correct dimensions

**Numeric Accuracy**
- State conversions accurate (Cartesian ↔ Keplerian)
- Time scale conversions accurate (UTC ↔ TDB ↔ TT)
- Tolerance: 1e-12 for dimensionless quantities
- Reference: validated against published ephemerides, analytical solutions

**Integration**
- Works with AstroStates state representations
- Compatible with ForwardDiff for AD
- Calcs properly update when domain objects change
- Composite calcs correctly aggregate values

**Regression**
- Fixed bugs have test coverage with issue references
- Known edge cases covered (singularities in orbital elements, time scale boundaries, etc.)

### Notes:
- Calc-specific edge cases in individual calc tests
- Performance benchmarks not required unless performance issues identified

## 6. Implementation Plan

- Refactor function names for attributes and test
  - Change `calc_is_settable()` → `is_settable()`
  - Change `calc_numvars()` → `length()`
  - Add `is_numeric()` (new function)
  - Update all tests to use new naming convention
- Implement Cartesian() and other state reps
- **TODO: Resolve ManeuverCalc constructor interface inconsistency**
  - Current: Struct has 3 fields (man, sc, var) but examples show both 2-arg and 3-arg constructors
  - Examples in maneuvercalc_deltavvector.jl show: `ManeuverCalc(dv, DeltaVVector())`
  - Docstring in AstroCallbacks.jl shows: `ManeuverCalc(m, sc, DeltaVVector())`
  - Decision needed: Should constructors take minimum args? Does DeltaVVector need spacecraft?
  - Impact: Update all examples and docs to match chosen interface
  - Principle: Constructors should take minimum required arg set