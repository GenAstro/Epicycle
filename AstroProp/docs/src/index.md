```@meta
CurrentModule = AstroProp
```

# AstroProp

AstroProp provides force models, orbital propagators, and stopping conditions
for modelling spacecraft motion. AstroProp provides interfaces to the extensive 
numerical integration libraries in Julia's OrdinaryDiffEq.jl.  AstroProp is tested against the General Mission Analysis Tool (GMAT).

## Quick Start

The example below shows how to propagate a spacecraft using various stopping conditions:

```julia
using AstroEpochs, AstroStates, AstroFrames, AstroUniverse 
using AstroModels, AstroCallbacks, AstroProp, OrdinaryDiffEq

# Spacecraft — define time, state, and physical properties
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
    mass=1000.0,
    drag=SphericalDrag(c_d=2.2, drag_area=10.0),
    srp=SphericalSRP(c_r=1.8, srp_area=10.0),
)

# Propagator - define forces, integrator, and propagator
gravity = HarmonicGravity(earth; degree=4, order=0, model=Zonal())
drag    = AtmosphericDrag(earth; model=Exponential())
srp     = SolarRadiationPressure(earth; shadow=DualCone())
forces  = ForceModel(gravity, drag, srp)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate for 1 hour (3600 seconds)
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

# Propagate to an absolute time
target_time = Time("2015-09-22T12:00:00", TDB(), ISOT())
propagate!(prop, sat, StopAt(sat, target_time))

# Propagate to periapsis (r·v = 0, increasing crossing)
propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
println(get_state(sat, Keplerian()))

# Propagate backward for 2 hours using negative duration
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -7200.0); direction=:infer)

# Stop when |r| reaches 7000 km 
propagate!(prop, sat, StopAt(sat, PosMag(), 7000.0))
println(get_state(sat, SphericalRADEC()))       

# Propagate multiple spacecraft with multiple stopping conditions
sc1 = Spacecraft(mass=1000.0, drag=SphericalDrag(c_d=2.2, drag_area=10.0))
sc2 = Spacecraft(mass=1000.0, drag=SphericalDrag(c_d=2.2, drag_area=10.0))
stop_sc1_node = StopAt(sc1, PosZ(), 0.0)
stop_sc2_periapsis = StopAt(sc2, PosDotVel(), 0.0; direction=+1)
propagate!(prop, [sc1, sc2], stop_sc1_node, stop_sc2_periapsis)
```

## Function Syntax

### propagate

Propagates one or more spacecraft under specified forces to one or more stopping conditions.

**Syntax:**

```julia
sol = propagate!(propagator, spacecraft, stops...; direction=:forward, kwargs...)
```

**Parameters:**

- `propagator`: An `OrbitPropagator` containing the force model and integrator configuration
- `spacecraft`: A `Spacecraft` or `Vector{Spacecraft}` to propagate
- `stops...`: One or more `StopAt` stopping conditions (varargs)
- `direction`: (optional) `:forward` (default), `:backward`, or `:infer` - controls time integration direction
- `kwargs...`: Additional keyword arguments passed to the ODE solver

**Returns:**
- `sol`: ODE solution object from DifferentialEquations.jl

**Common usage patterns:**

```julia
# Single spacecraft, single stop
propagate!(prop, sat, stop)

# Single spacecraft, multiple stops
propagate!(prop, sat, stop1, stop2, stop3)

# Multiple spacecraft, multiple stops
propagate!(prop, [sat1, sat2], stop1, stop2)

# With direction keyword
propagate!(prop, sat, stop; direction=:backward)
```

See the following sections for detailed configuration:
- [Force Model Configuration](#Force-Model-Configuration) - Setting up gravitational forces
- [Integrator Selection](#Integrator-Selection) - Choosing integrators and tolerances
- [Stopping Conditions](#Stopping-Conditions) - State-based and time-based stop syntax

## Propagator Configuration

An `OrbitPropagator` combines a force model and integrator configuration to define how spacecraft motion is computed. This section covers the configuration of both components.

**Basic setup pattern:**

```julia
# 1. Define the gravitational forces
gravity = PointMassGravity(central_body, (perturbers...,))
forces = ForceModel(gravity)

# 2. Configure the numerical integrator
integ = IntegratorConfig(algorithm; dt=step, reltol=rtol, abstol=atol)

# 3. Create the propagator
prop = OrbitPropagator(forces, integ)
```

### Force Model Configuration

The force model defines the dynamics for propagation. You build it by composing forces in a `ForceModel` — gravity, atmospheric drag, solar radiation pressure, and any others the sections below describe. Each force is configured independently and summed to give the total acceleration.  

**Basic usage:**

```julia
# Earth-centered with perturbations from Moon and Sun
gravity = PointMassGravity(earth, (moon, sun))
forces = ForceModel(gravity)
```

**Components:**

- **Central body**: The primary gravitational body (e.g., `earth`, `mars`, `sun`)
- **Perturbing bodies**: Tuple of additional bodies whose gravity affects the trajectory (e.g., `(moon, sun)`)

**Common configurations:**

```julia
# LEO - Earth only (fast, low-fidelity)
gravity_leo = PointMassGravity(earth, ())

# LEO/MEO - Earth with Moon and Sun (standard accuracy)
gravity_standard = PointMassGravity(earth, (moon, sun))

# Interplanetary - Sun-centered with planetary perturbations
gravity_interplanetary = PointMassGravity(sun, (earth, mars, jupiter))
```
### Integrator Selection

AstroProp leverages Julia's DifferentialEquations.jl ecosystem, providing access to a wide range of high-performance numerical integrators. The choice of integrator and its parameters affects both the accuracy and speed of your propagation.

#### Common Integrators

**Recommended integrators for orbital mechanics:**

- **`Tsit5()`**: Tsitouras 5th order adaptive method. Good default choice for most applications with moderate accuracy requirements.
- **`Vern9()`**: Verner 9th order adaptive method. Higher accuracy for demanding applications like precision orbit determination.
- **`DP5()`**: Dormand-Prince 5th order method. Classic choice, similar performance to Tsit5.
- **`Vern7()`**: Verner 7th order method. Balance between Vern9 accuracy and Tsit5 speed.

For a complete list of available integrators, see the [DifferentialEquations.jl documentation](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/#Full-List-of-Methods).

#### Integrator Parameters

The `IntegratorConfig` accepts several key parameters that control integration behavior:

**`dt`** - Initial/suggested step size (seconds)
- Sets the initial time step for adaptive integrators
- Typical values: 10-60 seconds for LEO, 60-600 seconds for GEO, 86400.0 for interplanetary
- Smaller steps increase computation time but may improve accuracy near discontinuities

**`reltol`** - Relative error tolerance
- Controls accuracy relative to the magnitude of the state
- Typical values: `1e-9` to `1e-12` for high-precision work, `1e-6` to `1e-9` for general use
- Smaller values = higher accuracy but slower computation

**`abstol`** - Absolute error tolerance  
- Controls absolute error floor (important when state components are near zero)
- Typical values: `1e-9` to `1e-12` for position/velocity
- Should generally match or be slightly smaller than `reltol`

**Example configurations:**

```julia
# Fast propagation (lower accuracy)
integ_fast = IntegratorConfig(Tsit5(); dt=60.0, reltol=1e-6, abstol=1e-6)

# Standard propagation (good balance)
integ_standard = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)

# High-precision propagation
integ_precise = IntegratorConfig(Vern9(); dt=10.0, reltol=1e-12, abstol=1e-12)
```

!!! tip "Starting Point"
    If you're unsure, start with `Tsit5()` with `dt=10.0`, `reltol=1e-9`, and `abstol=1e-9`. Adjust based on your accuracy requirements and performance needs.

## Force Models

You build a force model by constructing individual forces and adding them to a `ForceModel`. The
`ForceModel` sums the forces to apply the total acceleration during numerical integration. The
sections below describe how to configure each force and which forces are available in the
open-source and Enterprise versions.

### Gravity

Gravity is provided by `PointMassGravity` and `HarmonicGravity`. `PointMassGravity` treats the
central body and any additional bodies — the Moon, the Sun, the planets — as point masses.
`HarmonicGravity` adds the central body's non-spherical gravity field, evaluated to the degree and
order you specify.

```julia
grav = PointMassGravity(earth, (moon, sun))                            # central body + third bodies
grav = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())  # zonal gravity field, J2–J5
```

You pick the gravity field with the `model` keyword. The open-source version includes `Zonal` — the
J2 through J5 zonal harmonics, which capture the dominant flattening of the Earth and cover most
low-Earth-orbit analysis.

To add the Sun and Moon alongside a spherical-harmonic Earth field, use `PointMassGravity` with
`include_center = false` so it contributes only those bodies — the Earth's gravity comes from
`HarmonicGravity`, and isn't counted twice:

```julia
forces = ForceModel(
    HarmonicGravity(earth; degree = 5, order = 0, model = Zonal()),
    PointMassGravity(earth, (moon, sun); include_center = false),   # Sun & Moon only
)
```

!!! note "Enterprise"
    The Enterprise version adds the full gravity fields, `EGM96` and `EGM2008`, evaluated to high
    degree and order for precision work. You select one the same way — just change `model`:

    ```julia
    using EpicycleEnterprise
    grav = HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())
    ```

### Atmospheric drag

`AtmosphericDrag` computes drag from the spacecraft's velocity relative to the rotating atmosphere.
It uses the drag coefficient and area you set on the spacecraft, and gets the local air density from
the atmosphere model you pick.

```julia
sc.drag = SphericalDrag(; c_d = 2.2, drag_area = 10.0)    # drag properties, on the spacecraft
drag    = AtmosphericDrag(earth; model = Exponential())    # atmosphere model, on the force
```

The open-source version includes the `Exponential` atmosphere — a smooth analytic density profile
that's fast and works well for early analysis.

!!! note "Enterprise"
    The Enterprise version adds `MSISE00` (NRLMSISE-00), the empirical atmosphere used for
    operational drag work. It responds to solar and geomagnetic activity, taken from
    `SpaceIndices` tables:

    ```julia
    using EpicycleEnterprise
    drag = AtmosphericDrag(mars; model = MSISE00())
    ```

### Solar radiation pressure

`SolarRadiationPressure` computes the push of sunlight on the spacecraft. It uses the reflectivity
and area you set on the spacecraft, and accounts for eclipses with the shadow model you pick:
`NoShadow`, a simple `Cylindrical` shadow, or `DualCone` — a realistic shadow with both umbra and
penumbra.

```julia
sc.srp = SphericalSRP(; c_r = 1.3, srp_area = 10.0)
srp    = SolarRadiationPressure(earth; shadow = DualCone())
```

Once you've built the forces you want, add them to a `ForceModel`:

```julia
forces = ForceModel(grav, drag, srp)
```

## Stopping Conditions

AstroProp supports two categories of stopping conditions: state-based and time-based.

### State-Based Stopping Conditions

State-based stops trigger when a calculated quantity (position, velocity, or derived value) crosses a target value. These use calculation variables from AstroCallbacks:

```julia
# Stop at periapsis (r·v = 0, velocity increasing)
StopAt(sat, PosDotVel(), 0.0; direction=+1)

# Stop at apoapsis (r·v = 0, velocity decreasing)  
StopAt(sat, PosDotVel(), 0.0; direction=-1)

# Stop when radius reaches 7000 km (any direction)
StopAt(sat, PosMag(), 7000.0; direction=0)

# Stop at ascending node (z = 0, increasing)
StopAt(sat, PosZ(), 0.0; direction=+1)
```

The `direction` parameter specifies which zero-crossing triggers the stop:
- `+1`: Trigger when value is increasing (positive derivative)
- `-1`: Trigger when value is decreasing (negative derivative)  
- `0`: Trigger on any crossing (default)

!!! note "Event Crossing Direction"
    The `direction` parameter on `StopAt` for state-based stops controls which side of the zero-crossing triggers the callback. This is different from the `direction` keyword on `propagate!()` which controls the time integration direction.

### Time-Based Stopping Conditions

Time-based stops allow propagation for a specified duration or until an absolute epoch:

**Elapsed Time Stops:**

```julia
# Propagate forward for 3600 seconds
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

# Propagate forward for 2.5 days
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 2.5))

# Propagate backward for 1 hour (negative duration)
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -3600.0); direction=:infer)
```

**Absolute Time Stops:**

```julia

# Propagate to a future epoch
sat = Spacecraft(
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
)
target = Time("2015-09-22T12:00:00", UTC(), ISOT())
propagate!(prop, sat, StopAt(sat, target))

# Propagate backward to a past epoch
past = Time("2015-09-20T12:00:00", UTC(), ISOT())
propagate!(prop, sat, StopAt(sat, past); direction=:infer)
```

!!! warning "Time-Based Direction Parameter"
    Time-based stopping conditions must use `direction=0` (the default). The event crossing direction concept does not apply to time-based stops. Use the `direction` keyword on `propagate!()` to control backward vs forward propagation.

## Direction Keywords

The `propagate!()` function accepts a `direction` keyword to control time integration:

- **`:forward`** (default): Integrate forward in time. This is the most common case.
- **`:backward`**: Integrate backward in time explicitly.
- **`:infer`**: Automatically infer direction from the time-based stop condition.

**When to use `:infer`:**

The `:infer` keyword is particularly useful in optimization and when the propagation direction may vary:

```julia
# Optimization variable (can be positive or negative)
duration = optimize_parameter  # Could be -1000.0 or +1000.0

# Using :infer allows the sign to determine direction automatically
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), duration); direction=:infer)
```

**Direction Sign Semantics:**

For `PropDurationSeconds` and `PropDurationDays`, the sign of the duration encodes the propagation direction:
- Positive duration: Forward propagation
- Negative duration: Backward propagation

This design enables clean optimization code where the duration variable can explore both positive and negative values.

!!! tip "Optimization Compatibility"
    When using time-based stops with optimization, use `direction=:infer` and let the duration sign indicate direction. This avoids the need for manual if-tests to switch between forward and backward propagation.

**Conflict Detection:**

AstroProp validates that explicit directions don't contradict duration signs:

```julia
# These cause errors:
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -100.0); direction=:forward)  # Error!
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 2.0); direction=:backward)      # Error!

# These are valid:
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -100.0); direction=:backward) # ✓
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), -100.0); direction=:infer)    # ✓
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 100.0); direction=:forward)   # ✓
```

## Time Scale Handling

AstroProp automatically selects the appropriate time scale for propagation based on the central body of the force model:

- **Earth-centered propagation**: Use **Terrestrial Time (TT)**
- **All other bodies**: Use **Barycentric Dynamical Time (TDB)**

This ensures the correct dynamical time scale is used in the integration of the equations of motion. Time scale conversions are handled automatically:
- Input times (on the spacecraft or in `StopAt`) can be in any scale
- The propagator converts to the appropriate scale (TT or TDB) based on the force model's central body
- The final spacecraft time is updated in the same propagation scale

!!! note "Force Model Central Body"
    The integration time scale is determined by the central body in your dynamics model. You can express spacecraft states in any coordinate or time system and AstroProp will still use the appropriate dynamical time scale for the integration of the equations of motion under the hood. 

## Core Functions

```@docs
OrbitPropagator
IntegratorConfig
propagate!
StopAt
```

## API Reference

The core API is documented in the sections above; this reference sweeps up the remaining public
symbols. The `Filter` excludes the symbols already given a dedicated `@docs` block (here and on the
[Force Models](force_models.md) page) so nothing is documented twice.

```@autodocs
Modules = [AstroProp]
Order = [:type, :function, :macro, :constant]
Public = true
Filter = t -> !(t in (
    PointMassGravity, HarmonicGravity, AbstractGeopotential, Zonal,
    AtmosphericDrag, AbstractDensityModel, Exponential, SolarRadiationPressure,
    OrbitPropagator, IntegratorConfig, propagate!, StopAt,
))
```
# Index

```@index
```