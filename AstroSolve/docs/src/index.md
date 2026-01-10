# AstroSolve.jl

AstroSolve provides a framework for solving astrodynamics design problems through constrained optimization and event-driven architecture. AstroSolve builds an event sequence as a Directed Acyclic Graph (DAG), which will soon be fully differentiable using AD. The DAG model for optimization is inspired by the architecture used in NASA's Copernicus software.  

**Reference:** Jacob Williams, Robert Falck, and Izaak Beekman. "Application of Modern Fortran to Spacecraft Trajectory Design and Optimization." *AIAA/AAS Space Flight Mechanics Meeting*, Kissimmee, FL, 2018. [Available online](https://ntrs.nasa.gov/api/citations/20180000413/downloads/20180000413.pdf)

AstroSolve's architecture is based on a few key user-facing components:

| Component | Description |
|:----------|:------------|
| **Models** | Spacecraft, maneuvers, propagators, etc. that model the physics of the system |
| **SolverVariables** | Optimization variables that the solver can adjust during trajectory design |
| **Constraints** | Equality and inequality constraints that must be satisfied in the solution |
| **Events** | Discrete mission phases that execute actions (propagate, maneuver, etc.) and apply constraints |
| **Sequence** | DAG structure that defines event dependencies and execution order |

## Quick Start

The example below uses AstroSolve's Event Sequence architecture to solve a simple orbit raising problem using IPOPT.

```julia
using Epicycle

# ============================================================================
# Shared Resources
# ============================================================================

# Create spacecraft with default orbital state
sat = Spacecraft()

# Create force models, integrator, and propagator
gravity = PointMassGravity(earth, (moon, sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# ============================================================================
# Event 1: TOI - Transfer Orbit Insertion
# ============================================================================

# Define TOI maneuver
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
)

# Define solver variable for TOI delta-V
toi_var = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define TOI event struct with event function, solver variables, and constraints
toi_fun() = maneuver!(sat, toi) 
toi_event = Event(
    name = "TOI Maneuver", 
    event = toi_fun,
    vars = [toi_var],
    funcs = []
)

# ============================================================================
# Event 2: Propagate to Apoapsis
# ============================================================================

# Define constraint on position magnitude at apoapsis
pos_target = 55000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
)

# Define propagation event to apoapsis with position constraint
prop_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event = Event(
    name = "Propagate to Apoapsis", 
    event = prop_fun,
    funcs = [pos_con]
)

# ============================================================================
# Trajectory Optimization
# ============================================================================

# Create sequence and add events
seq = Sequence()
add_sequence!(seq, toi_event, prop_event) 

# Solve trajectory optimization using default settings (finite differences, IPOPT)
result = solve_trajectory!(seq)

# Write a report documenting sequence and solution
sequence_report(seq)
solution_report(seq, result)
```

## Core Functions

This section contains reference material for the key user facing elements of AstroSolve. Note: Constraint is currently in the AstroCallbacks module but will be moved to AstroSolve in a future release. 

```@docs
SolverVariable
Event
Sequence
add_events!
solve_trajectory!
sequence_report
solution_report
```
## API Reference

```@index
```

```@autodocs
Modules = [AstroSolve]
Order = [:type, :function, :macro, :constant]
```

---

## A More Complex Example

This example demonstrates GEO orbit insertion including plane change, using 8 sequential events (3 maneuvers, 5 propagations) with constraints applied at multiple trajectory points.  The solution employes a bi-elliptic transfer. 

```julia

using Epicycle

# ============================================================================
# Shared Resources
# ============================================================================

# Create spacecraft
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    time = Time("2000-01-01T11:59:28.000", UTC(), ISOT()), 
    name = "GeoSat-1"
)

# Create force models, integrator, and propagator
gravity = PointMassGravity(earth, ())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# ============================================================================
# Event 1: Propagate to Equatorial Plane Crossing
# ============================================================================

# Define propagation event to equatorial plane crossing
prop_to_z_crossing_1_fun() = propagate!(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_1_event = Event(
    name = "Prop to Z 1",
    event = prop_to_z_crossing_1_fun,
)

# ============================================================================
# Event 2: TOI - Transfer Orbit Insertion
# ============================================================================

# Define TOI maneuver
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 2.518,
    element2 = 0.0,
    element3 = 0.0,
)

# Define solver variable for TOI delta-V
toi_var = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi_v",
    lower_bound = [0.0, 0.0, 0.0],
    upper_bound = [8.0, 0.0, 0.0],
)

# Define TOI event struct with event function, solver variables, and constraints
toi_fun() = maneuver!(sat, toi)
toi_event = Event(
    name = "TOI",
    event = toi_fun,
    vars = [toi_var],
)

# ============================================================================
# Event 3: Propagate to Apoapsis
# ============================================================================

# Define constraint on radius at apoapsis
apogee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [85000.0],
    upper_bounds = [85000.0],
    scale = [1.0],
)

# Define propagation event to apoapsis with radius constraint
prop_to_apogee_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_to_apogee_event = Event(
    name = "Prop to Apoapsis",
    event = prop_to_apogee_fun,
    funcs = [apogee_radius_con],
)

# ============================================================================
# Event 4: Propagate to Perigee
# ============================================================================

prop_to_perigee_1_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_1_event = Event(
    name = "Prop to Perigee 1",
    event = prop_to_perigee_1_fun,
)

# ============================================================================
# Event 5: Propagate to Equatorial Plane Crossing Again
# ============================================================================

prop_to_z_crossing_2_fun() = propagate!(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_2_event = Event(
    name = "Prop to Z 2",
    event = prop_to_z_crossing_2_fun,
)

# ============================================================================
# Event 6: MCC - Mid-Course Correction
# ============================================================================

# Define MCC maneuver
mcc = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.559,
    element2 = 0.588,
    element3 = 0.0,
)

# Define solver variable for MCC delta-V
mcc_var = SolverVariable(
    calc = ManeuverCalc(mcc, sat, DeltaVVector()),
    name = "mcc_vn",
    lower_bound = [-1.0, -1.0, -0.001],
    upper_bound = [4.0, 1.0, 0.001],
)

# Define MCC event struct with event function and solver variables
mcc_fun() = maneuver!(sat, mcc)
mcc_event = Event(
    name = "MCC",
    event = mcc_fun,
    vars = [mcc_var],
)

# ============================================================================
# Event 7: Propagate to Perigee and Check Constraints
# ============================================================================

# Define constraints on inclination and perigee radius
inclination_con = Constraint(
    calc = OrbitCalc(sat, Inc()),
    lower_bounds = [deg2rad(2.0)],
    upper_bounds = [deg2rad(2.0)],
    scale = [1.0],
)

perigee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [42195.0],
    upper_bounds = [42195.0],
    scale = [1.0],
)

# Define propagation event to perigee with inclination and perigee radius constraints
prop_to_perigee_2_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_2_event = Event(
    name = "Prop to Perigee 2",
    event = prop_to_perigee_2_fun,
    funcs = [inclination_con, perigee_radius_con],
)

# ============================================================================
# Event 8: MOI - GEO Orbit Insertion
# ============================================================================

# Define MOI maneuver
moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.282,
    element2 = 0.0,
    element3 = 0.0,
)

# Define solver variable for MOI delta-V
moi_var = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi_v",
    lower_bound = [-1.0, -0.001, -0.001],
    upper_bound = [4.0, 0.001, 0.001],
)

# Define constraint on semi-major axis at GEO
final_sma_con = Constraint(
    calc = OrbitCalc(sat, SMA()),
    lower_bounds = [42166.90],
    upper_bounds = [42166.90],
    scale = [1.0],
)

# Define MOI event struct with event function, solver variables, and constraints
moi_fun() = maneuver!(sat, moi)
moi_event = Event(
    name = "MOI",
    event = moi_fun,
    vars = [moi_var],
    funcs = [final_sma_con],
)

# ============================================================================
# Trajectory Optimization
# ============================================================================

seq = Sequence()
add_sequence!(seq, prop_to_z_crossing_1_event, toi_event, prop_to_apogee_event,
              prop_to_perigee_1_event, prop_to_z_crossing_2_event, mcc_event,
              prop_to_perigee_2_event, moi_event)

# Set up IPOPT options
ipopt_options = Dict(
    "max_iter" => 1000,
    "tol" => 1e-6,
    "output_file" => "ipopt_geo_transfer$(rand(UInt)).out",
    "file_print_level" => 5,
    "print_level" => 5,
)
snow_options = Options(derivatives=ForwardFD(), solver=IPOPT(ipopt_options))

# Solve trajectory optimization with iteration recording enabled
result = solve_trajectory!(seq, snow_options; record_iterations=true)
sequence_report(seq)
solution_report(seq, result)

# Propagate about one day to see final orbit
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.1)) 

# Visualize with iterations
view = View3D()
add_spacecraft!(view, sat; show_iterations=true)
display_view(view)

```
