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

# Create spacecraft with default orbital state
sat = Spacecraft()

# Create the propagator with point mass gravity model
gravity = PointMassGravity(earth, (moon, sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Create an impulsive maneuver (only V component will vary)
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
)

# Define DeltaVVector of toi as a solver variable
# V component allowed to vary, N and B components constrained to zero
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define a constraint on position magnitude of spacecraft
# Target apoapsis altitude of 55000 km
pos_target = 55000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
)

# Create an event that applies the maneuver with toi as optimization variable
fun_toi() = maneuver(sat, toi) 
toi_event = Event(name = "TOI Maneuver", 
                  event = fun_toi,
                  vars = [var_toi],
                  funcs = [])

# Create propagation event to apoapsis with position constraint
fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_event = Event(name = "Propagate to Apoapsis", 
                   event = fun_prop_apo,
                   funcs = [pos_con])

# Build sequence: maneuver first, then propagate to apoapsis
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 

# Solve trajectory optimization using default settings (finite differences, IPOPT)
result = trajectory_solve(seq)

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
trajectory_solve
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

```julia

using Epicycle

# ========== Create Models =================================================================
# Create spacecraft
sat = Spacecraft(
    state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
    time = Time("2000-01-01T11:59:28.000", UTC(), ISOT())  # J2000 epoch from GMAT
)

# Create simple Earth point mass dynamics (no third bodies for now)
gravity = PointMassGravity(earth, ())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# Define maneuver models
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 1.518,
    element2 = 0.0,
    element3 = 0.0,
)

mcc = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.559,
    element2 = 0.588,
    element3 = 0.0,
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = -0.282,
    element2 = 0.0,
    element3 = 0.0,
)

# ========== Define Solver Variables ======================================================
# TOI: Only vary V component (Element1)
var_toi_v = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi_v",
    lower_bound = [-5.0, 0.0, 0.0],  # Only V component varies significantly
    upper_bound = [5.0, 0.0, 0.0],
)

# MCC: Vary V and N components (Element1 and Element2)
var_mcc_vn = SolverVariable(
    calc = ManeuverCalc(mcc, sat, DeltaVVector()),
    name = "mcc_vn",
    lower_bound = [-2.0, -2.0, -0.001],   # V and N components vary
    upper_bound = [2.0, 2.0, 0.001],
)

# MOI: Only vary V component (Element1)
var_moi_v = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi_v",
    lower_bound = [-2.0, -0.001, -0.001],  # Only V component varies significantly
    upper_bound = [2.0, 0.001, 0.001],
)

# ========== Define Constraints ===========================================================
# Constraint: Apoapsis radius = 85,000 km after TOI
apogee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [85000.0],
    upper_bounds = [85000.0],
    scale = [1.0],
)

# Constraint: Inclination = 2° after MCC
inclination_con = Constraint(
    calc = OrbitCalc(sat, Inc()),
    lower_bounds = [deg2rad(2.0)],
    upper_bounds = [deg2rad(2.0)],
    scale = [1.0],
)

# Constraint: Perigee radius = 42,195 km after MCC
perigee_radius_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [42195.0],
    upper_bounds = [42195.0],
    scale = [1.0],
)

# Constraint: Final SMA = 42,166.90 km (GEO) after MOI
final_sma_con = Constraint(
    calc = OrbitCalc(sat, SMA()),
    lower_bounds = [42166.90],
    upper_bounds = [42166.90],
    scale = [1.0],
)

# ========== Define events (variables, constraints, actions) ==============================
# Event 1: Propagate to Z=0 crossing (equatorial plane)
prop_to_z_crossing_1_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_1_event = Event(
    name = "Prop to Z 1",
    event = prop_to_z_crossing_1_fun,
)

# Event 2: Apply TOI maneuver
toi_fun() = maneuver(sat, toi)
toi_event = Event(
    name = "TOI",
    event = toi_fun,
    vars = [var_toi_v],
)

# Event 3: Propagate to apoapsis and check radius constraint
prop_to_apogee_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
prop_to_apogee_event = Event(
    name = "Prop to Apoapsis",
    event = prop_to_apogee_fun,
    funcs = [apogee_radius_con],
)

# Event 4: Propagate to perigee
prop_to_perigee_1_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_1_event = Event(
    name = "Prop to Perigee 1",
    event = prop_to_perigee_1_fun,
)

# Event 5: Propagate to Z=0 crossing again
prop_to_z_crossing_2_fun() = propagate(prop, sat, StopAt(sat, PosZ(), 0.0))
prop_to_z_crossing_2_event = Event(
    name = "Prop to Z 2",
    event = prop_to_z_crossing_2_fun,
)

# Event 6: Apply MCC maneuver
mcc_fun() = maneuver(sat, mcc)
mcc_event = Event(
    name = "MCC",
    event = mcc_fun,
    vars = [var_mcc_vn],
)

# Event 7: Propagate to perigee and check constraints
prop_to_perigee_2_fun() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=1))
prop_to_perigee_2_event = Event(
    name = "Prop to Perigee 2",
    event = prop_to_perigee_2_fun,
    funcs = [inclination_con, perigee_radius_con],
)

# Event 8: Apply MOI maneuver and check final SMA
moi_fun() = maneuver(sat, moi)
moi_event = Event(
    name = "MOI",
    event = moi_fun,
    vars = [var_moi_v],
    funcs = [final_sma_con],
)

# ========== Build the Sequence and Solve =================================================
seq = Sequence()
add_events!(seq, toi_event, [prop_to_z_crossing_1_event])
add_events!(seq, prop_to_apogee_event, [toi_event])
add_events!(seq, prop_to_perigee_1_event, [prop_to_apogee_event])
add_events!(seq, prop_to_z_crossing_2_event, [prop_to_perigee_1_event])
add_events!(seq, mcc_event, [prop_to_z_crossing_2_event])
add_events!(seq, prop_to_perigee_2_event, [mcc_event])
add_events!(seq, moi_event, [prop_to_perigee_2_event])

# Set up IPOPT options
ipopt_options = Dict(
    "max_iter" => 1000,
    "tol" => 1e-6,
    "output_file" => "ipopt_geo_transfer$(rand(UInt)).out",
    "file_print_level" => 5,
    "print_level" => 5,
)
snow_options = Options(derivatives=ForwardFD(), solver=IPOPT(ipopt_options))

result = trajectory_solve(seq, snow_options)
sequence_report(seq)
solution_report(seq, result)
```
