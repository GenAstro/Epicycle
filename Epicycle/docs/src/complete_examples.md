# Complete Examples

**Integrated mission simulations and workflows**

## Propagation Basics

```julia
using Epicycle

# Spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    #name="SC-StopAt",
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)

# Forces + integrator
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))
println(get_state(sat, Keplerian()))

# Propagate to apoapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
println(get_state(sat, Keplerian()))

# Stop when |r| reaches 7000 km 
propagate(prop, sat, StopAt(sat, PosMag(), 7000.0))
println(get_state(sat, SphericalRADEC()))       

# Propagate to x-position crossing (increasing)
sol = propagate(prop, sat, StopAt(sat, PosX(), 7.5; direction=+1))
println(get_state(sat, Cartesian()))

# Propagate multiple spacecraft with multiple stopping conditions
sc1 = Spacecraft(); sc2 = Spacecraft() 
stop_sc1_node = StopAt(sc1, PosZ(), 0.0)
stop_sc2_periapsis = StopAt(sc2, PosDotVel(), 0.0; direction=+1)
propagate(prop, [sc1,sc2], stop_sc1_node, stop_sc2_periapsis)
```

## Impulsive Maneuvers

```julia
using Epicycle

# Create a default spacecraft
sat1 = Spacecraft()

# Create an impulsive maneuver in the Inertial frame
deltav2 = ImpulsiveManeuver(
      axes = Inertial(),
      g0 = 9.80665,
      Isp = 250.0,
      element1 = 0.04,
      element2 = -0.3,
      element3 = 0.1
     )

# Apply the maneuver to the spacecraft
println("Initial mass: ", sat1.mass)
maneuver(sat1, deltav2)
println("Mass after Inertial maneuver: ", sat1.mass)
println("State after Inertial maneuver: \n", get_state(sat1, Cartesian()))

# Create an impulsive maneuver in the VNB frame
deltav1 = ImpulsiveManeuver(
      axes = VNB(),
      g0 = 9.80665,
      Isp = 250.0,
      element1 = 0.2,
      element2 = 0.1,
      element3 = -0.2
      )

# Apply the maneuver to the spacecraft
maneuver(sat1, deltav1)
println("Mass after VNB maneuver: ", sat1.mass)
println("State after VNB maneuver: \n", get_state(sat1, Cartesian()))
```

## Hohmann Transfer

```julia


using Epicycle

# Create spacecraft
sat = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
            )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat]
          )

# Create maneuver models for the hohmann transfer
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.2,
    element3 = 0.3
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.4,
    element2 = 0.5,
    element3 = 0.6
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat, DeltaVVector()),
    name = "toi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0],
)

# Define moi as a solver variable
var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat, DeltaVVector()),
    name = "moi",
    lower_bound = [-10.0, 0.0, 0.0],
    upper_bound = [10.0, 0.0, 0.0]
)

pos_target = 45000.0
pos_con = Constraint(
    calc = OrbitCalc(sat, PosMag()),
    lower_bounds = [pos_target],
    upper_bounds = [pos_target],
    scale = [1.0],
)

ecc_con = Constraint(
    calc = OrbitCalc(sat, Ecc()),
    lower_bounds = [0.0],
    upper_bounds = [0.0], 
    scale = [1.0],
)

# Create the TOI Event
toi_fun() = maneuver(sat, toi) 
toi_event = Event(name = "TOI", 
                  event = toi_fun, 
                  vars = [var_toi],
                  funcs = [])

# Create the prop to apopasis event

prop_apo_fun() = propagate(dynsys, integ, StopAtApoapsis(sat))
prop_event = Event(name = "Prop to Apoapsis", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver(sat, moi)
moi_event = Event(name = "MOI", 
                  event = moi_fun,
                  vars = [var_moi],
                  funcs = [pos_con, ecc_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

# Solve the trajectory and report results
result = trajectory_solve(seq)
sequence_report(seq)
solution_report(seq, result)
```

The output of this run:

```output
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:        4
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:        0

Total number of variables............................:        2
                     variables with only lower bounds:        0
                variables with lower and upper bounds:        2
                     variables with only upper bounds:        0
Total number of equality constraints.................:        2
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  0.0000000e+00 3.74e+04 0.00e+00   0.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  0.0000000e+00 1.77e+04 1.40e+02  -5.0 1.12e+01    -  4.69e-01 2.18e-01h  3
   2  0.0000000e+00 3.67e+03 1.25e-01  -0.6 1.55e-01    -  1.00e+00 1.00e+00h  1
   3  0.0000000e+00 2.35e+02 2.62e-02  -2.6 5.60e-02    -  1.00e+00 1.00e+00h  1
   4  0.0000000e+00 1.09e+00 3.48e-04  -4.5 4.41e-03    -  9.99e-01 1.00e+00h  1
   5  0.0000000e+00 2.40e-05 0.00e+00  -6.3 2.21e-05    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 5

                                   (scaled)                 (unscaled)
Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
Dual infeasibility......:   0.0000000000000000e+00    0.0000000000000000e+00
Constraint violation....:   7.2005827011268791e-07    2.3959924874361604e-05
Variable bound violation:   0.0000000000000000e+00    0.0000000000000000e+00
Complementarity.........:   4.7560577722374565e-07    4.7560577722374565e-07
Overall NLP error.......:   7.2005827011268791e-07    2.3959924874361604e-05


Number of objective function evaluations             = 7
Number of objective gradient evaluations             = 6
Number of equality constraint evaluations            = 10
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 6
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 0
Total seconds in IPOPT                               = 3.144

EXIT: Optimal Solution Found.

TRAJECTORY SEQUENCE SUMMARY
==================================================

Sequence Overview:
- Total Events: 3
- Variable Objects: 2 (6 optimization variables)
- Constraint Objects: 2 (2 constraint functions)
- Execution Order: ["TOI" → "Prop to Apoapsis" → "MOI"]

EVENT DETAILS:
--------------------

Event 1: "TOI"
├─ Variable Objects (1): 3 optimization variables
│  └─ toi: DeltaVVector() (ManeuverCalc) (3 components)
│     ├─  Component 1: ∈ [-10.0, 10.0]
│     ├─  Component 2: = 0.0
│     └─  Component 3: = 0.0
└─ Constraint Objects (0): None

Event 2: "Prop to Apoapsis"
├─ Variable Objects (0): None
└─ Constraint Objects (0): None

Event 3: "MOI"
├─ Variable Objects (1): 3 optimization variables
│  ├─ moi: DeltaVVector() (ManeuverCalc) (3 components)
│  │  ├─  Component 1: ∈ [-10.0, 10.0]
│  │  ├─  Component 2: = 0.0
│  │  └─  Component 3: = 0.0
└─ Constraint Objects (2): 2 constraint functions
   ├─ PosMag() (OrbitCalc) = 45000.0
   └─ Ecc() (OrbitCalc) = 0.0

STATEFUL OBJECTS:
--------------------
- ImpulsiveManeuver (×2)
- Spacecraft

==================================================

TRAJECTORY SOLUTION REPORT
==================================================

OPTIMIZATION STATUS:
- Converged: Solve_Succeeded
- Variable Objects: 2 (6 optimization variables)
- Constraint Objects: 2 (2 constraint functions)

OPTIMIZATION VARIABLES:
-------------------------
toi (DeltaVVector() (ManeuverCalc)):
  Component 1: 2.355732 (bounds: [-10.0, 10.0])
  Component 2: 0.0 (fixed at 0.0)
  Component 3: 0.0 (fixed at 0.0)
  Total ΔV: 2.355732

moi (DeltaVVector() (ManeuverCalc)):
  Component 1: 1.432696 (bounds: [-10.0, 10.0])
  Component 2: 0.0 (fixed at 0.0)
  Component 3: 0.0 (fixed at 0.0)
  Total ΔV: 1.432696

CONSTRAINT SATISFACTION:
-------------------------
Event "MOI":
  PosMag() (OrbitCalc): 45000.000024 (target: 45000.0)
  Ecc() (OrbitCalc): 0.0 (target: 0.0)

==================================================
```
## GEO Transfer

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

The output for this configuration is

```output
This is Ipopt version 3.14.19, running with linear solver MUMPS 5.8.1.

Number of nonzeros in equality constraint Jacobian...:       28
Number of nonzeros in inequality constraint Jacobian.:        0
Number of nonzeros in Lagrangian Hessian.............:        0

Total number of variables............................:        7
                     variables with only lower bounds:        0
                variables with lower and upper bounds:        7
                     variables with only upper bounds:        0
Total number of equality constraints.................:        4
Total number of inequality constraints...............:        0
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        0

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  0.0000000e+00 6.86e+04 0.00e+00   0.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  0.0000000e+00 5.34e+04 3.38e+00  -1.1 5.54e+00    -  4.81e-01 1.28e-01h  3
   2  0.0000000e+00 3.23e+04 1.09e+01  -0.8 1.46e+00    -  9.75e-01 2.50e-01h  3
   3  0.0000000e+00 4.01e+04 1.04e+00  -1.1 3.53e-01    -  7.31e-01 1.00e+00H  1
   4  0.0000000e+00 7.86e+03 2.52e+00  -1.0 9.37e-01    -  1.00e+00 5.00e-01h  2
   5  0.0000000e+00 5.46e+02 1.61e-01  -2.0 1.41e-01    -  1.00e+00 1.00e+00h  1
   6  0.0000000e+00 4.68e+00 2.37e-03  -3.7 2.42e-03    -  1.00e+00 1.00e+00h  1
   7  0.0000000e+00 3.82e-04 3.54e-05  -5.6 3.54e-05    -  1.00e+00 1.00e+00h  1
   8  0.0000000e+00 1.18e-07 3.02e-12 -11.0 3.92e-09    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 8

                                   (scaled)                 (unscaled)
Objective...............:   0.0000000000000000e+00    0.0000000000000000e+00
Dual infeasibility......:   3.0156910335045697e-12    3.0156910335045697e-12
Constraint violation....:   1.5534043499097701e-09    1.1797237675637005e-07
Variable bound violation:   0.0000000000000000e+00    0.0000000000000000e+00
Complementarity.........:   8.2697560462276952e-11    8.2697560462276952e-11
Overall NLP error.......:   1.5534043499097701e-09    1.1797237675637005e-07


Number of objective function evaluations             = 18
Number of objective gradient evaluations             = 9
Number of equality constraint evaluations            = 22
Number of inequality constraint evaluations          = 0
Number of equality constraint Jacobian evaluations   = 9
Number of inequality constraint Jacobian evaluations = 0
Number of Lagrangian Hessian evaluations             = 0
Total seconds in IPOPT                               = 2.412

EXIT: Optimal Solution Found.

TRAJECTORY SEQUENCE SUMMARY
==================================================

Sequence Overview:
- Total Events: 8
- Variable Objects: 3 (9 optimization variables)
- Constraint Objects: 4 (4 constraint functions)
- Execution Order: ["Prop to Z 1" → "TOI" → ... → "MOI"]

EVENT DETAILS:
--------------------

Event 1: "Prop to Z 1"
├─ Variable Objects (0): None
└─ Constraint Objects (0): None

Event 2: "TOI"
├─ Variable Objects (1): 3 optimization variables
│  └─ toi_v: DeltaVVector() (ManeuverCalc) (3 components)
│     ├─  Component 1: ∈ [-5.0, 5.0]
│     ├─  Component 2: = 0.0
│     └─  Component 3: = 0.0
└─ Constraint Objects (0): None

Event 3: "Prop to Apoapsis"
├─ Variable Objects (0): None
└─ Constraint Objects (1): 1 constraint functions
   └─ PosMag() (OrbitCalc) = 85000.0

Event 4: "Prop to Perigee 1"
├─ Variable Objects (0): None
└─ Constraint Objects (0): None

Event 5: "Prop to Z 2"
├─ Variable Objects (0): None
└─ Constraint Objects (0): None

Event 6: "MCC"
├─ Variable Objects (1): 3 optimization variables
│  └─ mcc_vn: DeltaVVector() (ManeuverCalc) (3 components)
│     ├─  Component 1: ∈ [-2.0, 2.0]
│     ├─  Component 2: ∈ [-2.0, 2.0]
│     └─  Component 3: ∈ [-0.001, 0.001]
└─ Constraint Objects (0): None

Event 7: "Prop to Perigee 2"
├─ Variable Objects (0): None
└─ Constraint Objects (2): 2 constraint functions
   ├─ Inc() (OrbitCalc) = 0.034907
   └─ PosMag() (OrbitCalc) = 42195.0

Event 8: "MOI"
├─ Variable Objects (1): 3 optimization variables
│  ├─ moi_v: DeltaVVector() (ManeuverCalc) (3 components)
│  │  ├─  Component 1: ∈ [-2.0, 2.0]
│  │  ├─  Component 2: ∈ [-0.001, 0.001]
│  │  └─  Component 3: ∈ [-0.001, 0.001]
└─ Constraint Objects (1): 1 constraint functions
   └─ SMA() (OrbitCalc) = 42166.9

STATEFUL OBJECTS:
--------------------
- ImpulsiveManeuver (×3)
- Spacecraft

==================================================

TRAJECTORY SOLUTION REPORT
==================================================

OPTIMIZATION STATUS:
- Converged: Solve_Succeeded
- Variable Objects: 3 (9 optimization variables)
- Constraint Objects: 4 (4 constraint functions)

OPTIMIZATION VARIABLES:
-------------------------
toi_v (DeltaVVector() (ManeuverCalc)):
  Component 1: 2.819829 (bounds: [-5.0, 5.0])
  Component 2: 0.0 (fixed at 0.0)
  Component 3: 0.0 (fixed at 0.0)
  Total ΔV: 2.819829

mcc_vn (DeltaVVector() (ManeuverCalc)):
  Component 1: 0.699575 (bounds: [-2.0, 2.0])
  Component 2: 0.895256 (bounds: [-2.0, 2.0])
  Component 3: 0.0 (bounds: [-0.001, 0.001])
  Total ΔV: 1.136173

moi_v (DeltaVVector() (ManeuverCalc)):
  Component 1: -0.480757 (bounds: [-2.0, 2.0])
  Component 2: -0.0 (bounds: [-0.001, 0.001])
  Component 3: -0.0 (bounds: [-0.001, 0.001])
  Total ΔV: 0.480757

CONSTRAINT SATISFACTION:
-------------------------
Event "Prop to Apoapsis":
  PosMag() (OrbitCalc): 85000.0 (target: 85000.0)

Event "Prop to Perigee 2":
  Inc() (OrbitCalc): 0.034907 (target: 0.034907)
  PosMag() (OrbitCalc): 42195.0 (target: 42195.0)

Event "MOI":
  SMA() (OrbitCalc): 42166.9 (target: 42166.9)

==================================================
```