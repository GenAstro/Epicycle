# AstroSolve.jl

AstroSolve provides a comprehensive framework for solving complex astrodynamics problems through constraint-based optimization and event-driven mission sequencing. Whether you're designing interplanetary transfers, station-keeping maneuvers, or multi-phase missions, AstroSolve offers the tools to model, solve, and optimize your trajectory problems.

## Key Features



## Quick Start

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


## API Reference

```@index
```

```@autodocs
Modules = [AstroSolve]
Order = [:type, :function, :macro, :constant]
```

---

*AstroSolve: Where complex missions become solvable problems.*
