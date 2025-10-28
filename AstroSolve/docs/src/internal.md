# Internal Details

This page provides implementation details, advanced usage patterns, and internal architecture information for AstroSolve.

## Architecture Overview

### Solver Framework Design

AstroSolve implements a constraint-based optimization framework specifically designed for astrodynamics problems. The architecture separates concerns into three main layers:

1. **Variable Layer** - `SolverVariable` objects encapsulate solver-controlled parameters
2. **Constraint Layer** - Functions that define mission requirements and objectives  
3. **Sequence Layer** - `Event` objects that organize mission phases with dependencies

### Event-Driven Mission Modeling

The event system allows complex missions to be modeled as a directed acyclic graph (DAG) where:

- **Events** represent mission phases (burns, coast arcs, flybys)
- **Dependencies** define the temporal or logical ordering of events
- **Variables** can be shared between events for global optimization
- **Constraints** enforce mission requirements at each event

## Advanced Usage Patterns

### Custom Variable Types

You can extend the solver framework with custom variable types by implementing the `AbstractCalc` interface:

```julia
# Custom calc for spacecraft mass optimization
struct MassCalc <: AbstractCalc
    spacecraft::Spacecraft
    var::MassVariable
end

# Define how to get/set the variable
get_calc(calc::MassCalc) = calc.spacecraft.mass
set_calc!(calc::MassCalc, value) = (calc.spacecraft.mass = value)

# Register with solver
mass_var = SolverVariable(calc=MassCalc(spacecraft, MassVariable()),
                         lower_bound=100.0, upper_bound=5000.0)
```

### Multi-Phase Mission Design

Complex missions can be assembled from modular events:

```julia
# Earth departure phase
earth_departure = Event(
    name="Earth Departure",
    event=() -> begin
        apply_maneuver!(spacecraft, departure_burn)
        propagate!(spacecraft, earth_escape_time)
    end,
    vars=[departure_dv, departure_time],
    funcs=[c3_constraint, launch_window_constraint]
)

# Deep space cruise
cruise_phase = Event(
    name="Cruise Correction",
    event=() -> begin
        propagate!(spacecraft, correction_time)
        apply_maneuver!(spacecraft, correction_burn)
    end,
    vars=[correction_dv],
    funcs=[trajectory_error_constraint]
)

# Mars arrival
mars_arrival = Event(
    name="Mars Orbit Insertion",
    event=() -> begin
        propagate!(spacecraft, mars_arrival_time)
        apply_maneuver!(spacecraft, insertion_burn)
    end,
    vars=[insertion_dv, arrival_time],
    funcs=[periapsis_constraint, inclination_constraint]
)

# Build complete mission
sequence = Sequence()
add_events!(sequence, [earth_departure, cruise_phase, mars_arrival])
```

### Solver Configuration

AstroSolve provides several options for configuring the underlying optimization:

```julia
# For large problems, use IPOPT
using Ipopt
solution = trajectory_solve!(sequence, solver=:ipopt, 
                           max_iter=1000, tolerance=1e-8)

# For smaller problems, use NLsolve
solution = trajectory_solve!(sequence, solver=:nlsolve,
                           method=:newton, autodiff=:forward)
```

## Performance Considerations

### Memory Management

AstroSolve uses stateful struct management to minimize memory allocations during iterative solving:

```julia
# Structs marked as stateful are automatically reset between iterations
is_astrosolve_stateful(::Type{Spacecraft}) = true
is_astrosolve_stateful(::Type{ImpulsiveManeuver}) = true

# Custom types can opt into this system
is_astrosolve_stateful(::Type{MyCustomType}) = true
```

### Variable Scaling

Proper variable scaling is crucial for convergence in astrodynamics problems:

```julia
# Position variables (km)
pos_var = SolverVariable(calc=pos_calc, 
                        scale=6378.0, shift=0.0)  # Earth radius scaling

# Velocity variables (km/s)  
vel_var = SolverVariable(calc=vel_calc,
                        scale=7.8, shift=0.0)     # Circular velocity scaling

# Time variables (days)
time_var = SolverVariable(calc=time_calc,
                         scale=86400.0, shift=0.0) # Seconds to days
```

### Constraint Formulation

Well-conditioned constraints improve solver robustness:

```julia
# Instead of: altitude = 400 km exactly
# Use: 399 km ≤ altitude ≤ 401 km (tolerance band)

# Instead of: arrival_time = target_time exactly  
# Use: |arrival_time - target_time| ≤ tolerance

# Normalize constraint magnitudes to similar scales
periapsis_constraint = (r_p - target_r_p) / target_r_p  # Relative error
```

## Integration Points

### AstroFun Integration

AstroSolve leverages the AstroFun calculation framework:

```julia
# OrbitCalc variables for orbital elements
sma_var = SolverVariable(calc=OrbitCalc(spacecraft, SMA()))
ecc_var = SolverVariable(calc=OrbitCalc(spacecraft, ECC()))

# ManeuverCalc variables for ΔV optimization
dv_var = SolverVariable(calc=ManeuverCalc(maneuver, DeltaVMag()))
dv_vec_var = SolverVariable(calc=ManeuverCalc(maneuver, DeltaVVector()))

# BodyCalc variables for gravitational parameters
mu_var = SolverVariable(calc=BodyCalc(body, GravParam()))
```

### AstroProp Integration

Trajectory propagation within mission events:

```julia
earth_departure = Event(
    event=() -> begin
        # Apply departure maneuver
        apply_maneuver!(spacecraft, departure_burn)
        
        # Propagate to sphere of influence
        prop_result = propagate!(spacecraft, soi_time,
                                forces=[point_mass_gravity])
        
        # Update spacecraft state
        set_posvel!(spacecraft, prop_result.final_state)
    end,
    vars=[departure_dv, departure_time],
    funcs=[escape_energy_constraint]
)
```

## Debugging and Diagnostics

### Convergence Issues

Common convergence problems and solutions:

```julia
# Check variable bounds
vars = manager.ordered_vars
for var in vars
    println("$(var.name): [$(var.lower_bound), $(var.upper_bound)]")
end

# Monitor constraint values during iteration
function debug_callback(x)
    constraint_vals = get_fun_values(manager, x)
    println("Constraint residuals: ", constraint_vals)
    return false  # Continue iteration
end

solution = trajectory_solve!(sequence, callback=debug_callback)
```

### State Inspection

Examine mission state at each event:

```julia
# Add debugging to events
debug_event = Event(
    name="Debug Checkpoint",
    event=() -> begin
        println("Spacecraft position: ", spacecraft.position)
        println("Spacecraft velocity: ", spacecraft.velocity)
        println("Spacecraft mass: ", spacecraft.mass)
    end,
    vars=[], funcs=[]
)
```

## Implementation Notes

### Topological Sorting

The event dependency system uses Kahn's algorithm for topological sorting, ensuring that mission events execute in the correct order while detecting circular dependencies.

### Variable Ordering

Variables are ordered deterministically to ensure consistent optimization behavior across runs. The ordering considers:

1. Event order in the topological sort
2. Variable declaration order within events
3. Variable type and name for tie-breaking

### Constraint Evaluation

Constraints are evaluated in sequence order, with automatic caching to avoid redundant calculations during gradient computation.

---

For additional implementation details, see the source code in the [AstroSolve repository](https://github.com/GenAstro/Epicycle/tree/main/AstroSolve).