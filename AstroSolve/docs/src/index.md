# AstroSolve.jl

**Advanced trajectory optimization and mission design for aerospace applications**

AstroSolve provides a comprehensive framework for solving complex astrodynamics problems through constraint-based optimization and event-driven mission sequencing. Whether you're designing interplanetary transfers, station-keeping maneuvers, or multi-phase missions, AstroSolve offers the tools to model, solve, and optimize your trajectory problems.

## Key Features

### 🎯 **Trajectory Optimization**
- **Constraint-based solving** with automatic variable scaling and bounds management
- **Multi-objective optimization** supporting position, velocity, and maneuver constraints
- **Robust convergence** using industry-standard NLsolve algorithms
- **Flexible variable types** including orbital elements, spacecraft states, and maneuver parameters

### 🚀 **Mission Sequencing**  
- **Event-driven architecture** for modeling complex multi-phase missions
- **Automatic dependency resolution** with topological sorting of mission events
- **State management** for propagation continuity across mission phases
- **Modular design** allowing easy assembly of trajectory segments

### ⚙️ **Solver Integration**
- **Professional-grade solvers** with IPOPT integration for large-scale problems
- **Smart initialization** with automatic variable scaling and shift management
- **Convergence diagnostics** and solver configuration options
- **Memory-efficient** stateful struct management for iterative solving

## Quick Start

```julia
using AstroSolve, AstroFun, AstroModels, AstroStates

# Define mission variables (orbital elements, maneuver magnitudes, etc.)
sma_var = SolverVariable(calc=OrbitCalc(spacecraft, SMA()), 
                        lower_bound=6700.0, upper_bound=42000.0)
dv_var = SolverVariable(calc=ManeuverCalc(maneuver, DeltaVMag()),
                       lower_bound=0.0, upper_bound=5.0)

# Create mission events with constraints
departure = Event(name="Earth Departure",
                 event=() -> apply_maneuver!(spacecraft, departure_burn),
                 vars=[dv_var],
                 funcs=[altitude_constraint])

arrival = Event(name="Mars Arrival", 
               event=() -> propagate!(spacecraft, mars_arrival_time),
               vars=[sma_var],
               funcs=[periapsis_constraint])

# Build and solve the trajectory sequence
sequence = Sequence()
add_events!(sequence, [departure, arrival])
solution = trajectory_solve!(sequence)
```

## Mission Design Workflow

AstroSolve follows a natural mission design progression:

1. **🎯 Define Variables** - Set up solver-controlled parameters (orbital elements, maneuver components, timing)
2. **📐 Specify Constraints** - Define mission requirements (altitude limits, arrival conditions, fuel budgets)  
3. **🔗 Build Sequence** - Assemble mission events with dependencies and state transitions
4. **🚀 Solve & Optimize** - Let AstroSolve find feasible solutions meeting all constraints
5. **📊 Analyze Results** - Extract optimal trajectories and performance metrics

## Applications

### **Interplanetary Missions**
Design transfers between planets with gravity assists, deep space maneuvers, and arrival constraints. AstroSolve handles the complex optimization of launch windows, flyby sequences, and arrival conditions.

### **Satellite Operations**
Optimize station-keeping maneuvers, orbit transfers, and constellation deployment. Model fuel-optimal strategies while maintaining operational constraints.

### **Lunar Missions**
Plan Earth-Moon transfers, lunar orbit insertion, and surface approach trajectories. Handle the complex three-body dynamics and mission timeline constraints.

### **Mission Analysis**
Perform trade studies on mission architecture, evaluate alternative trajectory strategies, and assess sensitivity to design parameters.

## Integration with Epicycle

AstroSolve seamlessly integrates with other Epicycle packages:

- **AstroStates** - Orbital state representations and conversions
- **AstroProp** - High-fidelity trajectory propagation with perturbations
- **AstroFun** - Constraint functions and mission analysis calculations
- **AstroModels** - Spacecraft and celestial body modeling
- **AstroMan** - Maneuver planning and execution

## API Reference

```@index
```

```@autodocs
Modules = [AstroSolve]
Order = [:type, :function, :macro, :constant]
```

---

*AstroSolve: Where complex missions become solvable problems.*
