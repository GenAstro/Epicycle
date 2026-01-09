# Spacecraft History

Trajectory history is automatically recorded during spacecraft propagation and organized into segments. This guide covers the history data structures and how to access trajectory data. 

`SpacecraftHistory` stores the time-ordered position and velocity data generated during propagation. History is organized into **segments** - continuous chunks of trajectory separated by discontinuities (e.g., maneuvers) or mission phases.

## Data Structures

### SpacecraftHistory

Container for all trajectory segments:

```julia
mutable struct SpacecraftHistory
    segments::Vector{HistorySegment}      # Solution trajectory
    iterations::Vector{HistorySegment}    # Solver iterations (diagnostic)
    record_segments::Bool                  # Enable solution recording
    record_iterations::Bool                # Enable iteration recording
end
```

**Two types of segments:**
- **`segments`**: Final solution trajectory (mission data)
- **`iterations`**: Solver convergence iterations (diagnostic data)

### HistorySegment

A continuous segment of trajectory data:

```julia
struct HistorySegment
    times::Vector{Time{Float64}}              # Time points
    states::Vector{CartesianState{Float64}}   # Position/velocity states
    coordinate_system::CoordinateSystem        # Reference frame
    name::String                               # Segment identifier
    metadata::Dict{String, Any}                # Extensible metadata
end
```

**Key points:**
- States always stored as `Float64` (automatic conversion from AD types)
- Parallel `times` and `states` vectors
- Coordinate system applies to entire segment
- Optional name for identification

## Accessing History Data

Spacecraft history is automatically populated during propagation and maneuvers. The examples below use a spacecraft `sat` with populated history created by this setup:

!!! note "Full Epicycle Ecosystem Required"
    The setup below requires AstroProp and AstroManeuvers. The data access examples below work with any populated `SpacecraftHistory`.

```julia
using Epicycle
using LinearAlgebra
    
# Create spacecraft
sat = Spacecraft()

# Setup propagator
gravity = PointMassGravity(earth, ())  # Only Earth gravity
forces  = ForceModel(gravity)
integ   = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
prop    = OrbitPropagator(forces, integ)

# Setup maneuver
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 2.518,
)

# Propagate and maneuver - this populates history
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))
maneuver(sat, toi)
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

# sat.history now contains 3 segments
```


### Basic Access

```julia
# Access segments
history = sat.history
segments = history.segments

# Number of segments
n = length(history)

# Check if empty
is_empty = isempty(history)

# Get specific segment
segment = history[1]              # First segment
segment = history.segments[2]     # Second segment
```

### Iteration

```julia
# Iterate over all segments
for segment in sat.history
    println("Segment: $(segment.name)")
    println("  Points: $(length(segment.times))")
end
```

### Segment Data

```julia
segment = history.segments[1]

# Access times and states
times = segment.times                    # Vector{Time{Float64}}
states = segment.states                  # Vector{CartesianState{Float64}}

# Get coordinate system
coord_sys = segment.coordinate_system

# Get name and metadata
name = segment.name
metadata = segment.metadata
```

### Extracting Position/Velocity

```julia
# Loop over states in a segment
for (t, state) in zip(segment.times, segment.states)
    x, y, z = state.position
    vx, vy, vz = state.velocity
    # Process data...
end
```

## Solution vs Iterations

### Solution Trajectory (segments)

The `segments` field contains the final mission trajectory:

```julia
# Final trajectory after propagation or optimization
for segment in sat.history.segments
    # This is the actual mission trajectory
end
```

This is what you use for:
- Visualization (plotting orbits)
- Ephemeris generation
- Mission analysis
- Reports

### Solver Iterations (iterations)

The `iterations` field contains diagnostic data from trajectory optimization:

```julia
# Check if iterations were recorded
if !isempty(sat.history.iterations)
    println("Solver iterations: ", length(sat.history.iterations))
    
    # Analyze convergence
    for iteration in sat.history.iterations
        # Each iteration is a full trajectory attempt
    end
end
```

Iterations are **opt-in** and used for:
- Debugging convergence problems
- Visualizing optimization process
- Performance analysis

See [AstroSolve documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/) for `trajectory_solve(...; record_iterations=true)`.

## Recording Control

History recording is controlled by flags:

```julia
history = sat.history

# Check recording status
history.record_segments      # true = recording solution
history.record_iterations    # true = recording iterations

# Flags are typically managed by solver
# Users rarely need to modify these directly
```

**Default behavior:**
- `record_segments = true`: Solution trajectory is recorded
- `record_iterations = false`: Iterations not recorded (saves memory)

## Common Patterns

### Check for Data

```julia
# Check if spacecraft has history
if !isempty(sat.history)
    println("History contains $(length(sat.history)) segments")
else
    println("No history recorded")
end
```

### Multi-Segment Missions

```julia
# After multi-phase mission
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))
maneuver(sat, toi)
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))

# Inspect phases
for (i, segment) in enumerate(sat.history)
    println("Segment $i: $(segment.name)")
    t_start = segment.times[1]
    t_end = segment.times[end]
    println("  Duration: $(t_end - t_start)")
end
```

### Coordinate System Awareness

```julia
# Segments may have different coordinate systems
for segment in sat.history
    origin = segment.coordinate_system.origin
    axes = segment.coordinate_system.axes
    println("Frame: $(origin.name) $(typeof(axes))")
end
```

## Integration with Other Modules

### Automatic Recording

History is populated automatically:

**AstroProp**: Creates new segment when `propagate()` is called
```julia
propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.5))
# sat.history now contains propagation data
```

**AstroManeuvers**: Creates new segment at maneuver application
```julia
maneuver(sat, toi)
# New segment created for post-maneuver trajectory
```

See respective module documentation for details on how they interact with history.

### Visualization

**View3D** (in Epicycle module) automatically renders all segments:
```julia
using Epicycle

view = View3D()
add_spacecraft!(view, sat)  # Renders entire history
display_view(view)
```

Each segment can be rendered with different colors for multi-phase missions.

## See Also

- [Spacecraft Guide](spacecraft.md) - Spacecraft construction and fields
- [AstroProp Documentation](https://genastro.github.io/Epicycle/AstroProp/dev/) - Orbit propagation
- [AstroSolve Documentation](https://genastro.github.io/Epicycle/AstroSolve/dev/) - Trajectory optimization
- [Reference](reference.md) - Complete API documentation
