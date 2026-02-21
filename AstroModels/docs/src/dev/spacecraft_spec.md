# Spacecraft Specification

## History

### 1. Overview & Motivation

The history field of Spacecraft stores propagated orbit state data for use in creating ephemeris files, 3D plots, reports etc. History segments are used to delineate discontinuities, and for separate mission events such as maneuver!(), and propagate!().

- History is stored as a vector of Time and matrix of cartesian state numbers but does not use the CartesianState struct. It does not scale well with new data like CoordinateSystem, segment name etc. Users access history data directly and the interface does not hide implementation so it is brittle.

- Implement a composed solution that provides an extendible interface for adding more history attributes like CoordinateSystem. Provide a better user interface that hides the implementation of the data/struct from the user. This means, accessing data for plotting, reporting etc. and a clean interface for integration with AstroCallbacks. 

### 2. Requirements

#### Data Requirements
- Store Time and CartesianState per data point within a segment
- Store CoordinateSystem per segment
- Support segment metadata (name, description, event type)
- Support multiple segments with discontinuities between them
- The history struct shall be extensible to handle potential fields and future metadata needs

#### Query/Access Requirements

- Access state at specific point
- Access state at specific time (with interpolation if needed)
- Access full segment by index or name
- Access data over time range (potentially spanning multiple segments)
- Export to arrays/matrices for plotting libraries
- Provide efficient iteration over time-state pairs
- Support querying segment boundaries and metadata

#### Non-functional Requirements
- Efficient append operations during propagation
- Memory efficient for long-duration propagations
- Type stable for performance-critical operations
- Minimal overhead for single-segment use cases

#### Interface Requirements
- Provide clean API for AstroCallbacks integration (adding history points during propagation)
- Support ephemeris file generation (access to time-ordered state data)
- Support 3D plotting (efficient bulk data access)
- Hide implementation details from user code (abstraction barrier)

### 3. Proposed Data Structure

#### Core Types

```julia
"""
    HistorySegment

Represents a continuous segment of spacecraft trajectory history.
"""
struct HistorySegment
    times::Vector{Time}                      # Ordered vector of time points
    states::Vector{CartesianState{Float64}}  # State vectors (always Float64)
    coordinate_system::CoordinateSystem      # Reference frame for this segment
    name::String                             # Optional segment identifier
    metadata::Dict{String, Any}              # Extensible metadata storage
end

"""
    SpacecraftHistory

Container for spacecraft trajectory history, organized into segments.
"""
struct SpacecraftHistory
    segments::Vector{HistorySegment}         # Solution trajectory segments
    iterations::Vector{HistorySegment}       # Solver iteration segments (diagnostic)
    record_segments::Bool                    # Enable solution segment recording
    record_iterations::Bool                  # Enable iteration segment recording
end
```

#### Design Rationale

**Storage Type Choice:**
- History stores Float64 exclusively for simplicity and type stability
- AD types (Dual numbers) are automatically cast to Float64 during recording
- History is for output/analysis (plots, ephemeris), not gradient computation
- Eliminates type parameter complexity while maintaining performance

**HistorySegment Fields:**
- `times` and `states` as parallel vectors for efficient storage and access
- `coordinate_system` at segment level (assumes single frame per segment)
- `name` for human-readable identification (empty string if unnamed)
- `metadata` Dict provides extensibility without modifying struct

**SpacecraftHistory:**
- Separates solution trajectory (`segments`) from diagnostic data (`iterations`)
- Recording flags control which segments are populated during propagation
- Default behavior: record solution only, skip iterations for efficiency
- Thin abstraction layer to hide implementation

#### Alternative Considerations

1. **Single vectors vs. Vector of tuples:** Parallel vectors chosen for performance and easier bulk access
2. **Mutable vs. Immutable:** Immutable structs chosen for safety; append operations create new segments or history
3. **Dict metadata vs. explicit fields:** Dict provides flexibility; can refactor to explicit fields if common patterns emerge
4. **Separate coordinate system per point:** Too granular; segment-level is sufficient for most use cases

### 4. API Design

#### Constructors

```julia (default: record solutions, not iterations)
# Create empty history
SpacecraftHistory()

# Create empty segment
HistorySegment(coord_system::CoordinateSystem; 
               name::String="", 
               metadata::Dict{String,Any}=Dict{String,Any}())

# Create segment from existing data
HistorySegment(times::Vector{Time}, 
               states::Vector{CartesianState{Float64}}, 
               coord_system::CoordinateSystem;
               name::String="",
               metadata::Dict{String,Any}=Dict{String,Any}())
```

#### Adding/Updating History

```julia
# Add a point to an existing segment (mutates segment's times/states vectors)
# Automatically converts Dual/BigFloat to Float64
push_state!(segment::HistorySegment, time::Time, state::CartesianState)

# Add a new segment to history (routes based on history.record_segments and history.record_iterations flags)
# - If record_iterations=true: adds to history.iterations
# - If record_segments=true and record_iterations=false: adds to history.segments
# - Otherwise: no-op (segment not recorded)
push_segment!(history::SpacecraftHistory, segment::HistorySegment)

# Start a new segment (convenience for propagation - creates and adds empty segment)
new_segment!(history::SpacecraftHistory, 
             coord_system::CoordinateSystem; 
             name::String="",
             metadata::Dict{String,Any}=Dict{String,Any}())
```

#### Query/Access Methods

```julia
# Direct field access (no getters needed)
segment.times                    # Vector{Time{Float64}}
segment.states                   # Vector{CartesianState{Float64}}
segment.coordinate_system        # CoordinateSystem
segment.name                     # String
segment.metadata                 # Dict{String,Any}

history.segments                 # Vector{HistorySegment} - solution trajectory
history.iterations               # Vector{HistorySegment} - solver iterations (if recorded)
history.record_segments          # Bool - is solution recording enabled?
history.record_iterations        # Bool - is iteration recording enabled?
history.record_iterations        # Bool - is iteration recording enabled?
history[i]                       # Get segment by index from solution trajectory
length(history)                  # Number of solution segments
isempty(history)                 # Check if solution trajectory empty
isempty(segment)                 # Check if segment empty

# Iteration (iterates over solution segments only)tory)                 # Check if empty
isempty(segment)                 # Check if segment empty

# Iteration
for segment in history
    # Process each segment
end
```

#### Export/Conversion

Export functions TBD when integrating with plotting/reporting systems

#### Integration Points

```julia
# For AstroCallbacks - append during propagation
# Callback adds (time, state) to current segment
function record_history!(history::SpacecraftHistory, time::Time, state::CartesianState)
    current_segment = last(history.segments)
    push_state!(current_segment, time, state)
end

# For propagate!() - start new segment
function propagate!(sc::Spacecraft, ...)
    new_segment!(sc.history, sc.coordinate_system, name="propagate_\$(timestamp)")
    # ... propagation loop calls record_history!
end

# For maneuver!() - mark discontinuity with new segment  
function maneuver!(sc::Spacecraft, ...)
    new_segment!(sc.history, sc.coordinate_system, name="maneuver_\$(timestamp)")
    # ... maneuver updates state
end
```

#### Design Notes

- **Immutability compromise:** While structs are immutable, the vectors they contain can be mutated for efficiency during propagation. `push_state!` modifies the underlying vector.
- **Naming convention:** `push_*!` for mutations, direct field access for data retrieval
- **Type inference:** Return types preserve `T` parameter for type stability
- **Segment access:** Direct indexing `history[i]` in Phase 1, name-based lookup in Phase 2
- **No getter functions:** Fields are public, use direct access (e.g., `segment.times` not `get_times(segment)`)

### 5. Usage Examples

This section explores how the history design integrates with other system components and enables common user workflows.

#### Plotting Integration

**Use Case:** User wants to visualize trajectory after propagation
```julia
# Simple case - plot entire trajectory
propagate!(sat, duration=3days)
view = View3D()
add_spacecraft!(view, sat)
display_view(view)
```
**Requirements from History:**
- Efficient iteration over segments and time-state pairs
- Extract position data (x, y, z) from CartesianState
- Support per-segment visualization (different colors for each segment)
- Query coordinate system for proper frame rendering

#### AstroCallbacks Integration

**Use Case:** History is recorded automatically during propagation
```julia
# History recording happens automatically
propagate!(sat, duration=1day)
# sat.history now contains trajectory data
```

**Requirements from History:**
- Fast append operations (`push_state!`) called frequently during integration
- Access to current segment for appending
- Automatic segment creation when propagate!() or maneuver!() is called

**Internal Implementation Note:**
History recording is handled internally by the propagation system via callbacks. Users don't need to configure history recording - it just works.

#### Ephemeris File Generation

**Use Case:** Export trajectory to SPK or other ephemeris format
```julia
# Create ephemeris configuration
ephem = EphemerisFile(
    filename = "mission.bsp",
    format = :spk,
    interpolation = :hermite
)

# Add spacecraft trajectory
add_spacecraft!(ephem, sat)

# Write the file
write_ephemeris(ephem)
```

**Requirements from History:**
- Time-ordered state data access per segment
- Coordinate system metadata for frame specifications
- Segment boundaries map naturally to SPK segments
- Support iteration for different interpolation schemes

**Requirements from History:**
- Time-ordered state data access per segment
- Coordinate system metadata for frame specifications
- Segment boundaries map naturally to SPK segments
- Support iteration for different interpolation schemes

#### 2D Plotting

**Use Case:** Create 2D plots of trajectory data
```julia
# Create plot configuration
plot2d = Plot2D(
    title = "Altitude Profile",
    xlabel = "Time (days)",
    ylabel = "Altitude (km)"
)

# Add quantities to plot
add_quantity!(plot2d, sat, :altitude)
add_quantity!(plot2d, sat, :velocity_magnitude)

# Display the plot
display_plot(plot2d)
```

**Requirements from History:**
- Efficient extraction of time series data
- Access to full state vectors for computing derived quantities
- Support for multiple spacecraft on same plot
- Handle segment boundaries (discontinuities in plots)

#### Report Generation

**Use Case:** Generate mission analysis reports
```julia
# Create report configuration
report = MissionReport(
    title = "Orbital Analysis",
    format = :pdf
)

# Add report sections
add_section!(report, :trajectory_summary, sat)
add_section!(report, :maneuver_performance, sat)
add_section!(report, :orbit_statistics, sat)

# Generate the report
generate_report(report, "mission_report.pdf")
```

**Requirements from History:**
- Query time ranges and statistics per segment
- Access segment metadata (names, event types)
- Extract state data for analysis functions
- Support tabular data export for report tables

#### Common User Workflows

**Workflow 1: Quick trajectory inspection (segments transparent)**
```julia
propagate!(sat, duration=2days)
plot(sat.history)  # User doesn't think about segments
```

**Workflow 2: Multi-phase mission (segments visible)**
```julia
propagate!(sat, duration=1day)
maneuver!(sat, dv=[0.1, 0, 0])
propagate!(sat, duration=1day)

# User wants to see each phase
for (i, segment) in enumerate(history.segments)
    println("Segment \$i: \$(segment.name), \$(length(segment.times)) points")
end
```

**Workflow 3: Custom analysis (full API access)**
```julia
# User has direct field access
segment1_times = history.segments[1].times
segment1_states = history.segments[1].states
# Custom processing...
```

#### Design Implications

- **Flattening is key:** Many operations (plot, export) should work across all segments by default
- **Segment awareness optional:** User shouldn't need to know about segments for simple cases
- **Performance matters for callbacks:** `push_state!` will be called thousands of times
- **Metadata enables rich features:** Names, coordinate systems, metadata enable better reports/plots
- **Standard interfaces:** Consider implementing standard Julia interfaces (Tables.jl, Plots.jl recipes)

### 6. Implementation Details

#### Storage Format
mutable struct SpacecraftHistory
    segments::Vector{HistorySegment}
    iterations::Vector{HistorySegment}
    record_segments::Bool
    record_iterations::Bool
end
```
- Mutable to allow flag modification by solver
- Separate vectors for solution vs iteration segments
- No type parameters - history always stores Float64
- Internal vectors
- Simple wrapper around vector of segments
- No type parameters - history always stores Float64
- Immutable struct, but internal vector can grow via `push!`

**HistorySegment struct:**
```julia
struct HistorySegment
    times::Vector{Time}
    states::Vector{CartesianState{Float64}}
    coordinate_system::CoordinateSystem
    name::String
    metadata::Dict{String, Any}
end
```
- Parallel vectors (`times`, `states`) for cache-friendly access
- Concrete Float64 type for full type stability
- Immutable struct protects invariants, but vectors can be mutated internally

**Memory layout:**
- Contiguous vectors enable SIMD and efficient bulk operations
- Vector{CartesianState{Float64}} stores actual data inline (not pointers)
- Metadata Dict only allocated when needed (empty dict is lightweight)

#### Indexing Strategy

**Direct field access (internal use only):**
```julia
# Fast path for performance-critical code
segment = getfield(history, :segments)[idx]
times_vec = getfield(segment, :times)
```

**Property access via getproperty (user interface):**
```julia
Base.getproperty(h::SpacecraftHistory, s::Symbol) = 
    if s === :num_segments
        length(getfield(h, :segments))
    elseif s === :segments
        getfield(h, :segments)
    else
        getfield(h, s)
    end

Base.getproperty(seg::HistorySegment, s::Symbol) = 
    if s === :num_points
        length(getfield(seg, :times))
    else
        getfield(seg, s)
    end
```

**Array-like indexing:**
```julia
# Access segment by index
Base.getindex(h::SpacecraftHistory, i::Int) = h.segments[i]
Base.getindex(h::SpacecraftHistory, r::UnitRange) = h.segments[r]

# Iteration support
Base.iterate(h::SpacecraftHistory, state=1) = 
    state > length(h.segments) ? nothing : (h.segments[state], state+1)

Base.length(h::SpacecraftHistory) = length(h.segments)
```

**Segment lookup by name (Phase 2):**
```julia
Base.getindex(h::SpacecraftHistory, name::String) = 
    findfirst(seg -> seg.name == name, h.segments)
```

#### Memory Management

**Preallocating segments:**
- During propagation, segments grow incrementally
- Consider sizehint! for expected trajectory length:
```julia
function new_segment!(history, coord_sys; expected_points=1000, kwargs...)
    segment = HistorySegment(coord_sys; kwargs...)
    sizehint!(segment.times, expected_points)
    sizehint!(segment.states, expected_points)
    push!(history.segments, segment)
end
```



**Append performance:**
```julia
function push_state!(segment::HistorySegment, time::Time, state::CartesianState)
    # Direct field access for performance
    times = getfield(segment, :times)
    states = getfield(segment, :states)
    # Convert to Float64 (handles Dual, BigFloat, etc.)
    state_f64 = to_float64(state)
    push!(times, time)
    push!(states, state_f64)
    return nothing
end
```
- Amortized O(1) append due to vector growth strategy
- Critical path during integration - no allocations except vector resize
- `to_float64` handles automatic conversion from AD types

**Type conversion helper:**
```julia
# Conversion from any CartesianState{T} to CartesianState{Float64}
to_float64(state::CartesianState{Float64}) = state  # No-op for Float64
to_float64(state::CartesianState) = CartesianState(Float64.(to_vector(state)))
```

**Memory footprint estimate:**
- Each Time: ~24 bytes
- Each CartesianState{Float64}: ~48 bytes (6 Float64s)
- Per point overhead: ~72 bytes
- 10,000 points ≈ 700 KB
- Segment metadata: negligible unless Dict is populated

**Cleanup considerations:**
- History is immutable from spacecraft perspective
- Clear history by creating new empty history
- Segments are garbage collected when history is replaced
- No manual memory management needed

#### Type Stability

**Concrete Float64 storage:**
```julia
# History always uses concrete Float64 type
history = SpacecraftHistory()
segment = HistorySegment(coord_sys)

# Automatic conversion from any numeric type
state_dual = CartesianState{Dual}([...])  # From AD propagation
push_state!(segment, time, state_dual)     # Converts to Float64 internally

# Result is always type-stable Float64
stored_state = segment.states[1]  # CartesianState{Float64}
```

**Type stability benefits:**
- No type parameters eliminates parametric complexity
- Vector{CartesianState{Float64}} is fully concrete
- AD overhead isolated to propagation, not storage
- metadata Dict is type-unstable by design (flexibility > performance here)
- Use getfield() in hot loops to avoid getproperty overhead

**P# Solver Iteration Recording (Optional Feature)

**Requirement:** Support recording solver iterations for trajectory optimization diagnostics.

**Use Case:** During trajectory optimization, the solver may evaluate thousands of candidate trajectories. By default, only the final converged solution is recorded in history. For debugging convergence issues or visualizing the optimization process, users can opt to record all iterations.

**Design:**
- Separate storage: `iterations` vector contains iteration segments, `segments` contains solution
- Two independent flags control recording:
  - `record_segments::Bool = true` - Controls solution segment recording
  - `record_iterations::Bool = false` - Controls iteration segment recording
- Routing logic in `push_segment!`:
  ```julia
  function push_segment!(history, segment)
      if history.record_iterations
          push!(history.iterations, segment)
      elseif history.record_segments
          push!(history.segments, segment)
      end
      # If both false: segment not recorded (no-op)
  end
  ```

**Solver Integration (AstroSolve):**
```julia
# In solve_trajectory!(seq, options; record_iterations=false)

# Default mode: record_iterations=false
sc.history.record_segments = false      # Don't record during optimization
sc.history.record_iterations = false
# ... optimizer runs, history not populated ...
sc.history.record_segments = true       # Record final solution only
# ... execute sequence one final time ...

# Diagnostic mode: record_iterations=true  
sc.history.record_segments = false      # Don't record in segments
sc.history.record_iterations = true     # Do record in iterations
# ... optimizer runs, all iterations stored in history.iterations ...
sc.history.record_iterations = false    # Switch to solution mode
sc.history.record_segments = true
# ... execute sequence final time, solution in history.segments ...
```

**Benefits:**
- Clean separation: solution data vs diagnostic data
- Memory efficient default: iterations not recorded unless requested
- Opt-in diagnostics: power users can visualize convergence behavior
- Natural access: `history.segments` always contains solution trajectory

**Rationale:**
- Alternative considered: Single `segments` vector with `is_iteration::Bool` field per segment
  - Rejected: Mixing solution and diagnostic data complicates filtering for ephemeris, plots
  - Rejected: Requires checking metadata/field for every operation
- Separate vectors make intent clear: iterations are purely diagnostic, segments are mission data

**Example Usage:**
```julia
# Debug mode - see all optimization iterations
result = solve_trajectory!(seq, options; record_iterations=true)
view = View3D()
add_spacecraft!(view, sat)
# Visualize shows all iteration trajectories overlaid (diagnostic view)

# Production mode - clean final solution only (default)
result = solve_trajectory!(seq, options)  # record_iterations=false
view = View3D()  
add_spacecraft!(view, sat)
# Visualize shows only final converged trajectory
```

###erformance validation:**
```julia
# Should show no allocations in steady state
@btime push_state!(\$segment, \$time, \$state)
```

### 7. Migration Strategy

**No backwards compatibility required** - breaking change acceptable since API is pre-release.

#### Package-Level Changes

**1. AstroModels (Core definitions)**
- Define `HistorySegment` and `SpacecraftHistory` structs
- Implement API functions (constructors, push_state!, push_segment!, queries)
- Update `Spacecraft` type: `history::Vector{...}` → `history::SpacecraftHistory`
- Remove old history accessor functions
- Export new types and functions

**2. AstroCallbacks (History recording)**
- Update callback functions to use `push_state!` instead of vector append
- Replace: `push!(sc.history[end], (time, vec(state)))` 
- With: `push_state!(current_segment(sc.history), time, state)`
- Add segment management if callbacks control propagation phases

**3. AstroProp (Propagation)**
- Call `new_segment!` at propagation start
- Update state recording to use `push_state!`
- Add segment naming/metadata for propagation phases
- Ensure coordinate system passed to segment constructor

**4. Epicycle (Visualization - View3D)**
- Update history access to iterate over segments
- Use `get_positions()` for 3D rendering per segment
- Iterate over segments if rendering with different styles
- Handle coordinate system transformations for multi-frame histories

**5. Examples & Tests**
- Update all example scripts using history
- Replace manual history indexing with API calls
- Add examples demonstrating segment-based workflows
- Update plotting examples to use `to_matrix()` or `get_positions()`

#### Implementation Order

```
Phase 1: Core (AstroModels)
├─ Define structs
├─ Implement API functions  
└─ Add tests

Phase 2: Recording (AstroCallbacks + AstroProp)
├─ Update propagation to create segments
├─ Update callbacks to push states
└─ Test basic propagation

Phase 3: Consumers (Epicycle + Examples)
├─ Update View3D rendering
├─ Update plotting utilities
├─ Update examples and docs
└─ Full integration tests

Phase 4: Cleanup
├─ Remove commented-out old code
├─ Update package documentation
└─ Verify all tests pass
```

#### Code Changes Required

**AstroModels/src/spacecraft.jl:**
```julia
# OLD
mutable struct Spacecraft{T}
    # ...
    history::Vector{Vector{Tuple{Time, Vector{Float64}}}}
end

# NEW  
mutable struct Spacecraft{T}
    # ...
    history::SpacecraftHistory
end

# Update constructor
Spacecraft(...) = Spacecraft(..., SpacecraftHistory(), ...)
```

**AstroCallbacks recording:**
```julia
# OLD
function record_callback(time, state, ...)
    push!(sc.history[end], (time, Vector{Float64}(vec(state))))
end

# NEW
function record_callback(time, state, ...)
    segment = last(sc.history.segments)
    push_state!(segment, time, CartesianState(state))
end
- ✓ Solver iteration history: Separate `iterations` vector with opt-in recording (Section 6.5)
- ✓ History recording control: `enabled` and `record_iterations` flags on SpacecraftHistory (Section 6.5)
```

**AstroProp propagation start:**
```julia
# OLD  
propagate!(sc, ...) = begin
    push!(sc.history, Vector{Tuple{Time,Vector{Float64}}}())
    # ...
end

# NEW
propagate!(sc, ...) = begin
    new_segment!(sc.history, sc.coordinate_system, 
                 name="propagate_\$(now())")
    # ...
end
```

**Epicycle View3D rendering:**
```julia
# OLD
positions = [hist[2] for seg in sc.history for hist in seg]

# NEW  
# Iterate per segment, handle coordinate systems
for segment in sc.history.segments
    positions = get_positions(segment)  # Nx3 matrix per segment
    # ... render with segment's coordinate system
end
```

#### File Locations

- `AstroModels/src/history.jl` - New file with struct definitions and API
- `AstroModels/src/spacecraft.jl` - Update Spacecraft type
- `AstroModels/src/AstroModels.jl` - Add exports
- `AstroCallbacks/src/recording.jl` - Update recording callbacks
- `AstroProp/src/propagate.jl` - Update propagation functions
- `Epicycle/src/view3d.jl` - Update visualization
- `examples/` - Update all example scripts

#### Migration Checklist

Strategy: 
Implement code, tests, and docs simultaneously as we go. 
1) Update AstroModels first, then address modules that integrate with AstroModels 
in this order AstroProp, Epicycle(View3D), AstroCallbacks

- [x] Implement HistorySegment, unit test, doc strings 
- [x] Implement SpacecraftHistory, unit test, and docstrings
- [x] Update Spacecraft type definition and update AstroModels test and docs
- [x] AstroProp integration
- [x] AstroManeuvers integration
- [x] AstroCallbacks integration (no changes needed)
- [x] Update Epicycle View3D access patterns
- [x] Update/add examples
- [ ] Write Spacecraft user guide (quick reference?, other user material)
- [ ] Run full test suite across all packages
- [ ] Remove old commented code

### 8. Testing Plan

Tests will be developed alongside implementation (TDD approach):
- Write unit test for each function as it's implemented
- Verify type stability, edge cases, and API contracts
- Add integration tests after each phase (propagate → record → query workflows)
- Benchmark performance-critical operations (push_state!, get_positions)

See Migration Checklist for test coverage strategy.

### 9. Open Questions

**Resolved:**
- ✓ Dual → Float64 conversion: Use `Float64.(to_vector(state))` via `to_float64()` helper (Section 6)
- ✓ Solver iteration history: Separate `iterations` vector with opt-in recording (Section 6.5)
- ✓ History recording control: `record_segments` and `record_iterations` flags on SpacecraftHistory (Section 6.5)

**To be addressed during implementation:**
- Error handling strategy (throw vs return codes)
- Time interpolation implementation (defer to Phase 2)
