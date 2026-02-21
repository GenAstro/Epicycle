# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    HistorySegment

Represents a continuous segment of spacecraft trajectory history.

A history segment stores time-ordered state data for a portion of a spacecraft's 
trajectory. Segments are used to delineate discontinuities (e.g., maneuvers) and 
to organize mission events (e.g., propagation phases).

# Fields
- `times::Vector{Time}`: Ordered vector of time points
- `states::Vector{CartesianState{Float64}}`: Corresponding state vectors (always Float64)
- `coordinate_system::CoordinateSystem`: Reference frame for this segment
- `name::String`: Optional segment identifier (empty string if unnamed)
- `metadata::Dict{String,Any}`: Extensible metadata storage for custom attributes

# Notes
- States are always stored as Float64 regardless of input type (e.g., Dual numbers)
- Immutable struct, but internal vectors can be mutated for efficiency
- `times` and `states` are parallel vectors of equal length

# Examples
```julia
using AstroModels, AstroUniverse, AstroEpochs, AstroFrames, AstroStates

# Create empty segment
coord_sys = CoordinateSystem(earth, ICRFAxes())
segment = HistorySegment(coord_sys, name="initial_orbit")

# Create segment from existing data
times = [Time("2024-01-01T00:00:00", TAI(), ISOT()), 
         Time("2024-01-01T01:00:00", TAI(), ISOT())]
states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
          CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
segment = HistorySegment(times, states, coord_sys, name="propagation_1")
```
"""
struct HistorySegment
    times::Vector{Time{Float64}}
    states::Vector{CartesianState{Float64}}
    coordinate_system::CoordinateSystem
    name::String
    metadata::Dict{String,Any}
end

"""
    HistorySegment(coord_system::CoordinateSystem; 
                   name::String="", 
                   metadata::Dict{String,Any}=Dict{String,Any}())

Create an empty history segment with specified coordinate system.
"""
function HistorySegment(coord_system::CoordinateSystem; 
                       name::String="",
                       metadata::Dict{String,Any}=Dict{String,Any}())
    HistorySegment(Vector{Time{Float64}}(), 
                   Vector{CartesianState{Float64}}(),
                   coord_system,
                   name,
                   metadata)
end

"""
    HistorySegment(times::Vector{<:Time}, 
                   states::Vector{CartesianState{Float64}}, 
                   coord_system::CoordinateSystem;
                   name::String="",
                   metadata::Dict{String,Any}=Dict{String,Any}())

Create a history segment from existing time and state data.
"""
function HistorySegment(times::Vector{T},
                       states::Vector{S},
                       coord_system::CoordinateSystem;
                       name::String="",
                       metadata::Dict{String,Any}=Dict{String,Any}()) where {T<:Time, S<:CartesianState}
    @assert length(times) == length(states) "times and states must have equal length"
    times_f64 = [to_float64(t) for t in times]
    states_f64 = [to_float64(s) for s in states]
    HistorySegment(times_f64, states_f64, coord_system, name, metadata)
end

# Base interface methods
Base.isempty(segment::HistorySegment) = isempty(segment.times)

# =============================================================================
# Utility Functions
# =============================================================================

"""
    _value_to_float64(x::Real)

Helper to safely convert any Real to Float64, handling Dual numbers from AD.
"""
function _value_to_float64(x::Real)
    # For regular numbers
    if x isa AbstractFloat || x isa Integer
        return Float64(x)
    # For ForwardDiff.Dual and other AD types with .value field
    elseif hasfield(typeof(x), :value)
        return Float64(x.value)
    else
        # Fallback: try to convert directly
        return Float64(x)
    end
end

"""
    to_float64(time::Time{Float64})
    to_float64(time::Time)

Convert a Time of any numeric type to Time{Float64}.

This function handles automatic conversion from automatic differentiation types 
(Dual numbers), BigFloat, and other numeric types to Float64 for history storage.

# Arguments
- `time::Time`: Time to convert (any numeric type parameter)

# Returns
- `Time{Float64}`: Time with Float64 components

# Notes
- No-op for Time{Float64} (returns input unchanged)
- Preserves time scale and format during conversion

# Examples
```julia
# From Dual numbers (automatic differentiation)
time_dual = Time{Dual}(...)
time_f64 = to_float64(time_dual)

# From Float64 (no-op)
time = Time("2024-01-01T00:00:00", TAI(), ISOT())
same_time = to_float64(time)  # Returns input unchanged
```
"""
to_float64(t::Time{Float64}) = t  # No-op for Float64
to_float64(t::Time) = Time(_value_to_float64(t.jd1), _value_to_float64(t.jd2), getfield(t, :scale), getfield(t, :format))

"""
    to_float64(state::CartesianState{Float64})
    to_float64(state::CartesianState)

Convert a CartesianState of any numeric type to CartesianState{Float64}.

This function handles automatic conversion from automatic differentiation types 
(Dual numbers), BigFloat, and other numeric types to Float64 for history storage.

# Arguments
- `state::CartesianState`: State to convert (any numeric type parameter)

# Returns
- `CartesianState{Float64}`: State with Float64 components

# Notes
- No-op for CartesianState{Float64} (returns input unchanged)
- Uses `to_vector()` and broadcasts Float64 conversion

# Examples
```julia
# From Dual numbers (automatic differentiation)
state_dual = CartesianState{Dual}([...])
state_f64 = to_float64(state_dual)

# From BigFloat
state_big = CartesianState{BigFloat}([...])
state_f64 = to_float64(state_big)

# From Float64 (no-op)
state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
same_state = to_float64(state)  
```
"""
to_float64(state::CartesianState{Float64}) = state  
to_float64(state::CartesianState) = CartesianState(_value_to_float64.(to_vector(state)))

# =============================================================================
# Display Methods
# =============================================================================

"""
    Base.show(io::IO, segment::HistorySegment)

Pretty-print a HistorySegment in a human-readable, multi-line summary.
"""
function Base.show(io::IO, segment::HistorySegment)
    # Segment name (or empty string if unnamed)
    name_str = isempty(segment.name) ? "" : segment.name
    println(io, "HistorySegment: \"", name_str, "\"")
    
    # Number of points
    n_points = length(segment.times)
    println(io, "  Points: ", n_points)
    
    # Time range (if data exists) - compact format
    if n_points > 0
        # Get compact time strings (first line only from Time's show)
        start_time = segment.times[1]
        end_time = segment.times[end]
        
        # Use sprint to get compact representation
        start_io = IOBuffer()
        print(start_io, start_time)
        start_lines = split(String(take!(start_io)), '\n')
        start_value = strip(split(start_lines[2], '=')[2])  
        
        end_io = IOBuffer()
        print(end_io, end_time)
        end_lines = split(String(take!(end_io)), '\n')
        end_value = strip(split(end_lines[2], '=')[2])
        
        # Format scale (handle symbol -> uppercase string)
        scale_str = uppercase(string(start_time.scale))
        
        println(io, "  Start Time: ", start_value, " ", scale_str)
        println(io, "  End Time: ", end_value, " ", scale_str)
    else
        println(io, "  Start Time: (no data)")
        println(io, "  End Time: (no data)")
    end
    
    # Coordinate system (just origin name and axes type name)
    origin = segment.coordinate_system.origin
    origin_name = hasfield(typeof(origin), :name) ? origin.name : string(origin)
    axes_name = split(string(typeof(segment.coordinate_system.axes)), ".")[end]  # Just type name
    println(io, "  Coordinate System: ", origin_name, ", ", axes_name)
    
    # Metadata count
    n_metadata = length(segment.metadata)
    print(io, "  Metadata: ", n_metadata, n_metadata == 1 ? " entry" : " entries")
end

#==============================================================================
# SpacecraftHistory
==============================================================================#

"""
    SpacecraftHistory

Container for spacecraft trajectory history, organized into segments.

A SpacecraftHistory stores the complete trajectory of a spacecraft as a collection 
of segments. Each segment represents a continuous portion of the 
trajectory, with segments typically separated by discontinuities (e.g., maneuvers) 
or mission phase boundaries.

# Fields
- `segments::Vector{HistorySegment}`: Solution trajectory segments
- `iterations::Vector{HistorySegment}`: Solver iteration segments (diagnostic data)
- `record_segments::Bool`: Enable solution segment recording (default: true)
- `record_iterations::Bool`: Enable iteration segment recording (default: false)

# Notes
- Mutable struct to allow solver to control recording flags
- Empty history is valid (zero segments)
- Provides iteration, indexing, and standard collection interfaces
- Iteration and indexing operate on `segments` only (solution trajectory)
- Access `iterations` directly for diagnostic data: `history.iterations`

# Examples
```julia
using AstroModels, AstroEpochs, AstroUniverse, AstroFrames

# Basic usage - record solution trajectory (default)
history = SpacecraftHistory()
coord_sys = CoordinateSystem(earth, ICRFAxes())

segment1 = HistorySegment(coord_sys, name="orbit_1")
push_segment!(history, segment1)

# Iterate over solution segments
for (i, segment) in enumerate(history)
    println("Segment ", i, ": ", segment.name, " (", length(segment.times), " points)")
end

# Access by index (operates on solution segments only)
first_segment = history[1]
println("History contains ", length(history), " solution segments")

# Advanced: Solver iteration recording for diagnostics
# Disable solution recording, enable iteration recording
history.record_segments = false
history.record_iterations = true

# During optimization, segments go to iterations vector
for iter in 1:10
    segment = HistorySegment(coord_sys, name="iteration_\$iter")
    push_segment!(history, segment)  # Routes to history.iterations
end

# Switch back for final solution
history.record_segments = true
history.record_iterations = false
final_segment = HistorySegment(coord_sys, name="solution")
push_segment!(history, final_segment)  # Routes to history.segments

# Result: 1 solution segment, 10 diagnostic iterations
println(length(history))             # 1 (solution only)
println(length(history.iterations))  # 10 (diagnostic data)
```
"""
mutable struct SpacecraftHistory
    segments::Vector{HistorySegment}
    iterations::Vector{HistorySegment}
    record_segments::Bool
    record_iterations::Bool
end

"""
    SpacecraftHistory()

Create an empty SpacecraftHistory with default settings.

# Defaults
- `record_segments = true`: Solution trajectory recording enabled
- `record_iterations = false`: Iteration recording disabled (diagnostic mode off)
"""
SpacecraftHistory() = SpacecraftHistory(Vector{HistorySegment}(), Vector{HistorySegment}(), true, false)

"""
    SpacecraftHistory(segments::Vector{HistorySegment})

Create a SpacecraftHistory from an existing vector of segments with default recording settings.

# Defaults
- `iterations`: Empty vector (no iteration data)
- `record_segments = true`: Solution trajectory recording enabled
- `record_iterations = false`: Iteration recording disabled
"""
SpacecraftHistory(segments::Vector{HistorySegment}) = SpacecraftHistory(segments, Vector{HistorySegment}(), true, false)

# Base interface methods
Base.length(h::SpacecraftHistory) = length(h.segments)
Base.isempty(h::SpacecraftHistory) = isempty(h.segments)
Base.getindex(h::SpacecraftHistory, i::Int) = h.segments[i]
Base.getindex(h::SpacecraftHistory, r::UnitRange) = h.segments[r]
Base.iterate(h::SpacecraftHistory, state=1) = 
    state > length(h.segments) ? nothing : (h.segments[state], state+1)
Base.firstindex(h::SpacecraftHistory) = 1
Base.lastindex(h::SpacecraftHistory) = length(h.segments)

"""
    Base.copy(history::SpacecraftHistory)

Create a deep copy of a SpacecraftHistory, copying all segments and iterations.
"""
function Base.copy(history::SpacecraftHistory)
    SpacecraftHistory(copy(history.segments), copy(history.iterations), 
                     history.record_segments, history.record_iterations)
end

"""
    Base.deepcopy_internal(history::SpacecraftHistory, dict::IdDict)

Deep copy implementation for SpacecraftHistory to ensure proper copying of mutable fields.
"""
function Base.deepcopy_internal(history::SpacecraftHistory, dict::IdDict)
    SpacecraftHistory(
        Base.deepcopy_internal(history.segments, dict),
        Base.deepcopy_internal(history.iterations, dict),
        history.record_segments,
        history.record_iterations
    )
end

"""
    push_segment!(history::SpacecraftHistory, segment::HistorySegment)

Add a segment to the history, routing based on recording flags.

# Routing Logic
- If `record_iterations=true`: adds to `history.iterations` (diagnostic data)
- Else if `record_segments=true`: adds to `history.segments` (solution trajectory)
- Else: no-op (segment not recorded)

# Notes
This routing is controlled by the solver during trajectory optimization:
- During optimization iterations: `record_iterations=true, record_segments=false`
- Final solution run: `record_iterations=false, record_segments=true`
"""
function push_segment!(history::SpacecraftHistory, segment::HistorySegment)
    if history.record_iterations
        push!(history.iterations, segment)
    elseif history.record_segments
        push!(history.segments, segment)
    end
    # If both false: no-op (segment not recorded)
    return nothing
end

"""
    Base.show(io::IO, history::SpacecraftHistory)

Display a SpacecraftHistory in a human-readable format.
"""
function Base.show(io::IO, history::SpacecraftHistory)
    n_segments = length(history)
    n_points = sum(length(seg.times) for seg in history.segments; init=0)
    
    println(io, "SpacecraftHistory with ", n_segments, 
            n_segments == 1 ? " segment" : " segments", 
            " (", n_points, " total points)")
    
    # Show iteration count if present
    n_iterations = length(history.iterations)
    if n_iterations > 0
        n_iter_points = sum(length(seg.times) for seg in history.iterations; init=0)
        println(io, "  Iterations: ", n_iterations, 
               n_iterations == 1 ? " segment" : " segments",
               " (", n_iter_points, " points, diagnostic data)")
    end
    
    if n_segments > 0
        for (i, segment) in enumerate(history)
            seg_points = length(segment.times)
            seg_name = isempty(segment.name) ? "(unnamed)" : segment.name
            println(io, "  [", i, "] ", seg_name, ": ", seg_points, 
                   seg_points == 1 ? " point" : " points")
        end
    end
end
