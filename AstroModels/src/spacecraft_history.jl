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
using AstroModels, AstroEpochs, AstroFrames, AstroStates

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
function HistorySegment(times::Vector{<:Time},
                       states::Vector{CartesianState{Float64}},
                       coord_system::CoordinateSystem;
                       name::String="",
                       metadata::Dict{String,Any}=Dict{String,Any}())
    @assert length(times) == length(states) "times and states must have equal length"
    times_f64 = [to_float64(t) for t in times]
    HistorySegment(times_f64, states, coord_system, name, metadata)
end

# Base interface methods
Base.isempty(segment::HistorySegment) = isempty(segment.times)

# =============================================================================
# Utility Functions
# =============================================================================

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
to_float64(t::Time) = Time(Float64(t.jd1), Float64(t.jd2), getfield(t, :scale), getfield(t, :format))

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
to_float64(state::CartesianState) = CartesianState(Float64.(to_vector(state)))

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
        start_value = strip(split(start_lines[2], '=')[2])  # "value = 2024-01-01..." -> "2024-01-01..."
        
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
of time-ordered segments. Each segment represents a continuous portion of the 
trajectory, with segments typically separated by discontinuities (e.g., maneuvers) 
or mission phase boundaries.

# Fields
- `segments::Vector{HistorySegment}`: Ordered list of trajectory segments

# Notes
- Immutable struct, but segments vector can be mutated via push_segment!
- Empty history is valid (zero segments)
- Provides iteration, indexing, and standard collection interfaces

# Examples
```julia
using AstroModels, AstroEpochs, AstroFrames

# Create empty history
history = SpacecraftHistory()

# Add segments during propagation
coord_sys = CoordinateSystem(earth, ICRFAxes())
segment1 = HistorySegment(coord_sys, name="orbit_1")
push_segment!(history, segment1)

# Iterate over segments
for (i, segment) in enumerate(history)
    println("Segment ", i, ": ", segment.name, " (", length(segment.times), " points)")
end

# Access by index
first_segment = history[1]

# Check if empty
if !isempty(history)
    println("History contains ", length(history), " segments")
end
```
"""
struct SpacecraftHistory
    segments::Vector{HistorySegment}
end

"""
    SpacecraftHistory()

Create an empty SpacecraftHistory with no segments.
"""
SpacecraftHistory() = SpacecraftHistory(Vector{HistorySegment}())

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
    push_segment!(history::SpacecraftHistory, segment::HistorySegment)

Add a segment to the history.
"""
function push_segment!(history::SpacecraftHistory, segment::HistorySegment)
    push!(history.segments, segment)
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
    
    if n_segments > 0
        for (i, segment) in enumerate(history)
            seg_points = length(segment.times)
            seg_name = isempty(segment.name) ? "(unnamed)" : segment.name
            println(io, "  [", i, "] ", seg_name, ": ", seg_points, 
                   seg_points == 1 ? " point" : " points")
        end
    end
end
