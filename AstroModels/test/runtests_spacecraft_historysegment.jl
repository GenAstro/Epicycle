# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using Test
using AstroModels
using AstroEpochs
using AstroFrames
using AstroStates
using AstroUniverse

# =============================================================================
# HistorySegment - Empty Constructor
# =============================================================================

@testset "SpacecraftHistory Tests" begin
@testset "HistorySegment - empty constructor basic" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys)
    
    @test segment.times isa Vector{Time{Float64}}
    @test isempty(segment.times)
    @test segment.states isa Vector{CartesianState{Float64}}
    @test isempty(segment.states)
    @test segment.coordinate_system === coord_sys
    @test segment.name == ""
    @test segment.metadata isa Dict{String,Any}
    @test isempty(segment.metadata)
end

@testset "HistorySegment - empty constructor with name" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys, name="test_segment")
    
    @test segment.name == "test_segment"
    @test isempty(segment.times)
    @test isempty(segment.states)
    @test isempty(segment.metadata)
end

@testset "HistorySegment - empty constructor with metadata" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    meta = Dict("mission" => "test", "phase" => 1)
    segment = HistorySegment(coord_sys, metadata=meta)
    
    @test segment.metadata["mission"] == "test"
    @test segment.metadata["phase"] == 1
    @test segment.name == ""
end

@testset "HistorySegment - empty constructor with all kwargs" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    meta = Dict("event" => "maneuver", "delta_v" => 0.5)
    segment = HistorySegment(coord_sys, name="maneuver_1", metadata=meta)
    
    @test segment.name == "maneuver_1"
    @test segment.metadata["event"] == "maneuver"
    @test segment.metadata["delta_v"] == 0.5
end

# =============================================================================
# HistorySegment - Full Constructor
# =============================================================================

@testset "HistorySegment - full constructor basic" begin
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT()),
                Time("2024-01-01T01:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    segment = HistorySegment(times, states, coord_sys)
    
    @test length(segment.times) == 2
    @test length(segment.states) == 2
    @test segment.times[1] == times[1]
    @test segment.times[2] == times[2]
    @test segment.states[1] == states[1]
    @test segment.states[2] == states[2]
    @test segment.coordinate_system === coord_sys
    @test segment.name == ""
    @test isempty(segment.metadata)
end

@testset "HistorySegment - full constructor with kwargs" begin
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    meta = Dict{String,Any}("source" => "propagation")
    
    segment = HistorySegment(times, states, coord_sys, name="orbit1", metadata=meta)
    
    @test segment.name == "orbit1"
    @test segment.metadata["source"] == "propagation"
end

@testset "HistorySegment - full constructor single point" begin
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    segment = HistorySegment(times, states, coord_sys)
    
    @test length(segment.times) == 1
    @test length(segment.states) == 1
end

@testset "HistorySegment - full constructor many points" begin
    n = 100
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT()) + i*60.0 for i in 0:n-1]
    states = [CartesianState([7000.0 + i, 0.0, 0.0, 0.0, 7.5, 0.0]) for i in 0:n-1]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    segment = HistorySegment(times, states, coord_sys, name="long_prop")
    
    @test length(segment.times) == n
    @test length(segment.states) == n
    @test segment.name == "long_prop"
end

@testset "HistorySegment - full constructor empty vectors" begin
    times = Vector{Time}()
    states = Vector{CartesianState{Float64}}()
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    segment = HistorySegment(times, states, coord_sys)
    
    @test isempty(segment.times)
    @test isempty(segment.states)
end

# =============================================================================
# HistorySegment - Validation
# =============================================================================

@testset "HistorySegment - mismatched lengths fail" begin
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    @test_throws AssertionError HistorySegment(times, states, coord_sys)
end

@testset "HistorySegment - mismatched lengths fail (reversed)" begin
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT()),
                Time("2024-01-01T01:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    @test_throws AssertionError HistorySegment(times, states, coord_sys)
end

# =============================================================================
# HistorySegment - Base Interface Methods
# =============================================================================

@testset "HistorySegment - isempty" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    
    # Empty segment
    empty_seg = HistorySegment(coord_sys)
    @test isempty(empty_seg)
    
    # Non-empty segment
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    non_empty_seg = HistorySegment(times, states, coord_sys)
    @test !isempty(non_empty_seg)
end

# =============================================================================
# HistorySegment - Immutability
# =============================================================================

@testset "HistorySegment - struct is immutable" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys, name="test")
    
    # Cannot reassign struct fields
    @test_throws ErrorException (segment.name = "new_name")
end

@testset "HistorySegment - internal vectors are mutable" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys)
    
    # Can mutate internal vectors
    time1 = Time("2024-01-01T00:00:00", TAI(), ISOT())
    state1 = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
    
    push!(segment.times, time1)
    push!(segment.states, state1)
    
    @test length(segment.times) == 1
    @test length(segment.states) == 1
    @test segment.times[1] == time1
    @test segment.states[1] == state1
end

@testset "HistorySegment - metadata dict is mutable" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys)
    
    # Can add to metadata dict
    segment.metadata["key1"] = "value1"
    segment.metadata["key2"] = 42
    
    @test segment.metadata["key1"] == "value1"
    @test segment.metadata["key2"] == 42
end

# =============================================================================
# to_float64() - Float64 No-op
# =============================================================================

@testset "to_float64 - Float64 no-op" begin
    state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
    result = to_float64(state)
    
    @test result === state  # Same object
    @test result isa CartesianState{Float64}
end

@testset "to_float64 - Float64 values unchanged" begin
    state = CartesianState([7123.456, 234.567, 345.678, 1.234, 7.567, 0.089])
    result = to_float64(state)
    
    @test result === state
    @test result.position[1] === 7123.456
    @test result.velocity[3] === 0.089
end

# =============================================================================
# to_float64() - Type Conversions
# =============================================================================

@testset "to_float64 - Int conversion" begin
    state_int = CartesianState([7000, 0, 0, 0, 7, 0])
    result = to_float64(state_int)
    
    @test result isa CartesianState{Float64}
    @test result.position[1] == 7000.0
    @test result.velocity[2] == 7.0
    @test result.velocity[3] == 0.0
end

@testset "to_float64 - BigFloat conversion" begin
    state_big = CartesianState([BigFloat(7000), BigFloat(300), BigFloat(0),
                                BigFloat(0), BigFloat(7.5), BigFloat(0.03)])
    result = to_float64(state_big)
    
    @test result isa CartesianState{Float64}
    @test result.position[1] ≈ 7000.0
    @test result.position[2] ≈ 300.0
    @test result.velocity[2] ≈ 7.5
    @test result.velocity[3] ≈ 0.03
end

@testset "to_float64 - Float32 conversion" begin
    state_f32 = CartesianState(Float32.([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]))
    result = to_float64(state_f32)
    
    @test result isa CartesianState{Float64}
    @test result.position[1] ≈ 7000.0
    @test result.velocity[2] ≈ 7.5
end

@testset "to_float64 - mixed Int and Float" begin
    # CartesianState constructor should promote types
    state = CartesianState([7000, 0, 0, 0.0, 7.5, 0])
    result = to_float64(state)
    
    @test result isa CartesianState{Float64}
    @test result.position[1] == 7000.0
end

# =============================================================================
# to_float64() - Precision Preservation
# =============================================================================

@testset "to_float64 - preserves Float64 precision" begin
    # Use values that test precision
    state = CartesianState([7000.123456789, 300.987654321, 0.0,
                            0.0, 7.512345678, 0.034567890])
    result = to_float64(state)
    
    @test result.position[1] === 7000.123456789
    @test result.position[2] === 300.987654321
    @test result.velocity[2] === 7.512345678
    @test result.velocity[3] === 0.034567890
end

# =============================================================================
# Integration Tests
# =============================================================================

@testset "Integration - build segment with mixed types" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys, name="mixed_types")
    
    # Add Float64 state
    push!(segment.times, Time("2024-01-01T00:00:00", TAI(), ISOT()))
    push!(segment.states, to_float64(CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])))
    
    # Add Int state (converted)
    push!(segment.times, Time("2024-01-01T01:00:00", TAI(), ISOT()))
    push!(segment.states, to_float64(CartesianState([7100, 100, 0, 0, 7, 0])))
    
    # Add BigFloat state (converted)
    push!(segment.times, Time("2024-01-01T02:00:00", TAI(), ISOT()))
    push!(segment.states, to_float64(CartesianState([BigFloat(7200), BigFloat(200),
                                                        BigFloat(0), BigFloat(0),
                                                        BigFloat(7.2), BigFloat(0)])))
    
    @test length(segment.states) == 3
    @test all(s -> s isa CartesianState{Float64}, segment.states)
    @test segment.states[1].position[1] == 7000.0
    @test segment.states[2].position[1] == 7100.0
    @test segment.states[3].position[1] ≈ 7200.0
end

@testset "Integration - segment with coordinate system" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    
    segment = HistorySegment(times, states, coord_sys, name="test")
    
    @test segment.coordinate_system.origin === earth
    @test segment.coordinate_system.axes isa ICRFAxes
end

# =============================================================================
# Show Method Tests
# =============================================================================

@testset "Show - empty segment" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys)
    
    output = sprint(show, segment)
    
    @test contains(output, "HistorySegment: \"\"")
    @test contains(output, "Points: 0")
    @test contains(output, "Start Time: (no data)")
    @test contains(output, "End Time: (no data)")
    @test contains(output, "Coordinate System: Earth, ICRFAxes")
    @test contains(output, "Metadata: 0 entries")
end

@testset "Show - segment with name" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    segment = HistorySegment(coord_sys, name="test_orbit")
    
    output = sprint(show, segment)
    
    @test contains(output, "HistorySegment: \"test_orbit\"")
end

@testset "Show - segment with single point" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
    segment = HistorySegment(times, states, coord_sys)
    
    output = sprint(show, segment)
    
    @test contains(output, "Points: 1")
    @test contains(output, "Start Time: 2024-01-01T00:00:00.000")
    @test contains(output, "End Time: 2024-01-01T00:00:00.000")
    @test contains(output, "TAI")
end

@testset "Show - segment with multiple points" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT()),
             Time("2024-01-01T01:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
              CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
    segment = HistorySegment(times, states, coord_sys)
    
    output = sprint(show, segment)
    
    @test contains(output, "Points: 2")
    @test contains(output, "Start Time: 2024-01-01T00:00:00.000")
    @test contains(output, "End Time: 2024-01-01T01:00:00.000")
end

@testset "Show - segment with metadata" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    meta = Dict{String,Any}("mission" => "demo", "phase" => 1, "dv" => 0.5)
    segment = HistorySegment(coord_sys, metadata=meta)
    
    output = sprint(show, segment)
    
    @test contains(output, "Metadata: 3 entries")
end

@testset "Show - segment with one metadata entry" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    meta = Dict{String,Any}("mission" => "demo")
    segment = HistorySegment(coord_sys, metadata=meta)
    
    output = sprint(show, segment)
    
    @test contains(output, "Metadata: 1 entry")
end

@testset "Show - complete segment" begin
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    times = [Time("2024-01-01T00:00:00", TAI(), ISOT()),
             Time("2024-01-01T06:00:00", TAI(), ISOT())]
    states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
              CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
    meta = Dict{String,Any}("event" => "propagation")
    segment = HistorySegment(times, states, coord_sys, name="orbit_1", metadata=meta)
    
    output = sprint(show, segment)
    
    @test contains(output, "HistorySegment: \"orbit_1\"")
    @test contains(output, "Points: 2")
    @test contains(output, "Start Time: 2024-01-01T00:00:00.000")
    @test contains(output, "End Time: 2024-01-01T06:00:00.000")
    @test contains(output, "TAI")
    @test contains(output, "Coordinate System: Earth, ICRFAxes")
    @test contains(output, "Metadata: 1 entry")
end

end

