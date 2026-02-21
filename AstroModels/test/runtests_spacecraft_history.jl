# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

using Test
using AstroModels
using AstroEpochs
using AstroStates
using AstroFrames
using AstroUniverse

@testset "SpacecraftHistory Tests" begin

    #==========================================================================
    # Setup - Reusable test data
    ==========================================================================#
    
    coord_sys = CoordinateSystem(earth, ICRFAxes())
    t1 = Time("2024-01-01T00:00:00", TAI(), ISOT())
    t2 = Time("2024-01-01T01:00:00", TAI(), ISOT())
    state1 = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
    state2 = CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])
    
    #==========================================================================
    # Constructor Tests
    ==========================================================================#
    
    @testset "Empty Constructor" begin
        history = SpacecraftHistory()
        
        @test isa(history, SpacecraftHistory)
        @test length(history) == 0
        @test isempty(history)
        @test isa(history.segments, Vector{HistorySegment})
    end
    
    @testset "Constructor with Segments Vector" begin
        segment1 = HistorySegment(coord_sys, name="seg1")
        segment2 = HistorySegment(coord_sys, name="seg2")
        segments = [segment1, segment2]
        
        history = SpacecraftHistory(segments)
        
        @test length(history) == 2
        @test !isempty(history)
        @test history.segments === segments
    end
    
    #==========================================================================
    # Base Interface Tests
    ==========================================================================#
    
    @testset "Length and Empty" begin
        history = SpacecraftHistory()
        @test length(history) == 0
        @test isempty(history)
        
        seg1 = HistorySegment(coord_sys)
        push_segment!(history, seg1)
        @test length(history) == 1
        @test !isempty(history)
        
        seg2 = HistorySegment(coord_sys)
        push_segment!(history, seg2)
        @test length(history) == 2
        @test !isempty(history)
    end
    
    @testset "Indexing - Single Element" begin
        seg1 = HistorySegment(coord_sys, name="first")
        seg2 = HistorySegment(coord_sys, name="second")
        seg3 = HistorySegment(coord_sys, name="third")
        history = SpacecraftHistory([seg1, seg2, seg3])
        
        @test history[1] === seg1
        @test history[2] === seg2
        @test history[3] === seg3
        @test history[1].name == "first"
        @test history[end].name == "third"
    end
    
    @testset "Indexing - Range" begin
        segs = [HistorySegment(coord_sys, name="seg$i") for i in 1:5]
        history = SpacecraftHistory(segs)
        
        subset = history[2:4]
        @test length(subset) == 3
        @test subset[1].name == "seg2"
        @test subset[3].name == "seg4"
        
        all_segs = history[1:end]
        @test length(all_segs) == 5
    end
    
    @testset "First and Last Index" begin
        segs = [HistorySegment(coord_sys) for _ in 1:3]
        history = SpacecraftHistory(segs)
        
        @test firstindex(history) == 1
        @test lastindex(history) == 3
        @test history[firstindex(history)] === segs[1]
        @test history[lastindex(history)] === segs[3]
    end
    
    @testset "Iteration" begin
        seg1 = HistorySegment(coord_sys, name="a")
        seg2 = HistorySegment(coord_sys, name="b")
        seg3 = HistorySegment(coord_sys, name="c")
        history = SpacecraftHistory([seg1, seg2, seg3])
        
        # Collect via iteration
        collected = collect(history)
        @test length(collected) == 3
        @test collected[1] === seg1
        @test collected[2] === seg2
        @test collected[3] === seg3
        
        # Enumerate
        names = String[]
        for (i, segment) in enumerate(history)
            push!(names, segment.name)
            @test i in 1:3
        end
        @test names == ["a", "b", "c"]
        
        # Empty iteration
        empty_hist = SpacecraftHistory()
        count = 0
        for _ in empty_hist
            count += 1
        end
        @test count == 0
    end
    
    #==========================================================================
    # push_segment! Tests
    ==========================================================================#
    
    @testset "push_segment! - Basic" begin
        history = SpacecraftHistory()
        seg1 = HistorySegment(coord_sys, name="first")
        
        result = push_segment!(history, seg1)
        
        @test result === nothing
        @test length(history) == 1
        @test history[1] === seg1
    end
    
    @testset "push_segment! - Multiple" begin
        history = SpacecraftHistory()
        
        for i in 1:5
            seg = HistorySegment(coord_sys, name="seg$i")
            push_segment!(history, seg)
        end
        
        @test length(history) == 5
        @test history[1].name == "seg1"
        @test history[5].name == "seg5"
    end
    
    @testset "push_segment! - With Data" begin
        history = SpacecraftHistory()
        
        times = [t1, t2]
        states = [state1, state2]
        seg = HistorySegment(times, states, coord_sys, name="data_seg")
        push_segment!(history, seg)
        
        @test length(history) == 1
        @test length(history[1].times) == 2
        @test history[1].states[1] == state1
    end
    
    #==========================================================================
    # Mutability Tests
    ==========================================================================#
    
    @testset "Vector Mutability" begin
        seg1 = HistorySegment(coord_sys, name="original")
        history = SpacecraftHistory([seg1])
        
        # Can mutate internal vector
        seg2 = HistorySegment(coord_sys, name="added")
        push!(history.segments, seg2)
        
        @test length(history) == 2
        @test history[2] === seg2
    end
    
    #==========================================================================
    # Integration Tests
    ==========================================================================#
    
    @testset "Multi-Segment Workflow" begin
        history = SpacecraftHistory()
        
        # Segment 1: Initial orbit
        seg1 = HistorySegment(coord_sys, name="initial")
        push_segment!(history, seg1)
        
        # Segment 2: Propagation
        times = [t1, t2]
        states = [state1, state2]
        seg2 = HistorySegment(times, states, coord_sys, name="propagate")
        push_segment!(history, seg2)
        
        # Segment 3: After maneuver
        seg3 = HistorySegment(coord_sys, name="post_maneuver")
        push_segment!(history, seg3)
        
        @test length(history) == 3
        @test history[1].name == "initial"
        @test history[2].name == "propagate"
        @test history[3].name == "post_maneuver"
        @test length(history[2].times) == 2
    end
    
    @testset "Total Points Calculation" begin
        history = SpacecraftHistory()
        
        # Empty history
        total = sum(length(seg.times) for seg in history.segments; init=0)
        @test total == 0
        
        # Add segments with varying point counts
        seg1 = HistorySegment([t1], [state1], coord_sys)
        seg2 = HistorySegment([t1, t2], [state1, state2], coord_sys)
        seg3 = HistorySegment(coord_sys)  # Empty
        
        push_segment!(history, seg1)
        push_segment!(history, seg2)
        push_segment!(history, seg3)
        
        total = sum(length(seg.times) for seg in history.segments; init=0)
        @test total == 3
    end
    
    #==========================================================================
    # Show Method Tests
    ==========================================================================#
    
    @testset "Show - Empty History" begin
        history = SpacecraftHistory()
        str = sprint(show, history)
        
        @test contains(str, "SpacecraftHistory")
        @test contains(str, "0 segments")
        @test contains(str, "0 total points")
    end
    
    @testset "Show - Single Segment" begin
        seg = HistorySegment(coord_sys, name="test_segment")
        history = SpacecraftHistory([seg])
        str = sprint(show, history)
        
        @test contains(str, "1 segment")
        @test contains(str, "0 total points")
        @test contains(str, "test_segment")
        @test contains(str, "[1]")
    end
    
    @testset "Show - Multiple Segments" begin
        seg1 = HistorySegment([t1], [state1], coord_sys, name="seg1")
        seg2 = HistorySegment([t1, t2], [state1, state2], coord_sys, name="seg2")
        seg3 = HistorySegment(coord_sys)  # unnamed
        history = SpacecraftHistory([seg1, seg2, seg3])
        
        str = sprint(show, history)
        
        @test contains(str, "3 segments")
        @test contains(str, "3 total points")
        @test contains(str, "[1] seg1: 1 point")
        @test contains(str, "[2] seg2: 2 points")
        @test contains(str, "[3] (unnamed): 0 points")
    end
    
    @testset "Show - Pluralization" begin
        # Test singular forms
        seg_single = HistorySegment([t1], [state1], coord_sys)
        hist_single = SpacecraftHistory([seg_single])
        str_single = sprint(show, hist_single)
        @test contains(str_single, "1 segment")
        @test contains(str_single, "1 point")
        
        # Test plural forms
        seg_plural = HistorySegment([t1, t2], [state1, state2], coord_sys)
        hist_plural = SpacecraftHistory([seg_plural, seg_plural])
        str_plural = sprint(show, hist_plural)
        @test contains(str_plural, "2 segments")
        @test contains(str_plural, "4 total points")
    end
    
    @testset "Show - Unnamed Segments" begin
        unnamed1 = HistorySegment(coord_sys)
        unnamed2 = HistorySegment(coord_sys, name="")
        history = SpacecraftHistory([unnamed1, unnamed2])
        
        str = sprint(show, history)
        @test occursin(r"\[1\] \(unnamed\)", str)
        @test occursin(r"\[2\] \(unnamed\)", str)
    end
    
    #==========================================================================
    # Edge Cases
    ==========================================================================#
    
    @testset "Large Number of Segments" begin
        history = SpacecraftHistory()
        
        for i in 1:100
            seg = HistorySegment(coord_sys, name="seg$i")
            push_segment!(history, seg)
        end
        
        @test length(history) == 100
        @test history[1].name == "seg1"
        @test history[50].name == "seg50"
        @test history[100].name == "seg100"
    end
    
    @testset "Mixed Segment Types" begin
        # Empty, small, and large segments
        empty_seg = HistorySegment(coord_sys, name="empty")
        small_seg = HistorySegment([t1], [state1], coord_sys, name="small")
        
        large_times = [Time("2024-01-01T$(lpad(i,2,'0')):00:00", TAI(), ISOT()) for i in 0:23]
        large_states = [CartesianState([7000.0+i*10, 0.0, 0.0, 0.0, 7.5, 0.0]) for i in 0:23]
        large_seg = HistorySegment(large_times, large_states, coord_sys, name="large")
        
        history = SpacecraftHistory([empty_seg, small_seg, large_seg])
        
        @test length(history) == 3
        @test length(history[1].times) == 0
        @test length(history[2].times) == 1
        @test length(history[3].times) == 24
    end
    
    #==========================================================================
    # Spacecraft Integration Tests
    ==========================================================================#
    
    @testset "Spacecraft Default History" begin
        # Test 1: Empty history on new spacecraft
        sc = Spacecraft(
            state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time=Time("2024-01-01T00:00:00", TAI(), ISOT()),
            mass=1000.0
        )
        
        @test sc.history isa SpacecraftHistory
        @test length(sc.history) == 0
        @test isempty(sc.history)
        @test length(sc.history.segments) == 0
    end
    
    @testset "Spacecraft Multiple Segments" begin
        # Test 3: Multiple segments
        sc = Spacecraft(
            state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time=Time("2024-01-01T00:00:00", TAI(), ISOT()),
            mass=1000.0
        )
        
        # Add first segment
        coord_sys1 = CoordinateSystem(earth, ICRFAxes())
        seg1_times = [Time("2024-01-01T00:00:00", TAI(), ISOT())]
        seg1_states = [CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])]
        seg1 = HistorySegment(seg1_times, seg1_states, coord_sys1, name="propagate_1")
        push_segment!(sc.history, seg1)
        
        # Add second segment (different coordinate system - spacecraft can change frames)
        coord_sys2 = CoordinateSystem(earth, ICRFAxes())  # Could be different frame
        seg2_times = [Time("2024-01-01T01:00:00", TAI(), ISOT())]
        seg2_states = [CartesianState([7100.0, 100.0, 0.0, 0.0, 7.4, 0.1])]
        seg2 = HistorySegment(seg2_times, seg2_states, coord_sys2, name="maneuver_1")
        push_segment!(sc.history, seg2)
        
        # Add third segment
        seg3_times = [Time("2024-01-01T02:00:00", TAI(), ISOT())]
        seg3_states = [CartesianState([7200.0, 200.0, 0.0, 0.0, 7.3, 0.2])]
        seg3 = HistorySegment(seg3_times, seg3_states, coord_sys2, name="propagate_2")
        push_segment!(sc.history, seg3)
        
        @test length(sc.history) == 3
        @test sc.history[1].name == "propagate_1"
        @test sc.history[2].name == "maneuver_1"
        @test sc.history[3].name == "propagate_2"
        @test length(sc.history.segments) == 3
    end
    
    @testset "Spacecraft History Iteration" begin
        # Test 4: History iteration
        sc = Spacecraft(
            state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time=Time("2024-01-01T00:00:00", TAI(), ISOT()),
            mass=1000.0
        )
        
        # Add multiple segments
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        for i in 1:5
            seg_times = [Time("2024-01-01T$(lpad(i-1,2,'0')):00:00", TAI(), ISOT())]
            seg_states = [CartesianState([7000.0+i*100, 0.0, 0.0, 0.0, 7.5, 0.0])]
            seg = HistorySegment(seg_times, seg_states, coord_sys, name="segment_$i")
            push_segment!(sc.history, seg)
        end
        
        # Test iteration with for loop
        segment_count = 0
        segment_names = String[]
        for segment in sc.history
            segment_count += 1
            push!(segment_names, segment.name)
        end
        
        @test segment_count == 5
        @test length(segment_names) == 5
        @test segment_names == ["segment_1", "segment_2", "segment_3", "segment_4", "segment_5"]
        
        # Test iteration with collect
        segments_collected = collect(sc.history)
        @test length(segments_collected) == 5
        @test all(seg isa HistorySegment for seg in segments_collected)
        
        # Test iteration with enumerate
        for (i, segment) in enumerate(sc.history)
            @test segment.name == "segment_$i"
        end
    end
    
    #==========================================================================
    # Iteration Recording Tests
    ==========================================================================#
    
    @testset "Constructor with All Arguments" begin
        seg1 = HistorySegment(coord_sys, name="solution_1")
        seg2 = HistorySegment(coord_sys, name="iteration_1")
        
        history = SpacecraftHistory([seg1], [seg2], false, true)
        
        @test length(history.segments) == 1
        @test length(history.iterations) == 1
        @test history.record_segments == false
        @test history.record_iterations == true
        @test history.segments[1].name == "solution_1"
        @test history.iterations[1].name == "iteration_1"
    end
    
    @testset "Default Flag Values" begin
        # Empty constructor
        history1 = SpacecraftHistory()
        @test history1.record_segments == true
        @test history1.record_iterations == false
        @test isempty(history1.segments)
        @test isempty(history1.iterations)
        
        # Constructor with segments vector
        seg = HistorySegment(coord_sys, name="test")
        history2 = SpacecraftHistory([seg])
        @test history2.record_segments == true
        @test history2.record_iterations == false
        @test length(history2.segments) == 1
        @test isempty(history2.iterations)
    end
    
    @testset "push_segment! Routing - Segments Only (Default)" begin
        history = SpacecraftHistory()
        seg1 = HistorySegment(coord_sys, name="seg1")
        seg2 = HistorySegment(coord_sys, name="seg2")
        
        @test history.record_segments == true
        @test history.record_iterations == false
        
        push_segment!(history, seg1)
        push_segment!(history, seg2)
        
        @test length(history.segments) == 2
        @test length(history.iterations) == 0
        @test history.segments[1].name == "seg1"
        @test history.segments[2].name == "seg2"
    end
    
    @testset "push_segment! Routing - Iterations Only" begin
        history = SpacecraftHistory()
        history.record_segments = false
        history.record_iterations = true
        
        seg1 = HistorySegment(coord_sys, name="iter1")
        seg2 = HistorySegment(coord_sys, name="iter2")
        
        push_segment!(history, seg1)
        push_segment!(history, seg2)
        
        @test length(history.segments) == 0
        @test length(history.iterations) == 2
        @test history.iterations[1].name == "iter1"
        @test history.iterations[2].name == "iter2"
    end
    
    @testset "push_segment! Routing - Both Flags True (Iterations Priority)" begin
        history = SpacecraftHistory()
        history.record_segments = true
        history.record_iterations = true
        
        seg = HistorySegment(coord_sys, name="test")
        push_segment!(history, seg)
        
        # Iterations take priority in if/elseif logic
        @test length(history.segments) == 0
        @test length(history.iterations) == 1
        @test history.iterations[1].name == "test"
    end
    
    @testset "push_segment! Routing - Both Flags False (No-op)" begin
        history = SpacecraftHistory()
        history.record_segments = false
        history.record_iterations = false
        
        seg = HistorySegment(coord_sys, name="test")
        push_segment!(history, seg)
        
        # Nothing recorded
        @test length(history.segments) == 0
        @test length(history.iterations) == 0
    end
    
    @testset "Flag Modification - Mutable Struct" begin
        history = SpacecraftHistory()
        
        # Start with default (segments enabled)
        @test history.record_segments == true
        @test history.record_iterations == false
        
        # Modify flags
        history.record_segments = false
        history.record_iterations = true
        
        @test history.record_segments == false
        @test history.record_iterations == true
        
        # Verify routing changes
        seg = HistorySegment(coord_sys, name="test")
        push_segment!(history, seg)
        
        @test length(history.segments) == 0
        @test length(history.iterations) == 1
    end
    
    @testset "Solver-like Workflow" begin
        history = SpacecraftHistory()
        
        # Phase 1: Collect iterations (optimization)
        history.record_segments = false
        history.record_iterations = true
        
        for i in 1:3
            seg = HistorySegment(coord_sys, name="iteration_$i")
            push_segment!(history, seg)
        end
        
        @test length(history.iterations) == 3
        @test length(history.segments) == 0
        
        # Phase 2: Collect final solution
        history.record_segments = true
        history.record_iterations = false
        
        seg_final = HistorySegment(coord_sys, name="solution")
        push_segment!(history, seg_final)
        
        @test length(history.iterations) == 3
        @test length(history.segments) == 1
        @test history.segments[1].name == "solution"
    end
    
    @testset "Base.show with Iterations" begin
        history = SpacecraftHistory()
        
        # Add solution segments
        seg1 = HistorySegment(coord_sys, name="solution_1")
        push_segment!(history, seg1)
        
        # Add iteration data
        history.record_segments = false
        history.record_iterations = true
        seg2 = HistorySegment(coord_sys, name="iter_1")
        seg3 = HistorySegment(coord_sys, name="iter_2")
        push_segment!(history, seg2)
        push_segment!(history, seg3)
        
        # Show should display iteration count
        output = sprint(show, history)
        @test occursin("1 segment", output)  # Solution segments
        @test occursin("Iterations: 2 segments", output)  # Iteration data
        @test occursin("diagnostic data", output)
    end
    
    @testset "Base.show without Iterations" begin
        history = SpacecraftHistory()
        seg = HistorySegment(coord_sys, name="test")
        push_segment!(history, seg)
        
        output = sprint(show, history)
        @test occursin("1 segment", output)
        @test !occursin("Iterations:", output)  # No mention of iterations
    end
    
    @testset "Base.show Iteration Pluralization" begin
        history = SpacecraftHistory()
        history.record_segments = false
        history.record_iterations = true
        
        # Single iteration
        seg1 = HistorySegment(coord_sys, name="iter1")
        push_segment!(history, seg1)
        output1 = sprint(show, history)
        @test occursin("Iterations: 1 segment", output1)
        
        # Multiple iterations
        seg2 = HistorySegment(coord_sys, name="iter2")
        push_segment!(history, seg2)
        output2 = sprint(show, history)
        @test occursin("Iterations: 2 segments", output2)
    end
    
    @testset "Base.copy with Iterations" begin
        history = SpacecraftHistory()
        
        # Add segments
        seg1 = HistorySegment(coord_sys, name="seg1")
        push_segment!(history, seg1)
        
        # Add iterations
        history.record_segments = false
        history.record_iterations = true
        seg2 = HistorySegment(coord_sys, name="iter1")
        push_segment!(history, seg2)
        
        # Copy
        history_copy = copy(history)
        
        # Verify independent copies
        @test length(history_copy.segments) == 1
        @test length(history_copy.iterations) == 1
        @test history_copy.record_segments == false
        @test history_copy.record_iterations == true
        
        # Verify independence - modify original
        seg3 = HistorySegment(coord_sys, name="seg2")
        history.record_segments = true
        history.record_iterations = false
        push_segment!(history, seg3)
        
        @test length(history.segments) == 2
        @test length(history_copy.segments) == 1  # Copy unchanged
        
        # Verify independence - modify copy
        seg4 = HistorySegment(coord_sys, name="iter2")
        history_copy.record_segments = false
        history_copy.record_iterations = true
        push_segment!(history_copy, seg4)
        
        @test length(history_copy.iterations) == 2
        @test length(history.iterations) == 1  # Original unchanged
    end
    
    @testset "Indexing and Iteration Operate on Segments Only" begin
        history = SpacecraftHistory()
        
        # Add solution segments
        for i in 1:3
            seg = HistorySegment(coord_sys, name="seg_$i")
            push_segment!(history, seg)
        end
        
        # Add iterations
        history.record_segments = false
        history.record_iterations = true
        for i in 1:5
            seg = HistorySegment(coord_sys, name="iter_$i")
            push_segment!(history, seg)
        end
        
        # Verify length only counts segments
        @test length(history) == 3
        @test length(history.iterations) == 5
        
        # Verify indexing only accesses segments
        @test history[1].name == "seg_1"
        @test history[2].name == "seg_2"
        @test history[3].name == "seg_3"
        
        # Verify range indexing
        range_result = history[1:2]
        @test length(range_result) == 2
        @test range_result[1].name == "seg_1"
        @test range_result[2].name == "seg_2"
        
        # Verify iteration
        names = [seg.name for seg in history]
        @test names == ["seg_1", "seg_2", "seg_3"]
        
        # Verify first/last
        @test first(history).name == "seg_1"
        @test last(history).name == "seg_3"
        @test firstindex(history) == 1
        @test lastindex(history) == 3
    end
    
    @testset "Direct Iteration Access" begin
        history = SpacecraftHistory()
        history.record_segments = false
        history.record_iterations = true
        
        # Add iterations
        for i in 1:4
            seg = HistorySegment(coord_sys, name="iter_$i")
            push_segment!(history, seg)
        end
        
        # Access iterations directly
        @test length(history.iterations) == 4
        @test history.iterations[1].name == "iter_1"
        @test history.iterations[4].name == "iter_4"
        
        # Verify empty when not recording iterations
        history2 = SpacecraftHistory()
        @test history2.record_iterations == false
        seg = HistorySegment(coord_sys, name="test")
        push_segment!(history2, seg)
        @test isempty(history2.iterations)
    end
    
    @testset "Segments and Iterations Independence" begin
        history = SpacecraftHistory()
        
        # Add to segments
        seg1 = HistorySegment(coord_sys, name="seg1")
        push_segment!(history, seg1)
        
        # Switch to iterations
        history.record_segments = false
        history.record_iterations = true
        seg2 = HistorySegment(coord_sys, name="iter1")
        push_segment!(history, seg2)
        
        # Switch back to segments
        history.record_segments = true
        history.record_iterations = false
        seg3 = HistorySegment(coord_sys, name="seg2")
        push_segment!(history, seg3)
        
        # Verify both vectors maintained independently
        @test length(history.segments) == 2
        @test length(history.iterations) == 1
        @test history.segments[1].name == "seg1"
        @test history.segments[2].name == "seg2"
        @test history.iterations[1].name == "iter1"
    end
    
end
