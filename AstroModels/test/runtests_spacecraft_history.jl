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
    # Immutability Tests
    ==========================================================================#
    
    @testset "Struct Immutability" begin
        history = SpacecraftHistory()
        
        # Cannot reassign field
        @test_throws ErrorException history.segments = Vector{HistorySegment}()
    end
    
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
    
end
