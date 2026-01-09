using Test
using AstroSolve
using AstroCallbacks
using AstroManeuvers
using AstroModels
using AstroProp
using AstroEpochs
using AstroStates
using AstroUniverse
using AstroFrames
using OrdinaryDiffEq

"""
    setup_simple_target()

Create a simple trajectory optimization problem based on Ex_SimpleTarget.
Returns a sequence and the spacecraft for testing history recording.
"""
function setup_simple_target()
    # Create spacecraft with default orbital state
    sat = Spacecraft()
    
    # Create the propagator with point mass gravity model
    gravity = PointMassGravity(earth, (moon, sun))
    forces = ForceModel(gravity)
    integ = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
    prop = OrbitPropagator(forces, integ)
    
    # Create an impulsive maneuver (only V component will vary)
    toi = ImpulsiveManeuver(axes=VNB(), element1=0.1)
    
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
    toi_event = Event(
        name = "TOI Maneuver",
        event = fun_toi,
        vars = [var_toi],
        funcs = []
    )
    
    # Create propagation event to apoapsis with position constraint
    fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
    prop_event = Event(
        name = "Propagate to Apoapsis",
        event = fun_prop_apo,
        funcs = [pos_con]
    )
    
    # Build sequence: maneuver first, then propagate to apoapsis
    seq = Sequence()
    add_events!(seq, prop_event, [toi_event])
    
    return seq, sat
end

@testset "Iteration Recording Tests" begin
    
    @testset "Default: No Iteration Recording" begin
        seq, sat = setup_simple_target()
        
        # Verify initial state
        @test sat.history.record_segments == true
        @test sat.history.record_iterations == false
        @test isempty(sat.history.iterations)
        
        # Solve with default settings (record_iterations=false)
        result = trajectory_solve(seq)
        
        # Verify no iterations recorded
        @test isempty(sat.history.iterations)
        
        # Verify solution recorded (default record_segments=true)
        @test length(sat.history.segments) > 0
        
        # Verify solution has data points
        total_points = sum(length(seg.times) for seg in sat.history.segments)
        @test total_points > 0
    end
    
    @testset "With Iteration Recording Enabled" begin
        seq, sat = setup_simple_target()
        
        # Verify spacecraft identity through solve
        println("Spacecraft object_id before solve: ", objectid(sat))
        println("Flags before solve: record_segments=", sat.history.record_segments, 
                ", record_iterations=", sat.history.record_iterations)
        
        # Solve with iteration recording enabled
        result = trajectory_solve(seq; record_iterations=true)
        
        println("Spacecraft object_id after solve: ", objectid(sat))
        println("Flags after solve: record_segments=", sat.history.record_segments,
                ", record_iterations=", sat.history.record_iterations)
        println("Iterations count: ", length(sat.history.iterations))
        println("Segments count: ", length(sat.history.segments))
        
        # Verify iterations were recorded
        @test length(sat.history.iterations) > 0
        
        # Verify solution also recorded (original flags restored)
        @test length(sat.history.segments) > 0
        
        # Verify iteration segments have data
        if !isempty(sat.history.iterations)
            iter_points = sum(length(seg.times) for seg in sat.history.iterations; init=0)
            @test iter_points > 0
        end
        
        # Verify segments and iterations are independent
        @test sat.history.iterations !== sat.history.segments
    end
    
    @testset "Flag Save and Restore" begin
        seq, sat = setup_simple_target()
        
        # Save original flag values
        original_seg_flag = sat.history.record_segments
        original_iter_flag = sat.history.record_iterations
        
        # Solve with iteration recording
        result = trajectory_solve(seq; record_iterations=true)
        
        # Verify flags were restored to original values
        @test sat.history.record_segments == original_seg_flag
        @test sat.history.record_iterations == original_iter_flag
    end
    
    @testset "Custom Initial Flags - Both True" begin
        seq, sat = setup_simple_target()
        
        # Set custom flags (both true)
        sat.history.record_segments = true
        sat.history.record_iterations = true
        
        # Solve
        result = trajectory_solve(seq; record_iterations=false)
        
        # Verify flags restored
        @test sat.history.record_segments == true
        @test sat.history.record_iterations == true
    end
    
    @testset "No Solution Recording When record_segments=false" begin
        seq, sat = setup_simple_target()
        
        # Disable solution recording before solve
        sat.history.record_segments = false
        sat.history.record_iterations = false
        
        # Solve with iteration recording enabled
        result = trajectory_solve(seq; record_iterations=true)
        
        # Verify iterations were captured during solve
        @test length(sat.history.iterations) > 0
        
        # Verify NO solution recorded (original record_segments was false)
        @test isempty(sat.history.segments)
        
        # Verify flags restored to original (false)
        @test sat.history.record_segments == false
        @test sat.history.record_iterations == false
    end
    
    @testset "Iteration Count Reasonable" begin
        seq, sat = setup_simple_target()
        
        # Solve with iteration recording
        result = trajectory_solve(seq; record_iterations=true)
        
        n_iterations = length(sat.history.iterations)
        
        # Verify iteration count is reasonable (not zero, not thousands)
        @test n_iterations > 0
        @test n_iterations < 1000  # Sanity check - simple problem shouldn't need thousands
    end
    
    @testset "Multiple Solves - History Accumulation" begin
        seq, sat = setup_simple_target()
        
        # First solve
        result1 = trajectory_solve(seq)
        n_segments_1 = length(sat.history.segments)
        
        # Second solve (history accumulates)
        result2 = trajectory_solve(seq)
        n_segments_2 = length(sat.history.segments)
        
        # Verify history accumulated from both solves
        @test n_segments_2 > n_segments_1
        @test n_segments_2 >= 2 * n_segments_1  # At least doubled (may be more due to propagation steps)
    end
    
    @testset "Iteration Data Independence from Solution" begin
        seq, sat = setup_simple_target()
        
        # Solve with iteration recording
        result = trajectory_solve(seq; record_iterations=true)
        
        # Verify both vectors populated
        @test length(sat.history.iterations) > 0
        @test length(sat.history.segments) > 0
        
        # Verify they're independent vectors
        @test sat.history.iterations !== sat.history.segments
        
        # Modify iterations vector shouldn't affect segments
        original_seg_count = length(sat.history.segments)
        push!(sat.history.iterations, sat.history.iterations[1])
        @test length(sat.history.segments) == original_seg_count
    end
    
    @testset "Convergence and Solution Quality" begin
        seq, sat = setup_simple_target()
        
        # Solve
        result = trajectory_solve(seq; record_iterations=true)
        
        # Verify convergence info returned
        @test haskey(result, :info)
        
        # Verify constraint satisfaction (should be close to target)
        pos_target = 55000.0
        @test abs(result.constraints[1] - pos_target) < 1.0  # Within 1 km
    end
    
end
