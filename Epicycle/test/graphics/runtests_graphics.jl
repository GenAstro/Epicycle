using Test
using Epicycle
using EpicycleBase
using AstroModels
using AstroStates
using AstroEpochs
using AstroUniverse
using LinearAlgebra

@testset "Graphics Module Tests" begin
    
    # Helper: Create a test spacecraft with history
    function create_test_spacecraft(name::String="TestSat")
        coord_sys = CoordinateSystem(earth, ICRFAxes())

        # LEO orbit
        r = [7000.0, 0.0, 0.0]  # km
        v = [0.0, 7.5, 0.0]     # km/s
        state = CartesianState(r, v)
        epoch = Time("2024-01-01T00:00:00", TAI(), ISOT())
        sc = Spacecraft(state=state, time=epoch, mass=1000.0, name=name, coord_sys=coord_sys)
        
        # Create a history segment and add trajectory points
        segment = HistorySegment(coord_sys, name="test_orbit")
        for i in 1:10
            θ = 2π * i / 10
            r_hist = [7000*cos(θ), 7000*sin(θ), 0.0]
            v_hist = [-7.5*sin(θ), 7.5*cos(θ), 0.0]
            push!(segment.times, epoch + i*600.0)
            push!(segment.states, CartesianState(r_hist, v_hist))
        end
        push_segment!(sc.history, segment)
        
        return sc
    end
    
    @testset "View3D Construction" begin
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        
        # Basic construction
        view = View3D(coord_sys=coord_sys)
        @test view.coord_sys == coord_sys
        @test isempty(view.spacecraft)
        @test isempty(view.options)
        @test view._scene === nothing
        
        # Type stability
        @test view isa View3D
        @test view.spacecraft isa Vector{Spacecraft}
        @test view.options isa Dict
    end
    
    @testset "View3D show methods" begin
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        view = View3D(coord_sys=coord_sys)
        
        io = IOBuffer()
        show(io, view)
        output = String(take!(io))
        @test contains(output, "View3D")
        @test contains(output, "Earth")
        @test contains(output, "ICRFAxes")
        @test contains(output, "0 spacecraft")
        
        # After adding spacecraft
        sc = create_test_spacecraft()
        add_spacecraft!(view, sc)
        show(io, view)
        output = String(take!(io))
        @test contains(output, "1 spacecraft")
    end
    
    @testset "add_spacecraft! validation" begin
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        view = View3D(coord_sys=coord_sys)
        
        # Valid spacecraft with history - returns view for chaining
        sc = create_test_spacecraft("ValidSat")
        result = add_spacecraft!(view, sc)
        @test result === view
        @test length(view.spacecraft) == 1
        @test view.spacecraft[1] === sc
        
        # Default options stored correctly
        @test haskey(view.options, sc)
        @test view.options[sc][:show_iterations] == false
        
        # Custom options stored correctly
        sc2 = create_test_spacecraft("IterSat")
        add_spacecraft!(view, sc2; show_iterations=true)
        @test view.options[sc2][:show_iterations] == true
        
        # Empty history throws error
        sc_empty = Spacecraft(
            state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]), 
            time=Time("2024-01-01T00:00:00", TAI(), ISOT()), 
            mass=1000.0, 
            name="Empty", 
            coord_sys=coord_sys
        )
        @test_throws ErrorException add_spacecraft!(view, sc_empty)
        err = try
            add_spacecraft!(view, sc_empty)
        catch e
            e
        end
        @test occursin("history is empty", err.msg)
    end
    
    @testset "display_view validation" begin
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        view = View3D(coord_sys=coord_sys)
        
        # Empty view throws error with clear message
        @test_throws ErrorException display_view(view)
        err = try
            display_view(view)
        catch e
            e
        end
        @test occursin("No spacecraft added", err.msg)
        
        # Coordinate system mismatch throws error (checks first spacecraft only)
        other_coord = CoordinateSystem(earth, MJ2000Axes())
        sc_mismatch = Spacecraft(
            state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time=Time("2024-01-01T00:00:00", TAI(), ISOT()),
            mass=1000.0,
            name="Mismatch",
            coord_sys=other_coord
        )
        # Add history so it passes add_spacecraft! validation
        segment = HistorySegment(other_coord, name="test")
        push!(segment.times, Time("2024-01-01T00:00:00", TAI(), ISOT()))
        push!(segment.states, CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]))
        push_segment!(sc_mismatch.history, segment)
        add_spacecraft!(view, sc_mismatch)
        
        @test_throws ErrorException display_view(view)
        err = try
            display_view(view)
        catch e
            e
        end
        @test occursin("coordinate system", lowercase(err.msg))
    end
    
    @testset "Multiple spacecraft" begin
        coord_sys = CoordinateSystem(earth, ICRFAxes())
        view = View3D(coord_sys=coord_sys)
        
        sc1 = create_test_spacecraft("Sat1")
        sc2 = create_test_spacecraft("Sat2")
        sc3 = create_test_spacecraft("Sat3")
        
        # Adding preserves order
        add_spacecraft!(view, sc1)
        add_spacecraft!(view, sc2; show_iterations=true)
        add_spacecraft!(view, sc3; show_iterations=false)
        
        @test length(view.spacecraft) == 3
        @test view.spacecraft[1] === sc1
        @test view.spacecraft[2] === sc2
        @test view.spacecraft[3] === sc3
        
        # Each spacecraft has independent options
        @test view.options[sc1][:show_iterations] == false
        @test view.options[sc2][:show_iterations] == true
        @test view.options[sc3][:show_iterations] == false
        # Visual test moved to visual_regression.jl
    end
    
    @testset "CAD model loading failure" begin
        # Test error handling when CAD model file doesn't exist
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "BadModel",
            cad_model = CADModel(
                file_path = "nonexistent_file.obj",
                scale = 1.0,
                visible = true
            )
        )
        
        # Propagate so spacecraft has history
        gravity = PointMassGravity(earth, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 100.0))
        
        # Should warn but not error - falls back to marker
        view = View3D()
        add_spacecraft!(view, sat)
        
        @test_logs (:warn, r"Failed to load spacecraft model") display_view(view; size=(400, 300))
    end
    
    @testset "show_iterations warning when no iterations" begin
        # Test warning when user requests iterations but none exist
        sat = create_test_spacecraft("NoIters")
        
        view = View3D()
        add_spacecraft!(view, sat; show_iterations=true)
        
        @test_logs (:warn, r"show_iterations=true.*does not contain iterations") display_view(view; size=(400, 300))
    end
    
    @testset "Empty segment handling" begin
        # Test that empty segments are skipped during rendering
        sat = create_test_spacecraft("EmptySeg")
        
        # Insert an empty segment in the middle (not at the end)
        # First, save the existing segment
        first_seg = sat.history[1]
        
        # Clear history and re-add segments with empty one in middle
        empty!(sat.history.segments)
        push_segment!(sat.history, first_seg)
        
        # Add empty segment
        empty_seg = HistorySegment(sat.coord_sys, name="empty")
        push_segment!(sat.history, empty_seg)
        
        # Add another segment with data so spacecraft has valid final position
        final_seg = HistorySegment(sat.coord_sys, name="final")
        push!(final_seg.times, Time("2024-01-02T00:00:00", TAI(), ISOT()))
        push!(final_seg.states, CartesianState([7500.0, 0.0, 0.0, 0.0, 7.3, 0.0]))
        push_segment!(sat.history, final_seg)
        
        view = View3D()
        add_spacecraft!(view, sat)
        
        # Should render without error, skipping the empty segment
        @test display_view(view; size=(400, 300)) === nothing
    end
    
    # Body texture tests moved to visual_regression.jl
end
