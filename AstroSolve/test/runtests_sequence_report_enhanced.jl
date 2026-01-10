using Test
using Epicycle


    
    # Helper function to capture printed output for validation
    function capture_output(func, args...)
        mktemp() do path, io
            redirect_stdout(io) do
                func(args...)
            end
            close(io)
            return read(path, String)
        end
    end
    
    @testset "Basic Functionality Test" begin
        # Create test spacecraft
        sat = Spacecraft(
            state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
            time = Time("2000-01-01T11:59:28.000", UTC(), ISOT())
        )
        
        # Create simple sequence
        toi = ImpulsiveManeuver(axes=VNB(), element1=1.5, element2=0.0, element3=0.0)
        
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "test_maneuver",
            lower_bound = [-2.0, -1.0, -1.0],
            upper_bound = [2.0, 1.0, 1.0],
            shift = [0.0, 0.0, 0.0],
            scale = [1.0, 1.0, 1.0]
        )
        
        constraint_test = Constraint(
            calc = OrbitCalc(sat, SMA()),
            lower_bounds = [42000.0],
            upper_bounds = [42000.0]
        )
        
        test_event = Event(
            name = "Test Event",
            event = () -> maneuver!(sat, toi),
            vars = [var_toi],
            funcs = [constraint_test]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test that functions don't crash
        @test_nowarn report_sequence(seq)
        
        # Test solution report
        mock_result = (
            variables = [1.2, 0.1, -0.05],
            objective = [42001.5],
            constraints = [42001.5],
            info = true
        )
        @test_nowarn report_solution(seq, mock_result)
    end
    
    @testset "Output Content Validation" begin
        # Create test sequence
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        toi = ImpulsiveManeuver(axes=VNB(), element1=1.0, element2=0.0, element3=0.0)
        
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "test_var",
            lower_bound = [-5.0, -1.0, -1.0],
            upper_bound = [5.0, 1.0, 1.0],
            shift = [0.0, 0.0, 0.0],
            scale = [1.0, 1.0, 1.0]
        )
        
        sma_constraint = Constraint(
            calc = OrbitCalc(sat, SMA()),
            lower_bounds = [42000.0],
            upper_bounds = [42000.0]
        )
        
        test_event = Event(
            name = "Test Maneuver",
            event = () -> maneuver!(sat, toi),
            vars = [var_toi],
            funcs = [sma_constraint]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test sequence report content
        output = capture_output(report_sequence, seq)
        @test occursin("TRAJECTORY SEQUENCE SUMMARY", output)
        @test occursin("Total Events: 1", output)
        @test occursin("Variable Objects: 1 (3 optimization variables)", output)
        @test occursin("Constraint Objects: 1 (1 constraint functions)", output)
        @test occursin("Event 1: \"Test Maneuver\"", output)
        @test occursin("test_var: DeltaVVector() (ManeuverCalc) (3 components)", output)
        @test occursin("SMA() (OrbitCalc) = 42000.0", output)
        
        # Test solution report content
        mock_result = (
            variables = [1.5, 0.2, -0.1],
            constraints = [42001.0],
            info = true
        )
        
        sol_output = capture_output(report_solution, seq, mock_result)
        @test occursin("TRAJECTORY SOLUTION REPORT", sol_output)
        @test occursin("Converged: true", sol_output)
        @test occursin("test_var (DeltaVVector() (ManeuverCalc)):", sol_output)
        @test occursin("Component 1: 1.5", sol_output)
        @test occursin("Component 2: 0.2", sol_output)
        @test occursin("Component 3: -0.1", sol_output)
        @test occursin("Total ΔV:", sol_output)
    end
    
    @testset "Multi-Event Sequence Test" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        # Create propagator
        gravity = PointMassGravity(earth, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        
        # Events
        prop_event = Event(
            name = "Propagation",
            event = () -> propagate!(prop, sat, StopAt(sat, PosZ(), 0.0)),
        )
        
        toi = ImpulsiveManeuver(axes=VNB(), element1=1.0, element2=0.0, element3=0.0)
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "toi_dv",
            lower_bound = [-3.0, -1.0, -1.0],
            upper_bound = [3.0, 1.0, 1.0],
            shift = [0.0, 0.0, 0.0],
            scale = [1.0, 1.0, 1.0]
        )
        
        toi_event = Event(
            name = "TOI Maneuver",
            event = () -> maneuver!(sat, toi),
            vars = [var_toi]
        )
        
        seq = Sequence()
        add_events!(seq, toi_event, [prop_event])
        
        output = capture_output(report_sequence, seq)
        @test occursin("Total Events: 2", output)
        @test occursin("Variable Objects: 1 (3 optimization variables)", output)
        @test occursin("Event 1: \"Propagation\"", output)
        @test occursin("Event 2: \"TOI Maneuver\"", output)
        @test occursin("Execution Order: [\"Propagation\" → \"TOI Maneuver\"]", output)
    end
    
    @testset "Constraint Constructor Tests" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        calc = OrbitCalc(sat, SMA())
        
        # Test equality constraint (both bounds)
        con1 = Constraint(calc=calc, lower_bounds=[42000.0], upper_bounds=[42000.0])
        @test con1.lower_bounds == [42000.0]
        @test con1.upper_bounds == [42000.0]
        
        # Test lower bound only (upper should default to Inf)
        con2 = Constraint(calc=calc, lower_bounds=[40000.0])
        @test con2.lower_bounds == [40000.0]
        @test con2.upper_bounds == [Inf]
        
        # Test upper bound only (lower should default to -Inf)
        con3 = Constraint(calc=calc, upper_bounds=[45000.0])
        @test con3.lower_bounds == [-Inf]
        @test con3.upper_bounds == [45000.0]
        
        # Test error when no bounds specified
        @test_throws ArgumentError Constraint(calc=calc)
    end
    
    @testset "Empty Sequence Test" begin
        empty_seq = Sequence()
        @test_nowarn report_sequence(empty_seq)

        output = capture_output(report_sequence, empty_seq)
        @test occursin("Total Events: 0", output)
        @test occursin("Variable Objects: 0 (0 optimization variables)", output)
        @test occursin("Constraint Objects: 0 (0 constraint functions)", output)
    end

    @testset "Large Sequence Truncation Test" begin
        # Create 5 events to trigger execution order truncation (>3 events)
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        # Create 5 simple events
        events = Event[]
        for i in 1:5
            event = Event(
                name = "Event_$i",
                event = () -> nothing  # Simple no-op function
            )
            push!(events, event)
        end
        
        seq = Sequence()
        add_events!(seq, events[1], Event[])
        for i in 2:5
            add_events!(seq, events[i], [events[i-1]])
        end
        
        output = capture_output(report_sequence, seq)
        @test occursin("Total Events: 5", output)
        # Should show truncated execution order: Event_1 → Event_2 → ... → Event_5
        @test occursin("Event_1\" → \"Event_2\" → ... → \"Event_5", output)
    end

    @testset "Single Component Variable Test" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        # Create a single-component variable (scalar orbital element)
        var_single = SolverVariable(
            calc = OrbitCalc(sat, SMA()),
            name = "sma_var",
            lower_bound = [40000.0],
            upper_bound = [45000.0],
            shift = [0.0],
            scale = [1.0]
        )
        
        test_event = Event(
            name = "SMA Test",
            event = () -> nothing,
            vars = [var_single]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        output = capture_output(report_sequence, seq)
        @test occursin("sma_var: SMA() (OrbitCalc) ∈ [40000.0, 45000.0]", output)
    end

    @testset "Inequality Constraints Test" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        calc = OrbitCalc(sat, SMA())
        
        # Test lower bound only (≥)
        con_lower = Constraint(calc=calc, lower_bounds=[40000.0])
        
        # Test upper bound only (≤)  
        con_upper = Constraint(calc=calc, upper_bounds=[50000.0])
        
        test_event = Event(
            name = "Inequality Test",
            event = () -> nothing,
            funcs = [con_lower, con_upper]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        output = capture_output(report_sequence, seq)
        @test occursin("≥ 40000.0", output)
        @test occursin("≤ 50000.0", output)
    end

    @testset "Multiple Stateful Objects Test" begin
        # Create multiple spacecraft to trigger stateful object pluralization
        sat1 = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        sat2 = Spacecraft(
            state = CartesianState([8000.0, 0.0, 0.0, 0.0, 7.0, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        sat3 = Spacecraft(
            state = CartesianState([9000.0, 0.0, 0.0, 0.0, 6.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        # Create multiple maneuvers to get more stateful objects
        man1 = ImpulsiveManeuver(axes=VNB(), element1=1.0)
        man2 = ImpulsiveManeuver(axes=VNB(), element1=1.5)
        
        var1 = SolverVariable(
            calc = ManeuverCalc(man1, sat1, DeltaVVector()),
            name = "var1",
            lower_bound = [-2.0, -1.0, -1.0],
            upper_bound = [2.0, 1.0, 1.0]
        )
        
        var2 = SolverVariable(
            calc = ManeuverCalc(man2, sat2, DeltaVVector()),
            name = "var2", 
            lower_bound = [-3.0, -1.0, -1.0],
            upper_bound = [3.0, 1.0, 1.0]
        )
        
        event1 = Event(
            name = "Multi Event 1",
            event = () -> maneuver!(sat1, man1),
            vars = [var1]
        )
        
        event2 = Event(
            name = "Multi Event 2", 
            event = () -> maneuver!(sat2, man2),
            vars = [var2]
        )
        
        seq = Sequence()
        add_events!(seq, event1, Event[])
        add_events!(seq, event2, [event1])
        
        output = capture_output(report_sequence, seq)
        # Should show multiple stateful objects with counts
        @test occursin("STATEFUL OBJECTS:", output)
        # Should have multiple spacecraft and maneuvers, some with counts >1
        @test (occursin("×", output) || (occursin("Spacecraft", output) && occursin("ImpulsiveManeuver", output)))
    end

    @testset "Fixed Variable Test" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        # Create variable with fixed bounds (lower == upper)
        var_fixed = SolverVariable(
            calc = OrbitCalc(sat, SMA()),
            name = "fixed_var",
            lower_bound = [42000.0],
            upper_bound = [42000.0],  # Same as lower
            shift = [0.0],
            scale = [1.0]
        )
        
        test_event = Event(
            name = "Fixed Test",
            event = () -> nothing,
            vars = [var_fixed]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test in solution report 
        mock_result = (
            variables = [42000.0],
            constraints = Float64[],
            info = true
        )
        
        sol_output = capture_output(report_solution, seq, mock_result)
        @test occursin("(fixed at 42000.0)", sol_output)
    end

    @testset "Backwards Compatibility Test" begin
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2020-01-01T00:00:00", UTC(), ISOT())
        )
        
        var_test = SolverVariable(
            calc = OrbitCalc(sat, SMA()),
            name = "test_var",
            lower_bound = [42000.0],
            upper_bound = [42000.0]
        )
        
        test_event = Event(
            name = "Compat Test",
            event = () -> nothing,
            vars = [var_test]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test result without .constraints field (old format)
        old_result = (
            variables = [42001.0],
            objective = [42001.0],  # No .constraints field
            info = false
        )
        
        @test_nowarn report_solution(seq, old_result)

        sol_output = capture_output(report_solution, seq, old_result)
        @test occursin("Converged: false", sol_output)
    end

    @testset "Multi-Component Constraint Test" begin
        # Test multi-component constraints using PositionVector (3 components)
        sat = Spacecraft()
        
        # Create a 3-component position vector constraint 
        pos_con = Constraint(
            calc = OrbitCalc(sat, PositionVector()),
            lower_bounds = [-55000.0, 0.5, 0.3],
            upper_bounds = [-55000.0, 0.5, 0.3],
            scale = [1.0, 1.0, 1.0]
        )
        
        test_event = Event(
            name = "Multi-Component Test",
            event = () -> nothing,
            funcs = [pos_con]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test sequence report shows multi-component constraint
        output = capture_output(report_sequence, seq)
        @test occursin("PositionVector() (OrbitCalc) (3 components)", output)
        
        # Test solution report with multi-component constraint values
        mock_result = (
            variables = Float64[],
            constraints = [-55000.1, 0.51, 0.31],  # 3 constraint values
            info = true
        )
        
        sol_output = capture_output(report_solution, seq, mock_result)
        @test occursin("PositionVector() (OrbitCalc) (3 components):", sol_output)
        @test occursin("Component 1: -55000.1", sol_output)
        @test occursin("Component 2: 0.51", sol_output) 
        @test occursin("Component 3: 0.31", sol_output)
        @test occursin("target: -55000.0", sol_output)
    end

nothing