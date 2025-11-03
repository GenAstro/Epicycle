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
            event = () -> maneuver(sat, toi),
            vars = [var_toi],
            funcs = [constraint_test]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test that functions don't crash
        @test_nowarn sequence_report(seq)
        
        # Test solution report
        mock_result = (
            variables = [1.2, 0.1, -0.05],
            objective = [42001.5],
            constraints = [42001.5],
            info = true
        )
        @test_nowarn solution_report(seq, mock_result)
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
            event = () -> maneuver(sat, toi),
            vars = [var_toi],
            funcs = [sma_constraint]
        )
        
        seq = Sequence()
        add_events!(seq, test_event, Event[])
        
        # Test sequence report content
        output = capture_output(sequence_report, seq)
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
        
        sol_output = capture_output(solution_report, seq, mock_result)
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
            event = () -> propagate(prop, sat, StopAt(sat, PosZ(), 0.0)),
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
            event = () -> maneuver(sat, toi),
            vars = [var_toi]
        )
        
        seq = Sequence()
        add_events!(seq, toi_event, [prop_event])
        
        output = capture_output(sequence_report, seq)
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
        @test_nowarn sequence_report(empty_seq)
        
        output = capture_output(sequence_report, empty_seq)
        @test occursin("Total Events: 0", output)
        @test occursin("Variable Objects: 0 (0 optimization variables)", output)
        @test occursin("Constraint Objects: 0 (0 constraint functions)", output)
    end

nothing