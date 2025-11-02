using Test
using AstroSolve
using AstroFun
using AstroStates
using AstroEpochs
using AstroMan

@testset "Show Methods Tests" begin

    @testset "Event show method" begin
        
        @testset "Empty Event" begin
            # Test basic empty event
            event = Event()
            output = string(event)
            @test contains(output, "Event(<unnamed>; 0 vars, 0 funcs)")
        end
        
        @testset "Named Event without vars/funcs" begin
            event = Event(name = "test_event")
            output = string(event)
            @test contains(output, "Event(\"test_event\"; 0 vars, 0 funcs)")
        end
        
        @testset "Event with SolverVariable (DeltaVVector)" begin
            # Create a test spacecraft and maneuver
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
            calc = ManeuverCalc(man, sc, DeltaVVector())
            sv = SolverVariable(calc = calc, name = "delta_v", 
                               lower_bound = [-1.0, -1.0, -1.0], 
                               upper_bound = [1.0, 1.0, 1.0])
            
            # Event with one DeltaVVector variable (should show 3 vars)
            event = Event(name = "maneuver", vars = [sv])
            output = string(event)
            @test contains(output, "Event(\"maneuver\"; 3 vars, 0 funcs)")
        end
        
        @testset "Event with multiple SolverVariables" begin
            # Create test spacecraft
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            # Create DeltaVVector variable (3 scalar vars)
            man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
            calc_dv = ManeuverCalc(man, sc, DeltaVVector())
            sv_dv = SolverVariable(calc = calc_dv, name = "delta_v")
            
            # Create scalar variable (1 scalar var)
            calc_ecc = OrbitCalc(sc, Ecc())
            sv_ecc = SolverVariable(calc = calc_ecc, name = "eccentricity")
            
            # Event with both variables (should show 4 total vars)
            event = Event(name = "combined", vars = [sv_dv, sv_ecc])
            output = string(event)
            @test contains(output, "Event(\"combined\"; 4 vars, 0 funcs)")
        end
        
        @testset "Event with functions" begin
            # Create a mock constraint (this is simplified - real constraints would be more complex)
            # For testing purposes, we'll create an event with functions
            test_func() = 0.0  # Simple mock function
            
            event = Event(name = "with_func", funcs = [test_func])
            output = string(event)
            # Should count the function as 1 func
            @test contains(output, "Event(\"with_func\"; 0 vars, 1 funcs)")
        end
        
    end
    
    @testset "SolverVariable show method" begin
        
        @testset "DeltaVVector SolverVariable" begin
            # Create test spacecraft
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
            calc = ManeuverCalc(man, sc, DeltaVVector())
            sv = SolverVariable(
                calc = calc, 
                name = "test_delta_v",
                lower_bound = [-2.0, -1.0, -0.5],
                upper_bound = [2.0, 1.0, 0.5],
                shift = [0.1, 0.2, 0.3],
                scale = [1.5, 2.0, 2.5]
            )
            
            # Capture output
            io = IOBuffer()
            show(io, sv)
            output = String(take!(io))
            
            # Test key components of output
            @test contains(output, "SolverVariable(ManeuverCalc")
            @test contains(output, "numvars:     3")
            @test contains(output, "lower_bound: [-2.0, -1.0, -0.5]")
            @test contains(output, "upper_bound: [2.0, 1.0, 0.5]")
            @test contains(output, "shift:       [0.1, 0.2, 0.3]")
            @test contains(output, "scale:       [1.5, 2.0, 2.5]")
            @test contains(output, "name:        test_delta_v")
        end
        
        @testset "Scalar SolverVariable (Eccentricity)" begin
            # Create test spacecraft
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            calc = OrbitCalc(sc, Ecc())
            sv = SolverVariable(
                calc = calc,
                name = "orbital_eccentricity", 
                lower_bound = 0.0,
                upper_bound = 0.99,
                shift = 0.01,
                scale = 10.0
            )
            
            # Capture output
            io = IOBuffer()
            show(io, sv)
            output = String(take!(io))
            
            # Test key components
            @test contains(output, "SolverVariable(OrbitCalc")
            @test contains(output, "numvars:     1")
            @test contains(output, "lower_bound: [0.0]")
            @test contains(output, "upper_bound: [0.99]")
            @test contains(output, "shift:       [0.01]")
            @test contains(output, "scale:       [10.0]")
            @test contains(output, "name:        orbital_eccentricity")
        end
        
        @testset "SolverVariable with default values" begin
            # Create minimal SolverVariable to test defaults
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            calc = OrbitCalc(sc, Ecc())
            sv = SolverVariable(calc = calc)  # Using defaults
            
            # Capture output
            io = IOBuffer()
            show(io, sv)
            output = String(take!(io))
            
            # Test that output includes the expected structure
            @test contains(output, "SolverVariable(OrbitCalc")
            @test contains(output, "numvars:")
            @test contains(output, "lower_bound:")
            @test contains(output, "upper_bound:")
            @test contains(output, "shift:")
            @test contains(output, "scale:")
            @test contains(output, "name:")
        end
        
    end
    
    @testset "Edge cases and robustness" begin
        
        @testset "Event with empty name" begin
            event = Event(name = "")
            output = string(event)
            @test contains(output, "Event(<unnamed>")
        end
        
        @testset "Event display consistency" begin
            # This tests the specific fix for variable counting
            sc = Spacecraft(
                state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            )
            
            man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
            calc = ManeuverCalc(man, sc, DeltaVVector())
            sv = SolverVariable(calc = calc, name = "test_dv")
            
            event = Event(name = "consistency_test", vars = [sv])
            output = string(event)
            
            # Should show 3 vars (the actual scalar variables) not 1 var (the struct count)
            @test contains(output, "3 vars")
            @test !contains(output, "1 vars")  # This was the bug
        end
        
    end
    
end
nothing