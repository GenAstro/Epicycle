
using Test
using AstroSolve

# Dummy Event type for testing


e1 = Event()
e2 = Event()
e3 = Event()

seq = Sequence()

# e2 and e3 depend on e1
add_events!(seq, e2, [e1])
add_events!(seq, e3, [e1, e2])

@testset "adj_map structure" begin
    # e1 should have e2 and e3 as dependents
    @test Set(seq.adj_map[e1]) == Set([e2, e3])
    # e2 should have e3 as dependent
    @test seq.adj_map[e2] == [e3]
    # e3 should have no dependents
    @test seq.adj_map[e3] == []
end

using Test
using AstroSolve

# Create all nodes (events) as in the diagram
S3  = Event()
S4  = Event()
S5  = Event()
S6  = Event()
S7  = Event()
S8  = Event()
S9  = Event()
S10 = Event()
S11 = Event()
S12 = Event()
S13 = Event()
S14 = Event()
S15 = Event()
S16 = Event()
S40 = Event()
S41 = Event()
S42 = Event()
S45 = Event()
S46 = Event()
S47 = Event()
S56 = Event()
S57 = Event()
S58 = Event()
S59 = Event()
S60 = Event()
S61 = Event()
S62 = Event()
S79 = Event()
S80 = Event()

seq = Sequence()

# Add dependencies as in the diagram (S4 depends on S3, S40, S45, etc.)
add_events!(seq, S4, [S3])
add_events!(seq, S40, [S3])
add_events!(seq, S45, [S3])
add_events!(seq, S46, [S3])
add_events!(seq, S47, [S3])
add_events!(seq, S56, [S3])

add_events!(seq, S41, [S40])
add_events!(seq, S5, [S4])

add_events!(seq, S6, [S5])
add_events!(seq, S57, [S5])
add_events!(seq, S58, [S5])

add_events!(seq, S7, [S6])
add_events!(seq, S8, [S7])
add_events!(seq, S42, [S8])
add_events!(seq, S9, [S8])
add_events!(seq, S10, [S9])
add_events!(seq, S11, [S10])

add_events!(seq, S12, [S11])
add_events!(seq, S59, [S11])

add_events!(seq, S13, [S12])
add_events!(seq, S79, [S12])
add_events!(seq, S60, [S12])
add_events!(seq, S61, [S12])
add_events!(seq, S62, [S12])

add_events!(seq, S80, [S79])

add_events!(seq, S14, [S13])
add_events!(seq, S15, [S14])
add_events!(seq, S42, [S15])
add_events!(seq, S16, [S15])

@testset "adj_map structure (complex graph)" begin
    # S3 should have S4, S40, S45, S46, S47, S56 as dependents
    @test Set(seq.adj_map[S3]) == Set([S4, S40, S45, S46, S47, S56])
    # S4 should have S5 as dependent
    @test seq.adj_map[S4] == [S5]
    # S5 should have S6, S57, S58 as dependents
    @test Set(seq.adj_map[S5]) == Set([S6, S57, S58])
    # S6 should have S7 as dependent
    @test seq.adj_map[S6] == [S7]
    # S7 should have S8 as dependent
    @test seq.adj_map[S7] == [S8]
    # S8 should have S42 and S9 as dependents
    @test Set(seq.adj_map[S8]) == Set([S42, S9])
    # S9 should have S10 as dependent
    @test seq.adj_map[S9] == [S10]
    # S10 should have S11 as dependent
    @test seq.adj_map[S10] == [S11]
    # S11 should have S12 and S59 as dependents
    @test Set(seq.adj_map[S11]) == Set([S12, S59])
    # S12 should have S13, S79, S60, S61, S62 as dependents
    @test Set(seq.adj_map[S12]) == Set([S13, S79, S60, S61, S62])
    # S13 should have S14 as dependent
    @test seq.adj_map[S13] == [S14]
    # S14 should have S15 as dependent
    @test seq.adj_map[S14] == [S15]
    # S15 should have S42 and S16 as dependents
    @test Set(seq.adj_map[S15]) == Set([S42, S16])
    # S79 should have S80 as dependent
    @test seq.adj_map[S79] == [S80]
end

@testset "topological sort (complex graph)" begin
    sorted = topo_sort(seq)
    # For every event, all its dependencies must appear before it in the sorted list
    idx = Dict(ev => i for (i, ev) in enumerate(sorted))
    for (event, deps) in seq.adj_map
        for dep in deps
            @test idx[event] < idx[dep]  # event must come before its dependents
        end
    end
    # All events should be present in the sorted list
    all_events = Set(keys(seq.adj_map))
    for deps in values(seq.adj_map)
        for dep in deps
            push!(all_events, dep)
        end
    end
    @test Set(sorted) == all_events
end

sc_ctx = Spacecraft()

# Create Events with SolverVariables to test ordering unique variables
sv3  = SolverVariable(calc=OrbitCalc(sc_ctx, TA()),                                             name="sv3")
sv4  = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv4")
sv5  = SolverVariable(calc=OrbitCalc(sc_ctx, RAAN()),                                           name="sv5")
sv6  = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv6")
sv7  = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv7")
sv8  = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv8")
sv9  = SolverVariable(calc=OrbitCalc(sc_ctx, TA()),                                         name="sv9")

sv13 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv13")
sv14 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv14")
sv15 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv15")
sv16 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv16")

sv40 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv40")
sv41 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv41")
sv45 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv45")
sv46 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv46")
sv47 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv47")

sv48 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv48")
sv49 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv49")
sv50 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv50")

sv56 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv56")
sv57 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv57")
sv58 = SolverVariable(calc=ManeuverCalc(ImpulsiveManeuver(), sc_ctx, DeltaVVector()),           name="sv58")

#=
# Create Events with SolverVariables to test ordering unqique variables
sv3  = SolverVariable(obj=Spacecraft(), variable=PosVel(), name="sv3")
sv4  = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv4")
sv5  = SolverVariable(obj=Spacecraft(), variable=PosVel(), name="sv5")
sv6  = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv6")
sv7  = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv7")
sv8  = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv8")
sv9  = SolverVariable(obj=Spacecraft(), variable=PosVel(), name="sv9") 


sv13 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv13")
sv14 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv14")
sv15 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv15")
sv16 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv16")


sv40 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv40")
sv41 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv41")
sv45 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv45")
sv46 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv46")
sv47 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv47")

sv48 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv48")
sv49 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv49")
sv50 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv50")

sv56 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv56")
sv57 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv57")
sv58 = SolverVariable(obj=ImpulsiveManeuver(), variable=DeltaV(), name="sv58")
=#
# Create Events with names and vars
S3  = Event(name="S3",  vars=[sv3])
S4  = Event(name="S4",  vars=[sv4, sv5])
S5  = Event(name="S5",  vars=[sv5, sv57])
S6  = Event(name="S6",  vars=[sv6])
S7  = Event(name="S7",  vars=[sv7])
S8  = Event(name="S8",  vars=[sv8, sv9])
S13 = Event(name="S13", vars=[sv13,sv14,sv15,sv16])
S40 = Event(name="S40", vars=[sv40])
S41 = Event(name="S41", vars=[sv41])
S45 = Event(name="S45", vars=[sv45])
S46 = Event(name="S46", vars=[sv46])
S47 = Event(name="S47", vars=[sv46,sv47, sv48, sv49, sv50])
S56 = Event(name="S56", vars=[sv56])
S57 = Event(name="S57", vars=[sv57])
S58 = Event(name="S58", vars=[sv58])

@testset "order_unique_vars deduplication and order" begin
    # Simulate a sorted event list (choose any order you want to test)
    sorted_events = [S3, S4, S5, S6, S7, S8, S13, S40, S41, S45, S46, S47, S56, S57, S58]

    ordered_vars = order_unique_vars(sorted_events)

    # Build the expected order by hand, based on first appearance in sorted_events
    expected = [
        sv3,   # from S3
        sv4,   # from S4
        sv5,   # from S4 (appears again in S5, but should not be duplicated)
        sv57,  # from S5
        sv6,   # from S6
        sv7,   # from S7
        sv8,   # from S8
        sv9,   # from S8
        sv13,  # from S13
        sv14,  # from S13
        sv15,  # from S13
        sv16,  # from S13
        sv40,  # from S40
        sv41,  # from S41
        sv45,  # from S45
        sv46,  # from S46
        sv47,  # from S47
        sv48,  # from S47
        sv49,  # from S47
        sv50,  # from S47
        sv56,  # from S56
        sv58   # from S58
    ]

    @test ordered_vars == expected
end

@testset "Trajectory Solve Interface Tests" begin
    
    @testset "Default Options Test" begin
        # Use proven working configuration from Ex_SimpleTarget.jl
        sat = Spacecraft()
        
        # Create the propagator with point mass gravity model
        gravity = PointMassGravity(earth, (moon, sun))
        forces  = ForceModel(gravity)
        integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
        prop    = OrbitPropagator(forces, integ)
        
        # Create an impulsive maneuver (only V component will vary)
        toi = ImpulsiveManeuver(
            axes = VNB(),
            element1 = 0.1,
        )
        
        # Define DeltaVVector of toi as a solver variable
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "toi",
            lower_bound = [-10.0, 0.0, 0.0],
            upper_bound = [10.0, 0.0, 0.0],
        )
        
        # Define a constraint on position magnitude of spacecraft
        pos_target = 55000.0
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0],
        )
        
        # Create events
        fun_toi() = maneuver(sat, toi) 
        toi_event = Event(name = "TOI Maneuver", 
                          event = fun_toi,
                          vars = [var_toi],
                          funcs = [])
        
        fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
        prop_event = Event(name = "Propagate to Apoapsis", 
                           event = fun_prop_apo,
                           funcs = [pos_con])
        
        # Build sequence
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event]) 
        
        # Test default trajectory_solve(seq) - call only once
        result = nothing
        @test_nowarn result = trajectory_solve(seq)
        
        # Verify return structure
        @test haskey(result, :variables)
        @test haskey(result, :objective) 
        @test haskey(result, :constraints)
        @test haskey(result, :info)
        
        # Verify types and sizes
        @test isa(result.variables, Vector{Float64})
        @test isa(result.constraints, Vector{Float64})
        @test result.info isa Union{Bool, Symbol}
        @test length(result.variables) == 3  # 3-component DeltaV
        @test length(result.constraints) == 1  # 1 position constraint
    end
    
    @testset "Custom Options Test" begin
        # Use same proven working configuration with custom options
        sat = Spacecraft()
        
        gravity = PointMassGravity(earth, (moon, sun))
        forces  = ForceModel(gravity)
        integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
        prop    = OrbitPropagator(forces, integ)
        
        toi = ImpulsiveManeuver(
            axes = VNB(),
            element1 = 0.1,
        )
        
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "toi",
            lower_bound = [-10.0, 0.0, 0.0],
            upper_bound = [10.0, 0.0, 0.0],
        )
        
        pos_target = 55000.0
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0],
        )
        
        fun_toi() = maneuver(sat, toi) 
        toi_event = Event(name = "TOI Maneuver", 
                          event = fun_toi,
                          vars = [var_toi],
                          funcs = [])
        
        fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
        prop_event = Event(name = "Propagate to Apoapsis", 
                           event = fun_prop_apo,
                           funcs = [pos_con])
        
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event]) 
        
        # Create custom SNOW options - use same tolerance as default for reliability
        custom_ip_options = Dict(
            "max_iter" => 500,
            "tol" => 1e-6,  # Same as default, not tighter
            "file_print_level" => 0,
            "output_file" => "ipopt_custom_$(time_ns())_$(rand(UInt32)).out"
        )
        custom_options = Options(derivatives=ForwardFD(), solver=IPOPT(custom_ip_options))
        
        # Test custom trajectory_solve(seq, options) - call only once
        result = nothing
        @test_nowarn result = trajectory_solve(seq, custom_options)
        
        # Verify same return structure
        @test haskey(result, :variables)
        @test haskey(result, :objective)
        @test haskey(result, :constraints) 
        @test haskey(result, :info)
        @test length(result.variables) == 3
        @test length(result.constraints) == 1
    end
    
    @testset "Default Options Function Test" begin
        # Test AstroSolve.default_snow_options() directly (fully qualified)
        @test_nowarn options = AstroSolve.default_snow_options()
        
        options = AstroSolve.default_snow_options()
        
        # Verify it returns SNOW.Options
        @test isa(options, Options)
        
        # Verify it has expected components (basic structure test)
        @test isa(options.derivatives, ForwardFD)
        @test isa(options.solver, IPOPT)
        
        # Test that default options work with trajectory_solve using simple working config
        sat = Spacecraft()
        
        gravity = PointMassGravity(earth, (moon, sun))
        forces  = ForceModel(gravity)
        integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
        prop    = OrbitPropagator(forces, integ)
        
        toi = ImpulsiveManeuver(
            axes = VNB(),
            element1 = 0.05,  # Smaller maneuver for faster convergence
        )
        
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "default_test",
            lower_bound = [-5.0, 0.0, 0.0],
            upper_bound = [5.0, 0.0, 0.0],
        )
        
        pos_target = 50000.0  # Smaller target change
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0],
        )
        
        fun_toi() = maneuver(sat, toi) 
        toi_event = Event(name = "TOI Maneuver", 
                          event = fun_toi,
                          vars = [var_toi],
                          funcs = [])
        
        fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
        prop_event = Event(name = "Propagate to Apoapsis", 
                           event = fun_prop_apo,
                           funcs = [pos_con])
        
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event]) 
        
        # Use default options explicitly - call only once
        result = nothing
        @test_nowarn result = trajectory_solve(seq, AstroSolve.default_snow_options())
        
        @test result.info isa Union{Bool, Symbol}
    end
    
    @testset "Constraint Values Test" begin
        # Test constraint values using proven working config with multiple constraints
        sat = Spacecraft()
        
        gravity = PointMassGravity(earth, (moon, sun))
        forces  = ForceModel(gravity)
        integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
        prop    = OrbitPropagator(forces, integ)
        
        toi = ImpulsiveManeuver(
            axes = VNB(),
            element1 = 0.1,
        )
        
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "constraint_test",
            lower_bound = [-10.0, 0.0, 0.0],
            upper_bound = [10.0, 0.0, 0.0],
        )
        
        # Multiple constraints to test constraint vector
        pos_target = 55000.0
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0],
        )
        
        # Add a second constraint on SMA at apoapsis
        sma_con = Constraint(
            calc = OrbitCalc(sat, SMA()),
            lower_bounds = [26000.0],  # Reasonable SMA range
            upper_bounds = [30000.0],
            scale = [1.0],
        )
        
        fun_toi() = maneuver(sat, toi) 
        toi_event = Event(name = "TOI Maneuver", 
                          event = fun_toi,
                          vars = [var_toi],
                          funcs = [])
        
        fun_prop_apo() = propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
        prop_event = Event(name = "Propagate to Apoapsis", 
                           event = fun_prop_apo,
                           funcs = [pos_con, sma_con])
        
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event]) 
        
        result = trajectory_solve(seq)
        
        # Verify constraint values are returned
        @test length(result.constraints) == 2  # 2 constraints
        @test isa(result.constraints[1], Float64)
        @test isa(result.constraints[2], Float64)
        
        # Constraint values should be within reasonable bounds for this problem
        # (Note: exact values depend on convergence, but should be in expected ranges)
        @test result.constraints[1] > 20000.0  # Position magnitude reasonable
        @test result.constraints[2] > 20000.0  # SMA reasonable
    end
    
end

