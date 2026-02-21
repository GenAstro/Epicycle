using LinearAlgebra
using AstroFrames
using AstroProp
using AstroEpochs
using AstroStates
using AstroManeuvers
using AstroSolve
using AstroUniverse
using Test

# Create spacecraft
sat1 = Spacecraft(
            state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]), 
            time = Time("2020-09-21T12:23:12", TAI(), ISOT())
            )

# Create force models, integrator, and dynamics system
pm_grav = PointMassGravity(earth,(moon,sun))
forces = ForceModel(pm_grav)
integ = IntegratorConfig(DP8(); abstol = 1e-11, reltol = 1e-11, dt = 4000)

# Define which spacecraft to propagate and which force model to use
dynsys = DynSys(
          forces = forces, 
          spacecraft = [sat1]
          )

# Create maneuver models for the hohmann transfer
toi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.1,
    element2 = 0.2,
    element3 = 0.3
)

moi = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.4,
    element2 = 0.5,
    element3 = 0.6
)

# Define toi as a solver variable
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sat1, DeltaVVector()),
    name = "toi",
    lower_bound = [-1.0, -2.0, -3.0],
    upper_bound = [1.0, 2.0, 3.0],
    shift = [0.01, 0.02, 0.03],
    scale = [1.1, 1.2, 1.3]
)

# Define moi as a solver variable
var_moi = SolverVariable(
    calc = ManeuverCalc(moi, sat1, DeltaVVector()),
    name = "moi",
    lower_bound = [-4.0, -5.0, -6.0],
    upper_bound = [4.0, 5.0, 6.0],
    shift = [0.04, 0.05, 0.06],
    scale = [1.4, 1.5, 1.6]
)

# Create the MOI Event
toi_fun() = maneuver!(sat1, toi) 
toi_event = Event(name = "toi", event = toi_fun, vars = [var_toi])

# Create the prop to apopasis event
prop_apo_fun() = propagate!(dynsys, integ, StopAtApoapsis(sat1))
prop_event = Event(name = "prop_apo", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver!(sat1, moi)
moi_event = Event(event = moi_fun, vars = [var_moi])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

sm = SequenceManager(seq)
get_var_values(sm.ordered_vars)
get_var_shifts(sm.ordered_vars)
get_var_scales(sm.ordered_vars)
get_var_lower_bounds(sm.ordered_vars)
get_var_upper_bounds(sm.ordered_vars)

@testset "SequenceManager variable vector assembly" begin
    # Should be in the order: var_toi, var_moi
    @test sm.ordered_vars == [var_toi, var_moi]

    # Test shifts
    @test get_var_shifts(sm.ordered_vars) == vcat(var_toi.shift, var_moi.shift)

    # Test scales
    @test get_var_scales(sm.ordered_vars) == vcat(var_toi.scale, var_moi.scale)

    # Test lower bounds
    @test get_var_lower_bounds(sm.ordered_vars) == vcat(var_toi.lower_bound, var_moi.lower_bound)

    # Test upper bounds
    @test get_var_upper_bounds(sm.ordered_vars) == vcat(var_toi.upper_bound, var_moi.upper_bound)
end

@testset "set_var_values and get_var_values roundtrip" begin
    # Set the variable vector
    test_vec = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
    set_var_values(test_vec, sm.ordered_vars)
    # Now get the values back out and check they match
    @test get_var_values(sm.ordered_vars) == test_vec
end

pos_target = 45000.0
pos_con = Constraint(
    calc = OrbitCalc(sat1, PosMag()),
    lower_bounds= [pos_target],
    upper_bounds=[pos_target],
    scale = [1.0],
)

vel_target = sqrt(earth.mu / pos_target)
vel_con = Constraint(
    calc = OrbitCalc(sat1, VelMag()),
    lower_bounds = [vel_target],
    upper_bounds = [vel_target],
    scale = [1.0],
)

# Create the MOI Event
toi_fun() = maneuver!(sat1, toi) 
toi_event = Event(name = "toi", 
                  event = toi_fun, 
                  vars = [var_toi],
                  funcs = [])

# Create the prop to apopasis event
prop_apo_fun() = propagate!(dynsys, integ, StopAtApoapsis(sat1))
prop_event = Event(name = "prop_apo", event = prop_apo_fun)

# Create the TOI event. 
moi_fun() = maneuver!(sat1, moi)

moi_event = Event(event = moi_fun, vars = [var_moi],
                  funcs = [pos_con, vel_con])

# Build sequence and solve
seq = Sequence()
add_events!(seq, prop_event, [toi_event]) 
add_events!(seq, moi_event, [prop_event])

sm = SequenceManager(seq)

@testset "SequenceManager ordered_funcs and fun_sizes" begin
    # The only event with constraints is moi_event, which has [pos_con, vel_con]
    expected_funcs = [pos_con, vel_con]
    expected_sizes = [pos_con.numvars, vel_con.numvars]

    @test sm.ordered_funcs == expected_funcs
    @test sm.fun_sizes == expected_sizes
end

@testset "SequenceManager get_fun_values" begin
    # Compute the expected function values using the constraints directly
    expected_vals = vcat(func_eval(pos_con), func_eval(vel_con))
    # Use the SequenceManager method
    actual_vals = get_fun_values(sm)
    @test actual_vals == expected_vals
end

# Test that setting decision variables takes and events are applied with new values

include("config_hohmann_transfer.jl")

# Create maneuver models for the hohmann transfer
toi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.15,
    element2 = 0.25,
    element3 = 0.35
)

moi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.45,
    element2 = 0.55,
    element3 = 0.65
)


@testset "solver_fun! sets variables and applies events as expected" begin
    # Set up x to match toi2 and moi2
    x = [toi2.element1, toi2.element2, toi2.element3,
         moi2.element1, moi2.element2, moi2.element3]

    F = similar(get_fun_values(sm))

    # Save sat1's state before
    state_before = to_posvel(sat1)

    # Call solver_fun! (this sets vars and applies events to sat1)
    solver_fun!(F, x, sm)
    state_after = to_posvel(sat1)

    # Now apply the same sequence to sat2 directly
    maneuver!(sat2, toi2)
    propagate!(dynsys2, integ, StopAtApoapsis(sat2))
    maneuver!(sat2, moi2)
    state_sat2 = to_posvel(sat2)

    # The states should match after both sequences
    @test isapprox(state_after, state_sat2; atol=1e-14)

end


include("config_hohmann_transfer.jl")

# Create maneuver models for the hohmann transfer
toi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.15,
    element2 = 0.25,
    element3 = 0.35
)

moi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.45,
    element2 = 0.55,
    element3 = 0.65
)

@testset "solver_fun! resets stateful structs before event application" begin
    # Set up x to match toi2 and moi2
    x = [toi2.element1, toi2.element2, toi2.element3,
         moi2.element1, moi2.element2, moi2.element3]

    F = similar(get_fun_values(sm))

    # Save sat1's state before
    state_before = copy(to_posvel(sat1))

    # Call solver_fun! (this sets vars and applies events to sat1)
    solver_fun!(F, x, sm)
    state_after = copy(to_posvel(sat1))

    # Now reset stateful structs and check state matches original
    reset_stateful_structs!(sm)
    state_reset = copy(to_posvel(sat1))

    @test isapprox(state_before, state_reset; atol=1e-14)

end

include("config_hohmann_transfer.jl")

# Create maneuver models for the hohmann transfer
toi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.15,
    element2 = 0.25,
    element3 = 0.35
)

moi2 = ImpulsiveManeuver(
    axes = VNB(),
    element1 = 0.45,
    element2 = 0.55,
    element3 = 0.65
)

@testset "solver_fun! resets stateful structs between calls" begin
    # Initial input vector
    x1 = [toi2.element1, toi2.element2, toi2.element3,
          moi2.element1, moi2.element2, moi2.element3]
    F1 = similar(get_fun_values(sm))*0.0
    solver_fun!(F1, x1, sm)
    state1 = copy(to_posvel(sat1))
    F1_copy = copy(F1)

    # Different input vector
    x2 = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    F2 = similar(get_fun_values(sm))*0.0
    solver_fun!(F2, x2, sm)
    state2 = copy(to_posvel(sat1))

    # Call again with the original input vector
    F3 = similar(get_fun_values(sm))*0.0
    solver_fun!(F3, x1, sm)
    state3 = copy(to_posvel(sat1))

    # The outputs and states for the same input should match
    @test F1_copy == F3
    @test isapprox(state1, state3; atol=1e-14)
    # And the state for the different input should differ
    @test !isapprox(state1, state2; atol=1e-14)
end

# Example unit test for solver_callback
include("config_hohmann_transfer.jl")

@testset "solver_callback output type and length" begin
    sm = SequenceManager(seq)
    x0 = get_var_values(sm.ordered_vars)
    # Apply all events once to set up constraints/objectives
    for event in sm.sorted_events
        apply_event(event)
    end
    F_template = similar(get_fun_values(sm))

    solver_callback = (F, x) -> solver_fun!(F, x, sm)
    solver_callback(F_template, x0)

    # Type and length checks
    @test typeof(F_template) == Vector{Float64}
    @test length(F_template) == length(get_fun_values(sm))
    @test all(isfinite, F_template)
end