
using Test
using LinearAlgebra
using AstroFrames
using AstroProp
using AstroEpochs
using AstroStates
using AstroManeuvers
using AstroUniverse
using AstroCallbacks
using AstroSolve

@testset "Event Propagation Test" begin
    # Reset spacecraft to initial state for independent computation
    sat1_reset = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time = Time("2020-09-21T12:23:12", TAI(), ISOT())
    )
    sat2_reset = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time = Time("2020-09-21T12:23:12", TAI(), ISOT())
    )
    dynsys1 = DynSys(forces = forces, spacecraft = [sat1_reset])
    dynsys2 = DynSys(forces = forces, spacecraft = [sat2_reset])

    # --- Event infrastructure ---
    toi_event = ImpulsiveManeuver(
        axes = VNB(),
        element1 = 0.4,
        element2 = 0.0,
        element3 = 0.3
    )
    var_toi = SolverVariable(
        calc = ManeuverCalc(toi_event, sat1_reset, DeltaVVector()),
        name = "toi", 
        lower_bound = [-1.0, 0.0, 0.0], 
        upper_bound = [1.0, 0.0, 0.0])

    apply_toi = Event(event = () -> maneuver(sat1_reset, toi_event), vars = [var_toi])
    prop_to_moi() = propagate(dynsys1, integ, StopAtApoapsis(sat1_reset))
    prop_moi = Event(event = prop_to_moi)

    # Apply maneuver and propagate using event infrastructure
    apply_event(apply_toi)
    apply_event(prop_moi)
    #event_result = copy(CartesianState(sat1_reset.state).posvel)
    event_result = copy(to_posvel(sat1_reset))

    # --- Independent computation ---
    maneuver(sat2_reset, toi_event)
    propagate(dynsys2, integ, StopAtApoapsis(sat2_reset))
    truth_result = copy(to_posvel(sat2_reset))

    # --- Compare results ---
    @test isapprox(event_result, truth_result; rtol=1e-12, atol=1e-12)
end