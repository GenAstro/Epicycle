# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0


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
using AstroSolve: apply_event

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
    prop1 = OrbitPropagator(forces, integ)
    prop2 = OrbitPropagator(forces, integ)

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

    apply_toi = Event(event = () -> maneuver!(sat1_reset, toi_event), vars = [var_toi])
    prop_to_moi() = propagate!(prop1, sat1_reset, StopAt(sat1_reset, PosDotVel(), 0.0; direction = -1))
    prop_moi = Event(event = prop_to_moi)

    # Apply maneuver and propagate using event infrastructure
    apply_event(apply_toi)
    apply_event(prop_moi)
    #event_result = copy(CartesianState(sat1_reset.state).posvel)
    event_result = copy(to_posvel(sat1_reset))

    # --- Independent computation ---
    maneuver!(sat2_reset, toi_event)
    propagate!(prop2, sat2_reset, StopAt(sat2_reset, PosDotVel(), 0.0; direction = -1))
    truth_result = copy(to_posvel(sat2_reset))

    # --- Compare results ---
    @test isapprox(event_result, truth_result; rtol=1e-12, atol=1e-12)
end