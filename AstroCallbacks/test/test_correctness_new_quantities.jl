# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# The four quantities that had a calc variable and no function.
#
# The property under test is that each function returns what its calc returns.
# A constraint written either way has to mean the same thing, because these
# exist so the calc form can be retired.
# =============================================================================

using Test
using AstroCallbacks
using AstroManeuvers
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: label, is_settable, is_cyclic, set_quantity!
using LinearAlgebra: norm

# `get_calc` hands back a length-1 vector for a scalar quantity; the function
# form hands back the number. Compare like with like.
_q_scalar(x) = x isa AbstractVector ? only(x) : x

_q_epoch() = Time("2020-09-21T12:23:12", TAI(), ISOT())
_q_sat()   = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                        time  = _q_epoch())
# A hyperbolic state, so the outgoing asymptote is defined.
_q_hyp()   = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]),
                        time  = _q_epoch())

@testset "mean_long_sma matches the MeanSMA calc" begin
    sat = _q_sat()
    @test mean_long_sma(sat) ≈ _q_scalar(get_calc(OrbitCalc(sat, MeanSMA())))
    @test label(mean_long_sma) == "Long-period mean semi-major axis"
end

@testset "outgoing_rla matches the OutGoingRLA calc" begin
    sat = _q_hyp()
    @test outgoing_rla(sat) ≈ _q_scalar(get_calc(OrbitCalc(sat, OutGoingRLA())))
    # A right ascension wraps, and a stopping condition has to know that.
    @test is_cyclic(outgoing_rla)
end

@testset "delta_v_magnitude matches the DeltaVMag calc" begin
    sat = _q_sat()
    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = -0.3)
    @test delta_v_magnitude(man) ≈ _q_scalar(get_calc(ManeuverCalc(man, sat, DeltaVMag())))
    @test delta_v_magnitude(man) ≈ norm(delta_v(man))
    # Settable, by scaling the delta-V along its direction.
    @test is_settable(man, delta_v_magnitude)
end

@testset "gravitational_parameter matches the GravParam calc, and writes" begin
    @test gravitational_parameter(earth) ≈ get_calc(BodyCalc(earth, GravParam()))
    @test is_settable(earth, gravitational_parameter)

    mu0 = gravitational_parameter(earth)
    try
        set_quantity!(earth, gravitational_parameter; to = 3.9e5)
        @test gravitational_parameter(earth) ≈ 3.9e5
        gravitational_parameter!(earth; to = 3.95e5)
        @test gravitational_parameter(earth) ≈ 3.95e5
    finally
        gravitational_parameter!(earth; to = mu0)            # shared global body
    end
    @test gravitational_parameter(earth) ≈ mu0
end

# A module standing in for a script: a quantity defined there, with a subject type, is found there.
module _ScriptQuantities
    using AstroCallbacks, AstroModels
    import EpicycleBase
    drag_area(sc) = sc.drag.drag_area
    EpicycleBase.label(::typeof(drag_area)) = "Drag area"
    AstroCallbacks.subject_type(::typeof(drag_area)) = Spacecraft
end

@testset "quantities lists the quantities that read a subject type" begin
    # Truth: the quantities table in docs/src/index.md, by subject.
    @test quantities(ImpulsiveManeuver) == [delta_v, delta_v_magnitude]
    @test quantities(CelestialBody)     == [gravitational_parameter]

    on_sc = quantities(Spacecraft)
    @test state in on_sc && semi_major_axis in on_sc && position_dot_velocity in on_sc
    @test !(delta_v in on_sc) && !(gravitational_parameter in on_sc)
    @test length(on_sc) == 18

    # A bare state reads the orbital quantities, but it is not a spacecraft, so not `state`.
    on_coord = quantities(Coordinate)
    @test inclination in on_coord && !(state in on_coord)

    @test length(quantities()) == 21
    @test issorted(String.(nameof.(quantities())))

    # A script's own quantity, found in the script's module.
    @test quantities(Spacecraft, _ScriptQuantities) == [_ScriptQuantities.drag_area]
end
