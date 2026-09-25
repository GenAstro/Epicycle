# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Origins that carry their own state.
#
# `CoordinateSystem(sc, RIC())` names a spacecraft as the origin, and RIC axes
# are built from an orbit — that orbit, the origin's. Two things used to be
# wrong about it:
#
#   1. the reference orbit was not taken from the origin, so the caller had to
#      pass the origin's own state back in as `reference_state`. Redundant, and
#      worse, nothing checked that the two agreed — you could build the frame
#      from one orbit and centre it on another.
#
#   2. the origin offset went to the ephemeris, which has no entry for a
#      spacecraft. `translate_state(::Spacecraft, ::CelestialBody, ::Float64)`
#      raised a `MethodError` from two calls down.
#
# The tests use a local subject type rather than a `Spacecraft`, because the
# thing being tested is the interface — anything implementing `state_of`,
# `frame_of` and `epoch_of` gets this — and because AstroModels depends on this
# package, not the other way around.
# =============================================================================

using AstroFrames
using AstroUniverse
using AstroStates
using AstroEpochs
using EpicycleBase: AbstractPoint
using LinearAlgebra
using Test

const _SUBJ_EPOCH = Time(2458849.5, 0.0, :tdb, :jd)

"""A minimal origin that carries its own state, as a spacecraft does."""
struct _Chief <: AbstractPoint
    state::CartesianState
    cs::CoordinateSystem
    epoch::Time
end

AstroFrames.state_of(c::_Chief) = c.state
AstroFrames.frame_of(c::_Chief) = c.cs
AstroFrames.epoch_of(c::_Chief) = c.epoch

_subject(v, epoch = _SUBJ_EPOCH) =
    _Chief(CartesianState(v), CoordinateSystem(earth, ICRF()), epoch)

const _CHIEF = _subject([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0])

@testset "carrying a state is what makes an origin one" begin
    # Implementing `state_of` is the whole opt-in — there is no second list to
    # register in. A body does not implement it, so it is still an ephemeris
    # lookup.
    @test AstroFrames.carries_own_state(_CHIEF)
    @test !AstroFrames.carries_own_state(earth)
    @test !AstroFrames.carries_own_state(mars)
end

@testset "the reference orbit comes from the origin" begin
    # No params. The frame is defined by the origin's orbit and the origin is
    # right there.
    for axes in (RIC(), LVLH())
        c = Coordinate(_CHIEF, CoordinateSystem(_CHIEF, axes))
        @test all(iszero, to_vector(state_of(c)))     # the origin is at its own origin
    end
end

@testset "a deputy reads as radial, in-track, cross-track" begin
    # The check that the whole path is right end to end: the reference orbit is
    # the chief's, the offset is applied in the target axes, and the axes are
    # in the order RIC names.
    base = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
    for (offset, expected, what) in (([10.0, 0.0, 0.0], [10.0, 0.0, 0.0], "radial"),
                                     ([0.0, 10.0, 0.0], [0.0, 10.0, 0.0], "in-track"),
                                     ([0.0, 0.0, 10.0], [0.0, 0.0, 10.0], "cross-track"))
        deputy = _subject(vcat(base[1:3] .+ offset, base[4:6]))
        r = to_vector(state_of(Coordinate(deputy, CoordinateSystem(_CHIEF, RIC()))))[1:3]
        @testset "$(what)" begin
            @test r ≈ expected atol = 1e-9
        end
    end
end

@testset "an explicit reference orbit is not overridden" begin
    # Passing one yourself still works and still wins. Filling in from the
    # origin is a default, not a policy.
    #
    # It has to be read off a deputy: the chief about its own origin is at the
    # origin, and rotating zero by either set of axes is still zero.
    deputy = _subject([7010.0, 0.0, 0.0, 0.0, 7.546, 0.0])

    # The same orbit a quarter turn on, so its radial axis is the chief's
    # in-track and the 10 km radial offset lands in a different component.
    quarter = [0.0, 7000.0, 0.0, -7.546, 0.0, 0.0]

    auto = to_vector(state_of(Coordinate(deputy, CoordinateSystem(_CHIEF, RIC()))))[1:3]
    mine = to_vector(state_of(Coordinate(deputy, CoordinateSystem(_CHIEF, RIC()),
                                         (; reference_state = quarter))))[1:3]

    @test auto ≈ [10.0, 0.0, 0.0] atol = 1e-9      # from the origin
    @test !isapprox(mine, auto; atol = 1e-6)       # the supplied one was used
    @test norm(mine) ≈ 10.0 atol = 1e-9            # same offset, different axes
end

@testset "a body origin cannot supply one, and says so" begin
    # Nothing is filled in, and the frame reports what it needs — the behaviour
    # before any of this, unchanged.
    message = try
        Coordinate(_CHIEF, CoordinateSystem(earth, RIC()))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("reference_state", message)
    @test occursin("defined by a reference orbit", message)
end

@testset "VNB works from the origin too, unaccelerated" begin
    # The state comes from the origin; the acceleration does not, because no
    # origin carries one. Omitting it is not a gap — it says the reference is
    # not being accelerated, and VNB about such a reference does not turn.
    c = Coordinate(_CHIEF, CoordinateSystem(_CHIEF, VNB()))
    @test all(iszero, to_vector(state_of(c)))

    # Supplying it changes the rate and nothing else.
    μ = get_gravparam(earth)
    r = to_vector(state_of(_CHIEF))[1:3]
    accel = (-μ / norm(r)^3) .* r

    deputy = _subject([7010.0, 0.0, 0.0, 0.0, 7.546, 0.0])
    without = to_vector(state_of(Coordinate(deputy, CoordinateSystem(_CHIEF, VNB()))))
    with    = to_vector(state_of(Coordinate(deputy, CoordinateSystem(_CHIEF, VNB()),
                                            (; reference_accel = accel))))

    @test without[1:3] ≈ with[1:3] atol = 1e-12      # same axes, same position
    @test !isapprox(without[4:6], with[4:6]; atol = 1e-9)   # different rate
end

@testset "an origin's state belongs to one epoch" begin
    # The origin carries its state at its own epoch. Using it at another would
    # build the frame from an orbit the origin has since left, and be wrong by
    # however far it moved — silently. Made loud.
    later = _subject([7010.0, 0.0, 0.0, 0.0, 7.546, 0.0],
                     Time(2458850.5, 0.0, :tdb, :jd))     # one day on

    message = try
        Coordinate(later, CoordinateSystem(_CHIEF, RIC()))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("epoch", message)
    @test occursin("86400", replace(message, "." => ""))   # names the gap, in seconds
    @test occursin("reference_state", message)             # and the way round it
end

@testset "origin_translation reads a carried state, not the ephemeris" begin
    # The offset is the source origin seen from the target. With the chief as
    # target and Earth as source, that is minus the chief's own state.
    Δ = origin_translation(earth, _CHIEF, ICRF(), _SUBJ_EPOCH)
    @test Δ ≈ -to_vector(state_of(_CHIEF)) atol = 1e-9

    # And the reverse direction negates, as for any pair of origins.
    Δ_rev = origin_translation(_CHIEF, earth, ICRF(), _SUBJ_EPOCH)
    @test Δ_rev ≈ -Δ atol = 1e-9
end
