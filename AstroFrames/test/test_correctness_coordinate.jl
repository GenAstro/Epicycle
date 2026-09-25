# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Converting without a spacecraft.
#
# The astronomy case: a state, a frame, an epoch, and nothing from AstroModels
# anywhere. Everything below is written against AstroFrames and AstroStates
# only, which is the property under test as much as the numbers are.
# =============================================================================

using AstroFrames
using AstroUniverse
using AstroEpochs
using AstroStates
using LinearAlgebra
using Test

const _JD_C  = 2458849.5
const _EP_C  = Time(_JD_C, 0.0, :tdb, :jd)
const _VEC_C = [7000.0, 0.0, 0.0, 0.0, 7.5460491, 0.0]
const _EQ    = CoordinateSystem(earth, MJ2000Eq())

"""A type AstroFrames has never heard of, made convertible by three methods."""
struct _TestSubject end
AstroFrames.state_of(::_TestSubject) = CartesianState(_VEC_C)
AstroFrames.frame_of(::_TestSubject) = _EQ
AstroFrames.epoch_of(::_TestSubject) = _EP_C

@testset "Coordinate construction" begin

    @testset "from a plain vector, which is the astronomy case" begin
        c = Coordinate(_VEC_C, _EQ, _EP_C)
        @test state_of(c) isa CartesianState
        @test frame_of(c) === _EQ
        @test epoch_of(c) === _EP_C
        @test to_vector(state_of(c)) ≈ _VEC_C
    end

    @testset "from a state object" begin
        c = Coordinate(CartesianState(_VEC_C), _EQ, _EP_C)
        @test to_vector(state_of(c)) ≈ _VEC_C
    end

    @testset "it holds a state, not raw numbers" begin
        # The design point: composing an AstroStates state means every other
        # representation comes along, rather than a second state type existing.
        μ = get_gravparam(earth)
        kep = KeplerianState(CartesianState(_VEC_C), μ)
        c = Coordinate(kep, _EQ, _EP_C)
        @test state_of(c) isa KeplerianState

        # And it still converts, because the μ comes from its own origin.
        @test to_vector(CartesianState(c)) ≈ _VEC_C rtol = 1e-9
    end

    @testset "it shows something a user can read" begin
        out = sprint(show, MIME"text/plain"(), Coordinate(_VEC_C, _EQ, _EP_C))
        @test occursin("Coordinate", out)
        @test occursin("Earth", out)
        @test occursin("MJ2000Eq", out)
    end
end

@testset "Coordinate conversion" begin

    @testset "axes only, against the operator tier" begin
        # Must agree exactly with doing it by hand, because it is the same
        # code path with the epoch and source frame filled in.
        c  = Coordinate(_VEC_C, _EQ, _EP_C)
        ec = CoordinateSystem(earth, MJ2000Ec())

        by_hand = axes_rotation(MJ2000Eq(), MJ2000Ec(), _EP_C) * _VEC_C
        @test to_vector(CartesianState(c, ec)) ≈ by_hand

        c2 = Coordinate(c, ec)
        @test frame_of(c2) === ec
        @test epoch_of(c2) === _EP_C          # the epoch is carried, not reset
        @test to_vector(state_of(c2)) ≈ by_hand
    end

    @testset "origin change is applied, not skipped" begin
        # The failure this guards is a silent error the size of the origin
        # separation — 384 000 km, not a rounding difference.
        c  = Coordinate(_VEC_C, _EQ, _EP_C)
        moon_eq = CoordinateSystem(moon, MJ2000Eq())

        r_earth = to_vector(CartesianState(c, _EQ))[1:3]
        r_moon  = to_vector(CartesianState(c, moon_eq))[1:3]
        @test norm(r_moon - r_earth) > 300_000
        @test 350_000 < norm(r_moon) < 420_000
    end

    @testset "round trip" begin
        c   = Coordinate(_VEC_C, _EQ, _EP_C)
        itrf = CoordinateSystem(earth, ITRF())
        back = Coordinate(Coordinate(c, itrf), _EQ)
        @test to_vector(state_of(back)) ≈ _VEC_C rtol = 1e-12
    end

    @testset "a frame defined by another orbit takes params" begin
        c = Coordinate(_VEC_C, _EQ, _EP_C)
        chief = [7100.0, 50.0, 0.0, -0.05, 7.49, 0.0]
        ric = CoordinateSystem(earth, RIC())

        out = to_vector(CartesianState(c, ric, (; reference_state = chief)))
        by_hand = axes_rotation(MJ2000Eq(), RIC(), _EP_C, (; reference_state = chief)) * _VEC_C
        @test out ≈ by_hand
    end

    @testset "converting to the frame it is already in changes nothing" begin
        c = Coordinate(_VEC_C, _EQ, _EP_C)
        @test to_vector(CartesianState(c, _EQ)) ≈ _VEC_C rtol = 1e-14
    end
end

@testset "the subject interface" begin

    @testset "a user's own type works with three methods" begin
        # The claim the design makes: implement three methods and every
        # conversion works. Nothing else about this type is known here.
        # `Coordinate(subject, target)` is the generic path — it takes any
        # subject. A `CartesianState(x, cs)` method belongs to whoever owns `x`,
        # which is why AstroModels defines the Spacecraft one and a user would
        # define theirs.
        c = Coordinate(_TestSubject(), CoordinateSystem(earth, MJ2000Ec()))
        by_hand = axes_rotation(MJ2000Eq(), MJ2000Ec(), _EP_C) * _VEC_C
        @test to_vector(state_of(c)) ≈ by_hand
        @test epoch_of(c) === _EP_C

        # And an origin change works on it too, with nothing else declared.
        c2 = Coordinate(_TestSubject(), CoordinateSystem(moon, MJ2000Eq()))
        @test 350_000 < norm(to_vector(state_of(c2))[1:3]) < 420_000
    end

    @testset "a type that has not implemented it says so by name" begin
        # Not a MethodError several calls down.
        e = try
            Coordinate("not a subject", _EQ)
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("state_of", e.msg)
        @test occursin("Coordinate", e.msg)
    end
end
