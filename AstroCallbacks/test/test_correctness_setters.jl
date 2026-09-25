# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Writing quantities.
#
# A settable quantity has a setter named for it, `state!(sat; to = x)`, and a `set_quantity!`
# method a solver calls because it holds the quantity as a value. Both reach the same code. An
# element of a representation is written by `set_element!`, which converts the state into that
# representation in the coordinate system given, replaces the one element and converts back, and
# stores the result in the subject's own frame through `set_state!`.
#
# Truth: **round trip**. What is written reads back, in the coordinate system it was written in,
# and the rest of its representation reads back unchanged. The state stays in the subject's own
# frame throughout, so a write in another frame is checked by reading in both.

using Test
using LinearAlgebra
using AstroCallbacks
using AstroCallbacks: set_element!, set_state!
using AstroStates
using AstroEpochs
using AstroFrames
using AstroModels
using AstroUniverse
using AstroManeuvers
using EpicycleBase

const _W_EPOCH = Time(2458849.5, 0.0, TDB(), JD())

_w_sat() = Spacecraft(CartesianState([7000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]), _W_EPOCH;
                      coord_sys = EarthMJ2000Eq)

_kep(sat, cs) = [element(KeplerianState, sat, cs, f) for f in (:sma, :ecc, :inc, :raan, :aop, :ta)]

@testset "state! replaces the state, and the bridge reaches it" begin
    sat = _w_sat()
    x   = [7100.0, 0.0, 0.0, 0.0, 7.4, 0.5]
    @test state!(sat; to = x) === sat
    @test state(sat) ≈ x

    y = [7200.0, 100.0, 0.0, 0.0, 7.3, 0.4]
    set_quantity!(sat, state; to = y)
    @test state(sat) ≈ y
    @test is_settable(sat, state)
end

@testset "set_calc! writes through the bridge, with the reader's own arguments" begin
    sat = _w_sat()
    c   = Calc(state, sat)
    x   = [6900.0, 0.0, 0.0, 0.0, 7.6, 0.0]
    set_calc!(c, x)
    @test c() ≈ x
end

@testset "set_element! writes one element and holds the rest" begin
    sat    = _w_sat()
    before = _kep(sat, EarthMJ2000Eq)

    @test set_element!(KeplerianState, sat, EarthMJ2000Eq, :sma; to = 8000.0) === sat
    after = _kep(sat, EarthMJ2000Eq)
    @test after[1] ≈ 8000.0 rtol = 1e-12
    @test after[2:6] ≈ before[2:6] rtol = 1e-10
    @test frame_of(sat) === EarthMJ2000Eq          # stored where it was
end

@testset "set_element! in another frame holds the rest there, and stores in the own frame" begin
    # The spacecraft keeps an equatorial state. Its ecliptic inclination is set to 30°, so the
    # other ecliptic elements must hold, and the equatorial ones must move.
    sat      = _w_sat()
    ec0      = _kep(sat, EarthMJ2000Ec)
    eq0      = _kep(sat, EarthMJ2000Eq)

    set_element!(KeplerianState, sat, EarthMJ2000Ec, :inc; to = deg2rad(30.0))
    ec1 = _kep(sat, EarthMJ2000Ec)
    @test ec1[3] ≈ deg2rad(30.0) rtol = 1e-12
    @test ec1[[1, 2, 4, 5, 6]] ≈ ec0[[1, 2, 4, 5, 6]] rtol = 1e-10
    @test frame_of(sat) === EarthMJ2000Eq
    @test !(_kep(sat, EarthMJ2000Eq)[3] ≈ eq0[3])  # the equatorial inclination changed
end

@testset "set_state! stores a state given in another frame in the subject's own" begin
    sat  = _w_sat()
    cart = CartesianState(Coordinate(sat, EarthMJ2000Ec, NamedTuple()))
    x0   = state(sat)
    set_state!(sat, cart, EarthMJ2000Ec)          # the same state, expressed elsewhere
    @test state(sat) ≈ x0 rtol = 1e-12
    @test frame_of(sat) === EarthMJ2000Eq
end

# --- the Keplerian setters ------------------------------------------------------------------
#
# The invariant is the one test_correctness_element_setters.jl holds the legacy setters to:
# setting one element changes that element and nothing else. Written over all six, so a
# permutation error in any of them fails here.

const _KEP_SETTERS = ((semi_major_axis!,       1, 8000.0),
                      (eccentricity!,          2, 0.05),
                      (inclination!,           3, 0.9),
                      (raan!,                  4, 1.2),
                      (argument_of_periapsis!, 5, 2.0),
                      (true_anomaly!,          6, 0.4))

@testset "each Keplerian setter writes its element and holds the other five" begin
    for (set!, k, value) in _KEP_SETTERS
        @testset "$(nameof(set!))" begin
            sat    = _w_sat()
            before = _kep(sat, EarthMJ2000Eq)
            @test set!(sat; to = value) === sat
            after = _kep(sat, EarthMJ2000Eq)
            @test after[k] ≈ value rtol = 1e-10
            others = setdiff(1:6, k)
            @test after[others] ≈ before[others] rtol = 1e-9
        end
    end
end

@testset "a Keplerian setter in another frame holds the rest there" begin
    sat = _w_sat()
    ec0 = _kep(sat, EarthMJ2000Ec)
    raan!(sat, EarthMJ2000Ec; to = 1.2)
    ec1 = _kep(sat, EarthMJ2000Ec)
    @test ec1[4] ≈ 1.2 rtol = 1e-10
    @test ec1[[1, 2, 3, 5, 6]] ≈ ec0[[1, 2, 3, 5, 6]] rtol = 1e-9
    @test frame_of(sat) === EarthMJ2000Eq
end

@testset "the bridge is the setter, and settability is the spacecraft's" begin
    a = _w_sat(); b = _w_sat()
    semi_major_axis!(a, EarthMJ2000Ec; to = 7500.0)
    set_quantity!(b, semi_major_axis, EarthMJ2000Ec; to = 7500.0)
    @test state(a) == state(b)

    @test is_settable(a, semi_major_axis)
    @test is_settable(a, semi_major_axis, EarthMJ2000Ec)
    @test !is_settable(a, semi_major_axis, 1.0)                 # not a coordinate system
    coord = Coordinate(a, EarthMJ2000Eq, NamedTuple())
    @test !is_settable(coord, semi_major_axis)                   # a Coordinate is a value

    # And a solver writes it through a Calc with the reader's own arguments.
    c = Calc(semi_major_axis, a, EarthMJ2000Ec)
    set_calc!(c, 7600.0)
    @test c() ≈ 7600.0 rtol = 1e-10
end

@testset "elements that describe no orbit are refused, and the state is not touched" begin
    hyperbolic() = Spacecraft(CartesianState([7000.0, 0.0, 0.0, 0.0, 12.0, 0.0]), _W_EPOCH;
                              coord_sys = EarthMJ2000Eq)
    cases = (("negative eccentricity",     _w_sat,     eccentricity!,    -0.1),
             ("eccentricity of one",       _w_sat,     eccentricity!,     1.0),
             ("ellipse, negative sma",     _w_sat,     semi_major_axis!, -7000.0),
             ("hyperbola, positive sma",   hyperbolic, semi_major_axis!,  7000.0),
             ("inclination beyond π",      _w_sat,     inclination!,      4.0),
             ("negative inclination",      _w_sat,     inclination!,     -0.5),
             ("beyond the asymptote",      hyperbolic, true_anomaly!,     3.0))
    for (name, make, set!, value) in cases
        @testset "$name" begin
            sat = make()
            x0  = copy(state(sat))
            @test_throws ArgumentError set!(sat; to = value)
            @test state(sat) == x0
        end
    end
end

@testset "epoch! moves the epoch and keeps the state's numbers" begin
    sat = _w_sat()
    x0  = copy(state(sat))
    t1  = Time(2458850.5, 0.0, TDB(), JD())
    @test epoch!(sat; to = t1) === sat
    @test epoch(sat) == t1
    @test state(sat) == x0

    t2 = Time(2458851.5, 0.0, TDB(), JD())
    set_quantity!(sat, epoch; to = t2)
    @test epoch(sat) == t2
    @test is_settable(sat, epoch)
end

# --- the Cartesian setters ------------------------------------------------------------------

const _CART_SETTERS = ((position_x!,      1:1, 7100.0),
                       (position_y!,      2:2, -500.0),
                       (position_z!,      3:3, 1500.0),
                       (position_vector!, 1:3, [7100.0, -500.0, 1500.0]),
                       (velocity_vector!, 4:6, [1.0, 7.0, 0.5]))

_cart(sat, cs) = collect(to_vector(CartesianState(Coordinate(sat, cs, NamedTuple()))))

@testset "each Cartesian setter writes its components and holds the rest, in any frame" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec), (set!, idx, value) in _CART_SETTERS
        @testset "$(nameof(set!)) in $(cs === EarthMJ2000Eq ? "Eq" : "Ec")" begin
            sat    = _w_sat()
            before = _cart(sat, cs)
            @test set!(sat, cs; to = value) === sat
            after = _cart(sat, cs)
            @test after[idx] ≈ vcat(value) rtol = 1e-12
            others = setdiff(1:6, idx)
            @test after[others] ≈ before[others] rtol = 1e-12
            @test frame_of(sat) === EarthMJ2000Eq
        end
    end
end

@testset "the Cartesian bridges, and what they refuse" begin
    a = _w_sat(); b = _w_sat()
    position_z!(a, EarthMJ2000Ec; to = 900.0)
    set_quantity!(b, position_z, EarthMJ2000Ec; to = 900.0)
    @test state(a) == state(b)
    @test is_settable(a, velocity_vector, EarthMJ2000Ec)

    sat = _w_sat(); x0 = copy(state(sat))
    @test_throws ArgumentError position_vector!(sat; to = [1.0, 2.0])
    @test_throws ArgumentError velocity_vector!(sat; to = [NaN, 0.0, 0.0])
    @test state(sat) == x0
end

# --- the mean semi-major axis, the magnitudes, and the delta-V magnitude --------------------

_leo() = Spacecraft(CartesianState(KeplerianState(7000.0, 0.01, 0.9, 0.3, 0.2, 0.1), 398600.4418),
                    _W_EPOCH; coord_sys = EarthMJ2000Eq)
_mean(sat) = [element(BrouwerMeanLongState, sat, EarthMJ2000Eq, f)
              for f in (:sma, :ecc, :inc, :raan, :aop, :ma)]

@testset "mean_long_sma! writes the mean semi-major axis and holds the other mean elements" begin
    # The mean-to-osculating conversion iterates to a tolerance, so this is not machine precision.
    # Measured: sma within 2.4e-6 km; ecc, inc and raan within 2e-9; aop and ma within 6e-7, and
    # opposite, because periapsis is poorly defined on a near-circular orbit.
    sat = _leo()
    m0  = _mean(sat)
    @test mean_long_sma!(sat; to = 7100.0) === sat
    m1 = _mean(sat)
    @test m1[1] ≈ 7100.0 atol = 1e-4
    @test m1[2:4] ≈ m0[2:4] atol = 1e-7
    @test m1[5:6] ≈ m0[5:6] atol = 1e-5
    @test m1[5] + m1[6] ≈ m0[5] + m0[6] atol = 1e-8        # the argument of latitude holds tightly

    b = _leo()
    set_quantity!(b, mean_long_sma; to = 7100.0)
    @test state(b) == state(sat)

    # Outside the theory's domain the write is refused.
    x0 = copy(state(sat))
    @test_throws ArgumentError mean_long_sma!(sat; to = 2000.0)     # periapsis below 3000 km
    @test state(sat) == x0
end

@testset "a magnitude setter scales its vector and holds the direction and the other vector" begin
    for cs in (EarthMJ2000Eq, EarthMJ2000Ec)
        sat = _w_sat()
        r0, v0 = position_vector(sat, cs), velocity_vector(sat, cs)
        position_magnitude!(sat, cs; to = 8000.0)
        r1, v1 = position_vector(sat, cs), velocity_vector(sat, cs)
        @test norm(r1) ≈ 8000.0 rtol = 1e-12
        @test r1 / norm(r1) ≈ r0 / norm(r0) rtol = 1e-12
        @test v1 ≈ v0 rtol = 1e-12

        sat = _w_sat()
        r0, v0 = position_vector(sat, cs), velocity_vector(sat, cs)
        velocity_magnitude!(sat, cs; to = 7.0)
        r1, v1 = position_vector(sat, cs), velocity_vector(sat, cs)
        @test norm(v1) ≈ 7.0 rtol = 1e-12
        @test v1 / norm(v1) ≈ v0 / norm(v0) rtol = 1e-12
        @test r1 ≈ r0 rtol = 1e-12
        @test frame_of(sat) === EarthMJ2000Eq
    end

    a = _w_sat(); b = _w_sat()
    velocity_magnitude!(a, EarthMJ2000Ec; to = 7.2)
    set_quantity!(b, velocity_magnitude, EarthMJ2000Ec; to = 7.2)
    @test state(a) == state(b)
end

@testset "a magnitude with no direction, or a negative one, is refused" begin
    still = Spacecraft(CartesianState([7000.0, 0.0, 0.0, 0.0, 0.0, 0.0]), _W_EPOCH;
                       coord_sys = EarthMJ2000Eq)
    @test_throws ArgumentError velocity_magnitude!(still; to = 7.5)
    @test state(still) == [7000.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    sat = _w_sat(); x0 = copy(state(sat))
    @test_throws ArgumentError position_magnitude!(sat; to = -1.0)
    @test_throws ArgumentError velocity_magnitude!(sat; to = -1.0)
    @test state(sat) == x0
end

@testset "delta_v_magnitude! scales the delta-V and holds its direction" begin
    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = -0.3)
    d0  = delta_v(man)
    @test delta_v_magnitude!(man; to = 0.5) === man
    @test delta_v_magnitude(man) ≈ 0.5 rtol = 1e-14
    @test delta_v(man) / 0.5 ≈ d0 / norm(d0) rtol = 1e-14

    set_quantity!(man, delta_v_magnitude; to = 0.25)
    @test delta_v_magnitude(man) ≈ 0.25 rtol = 1e-14

    zero_burn = ImpulsiveManeuver(axes = VNB())
    @test_throws ArgumentError delta_v_magnitude!(zero_burn; to = 0.1)
    @test_throws ArgumentError delta_v_magnitude!(man; to = -0.1)
    @test delta_v_magnitude(man) ≈ 0.25 rtol = 1e-14
end
