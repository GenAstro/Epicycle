# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# `element` is a quantity like any other (#11).
#
# Every consumer of a quantity, history, Calc, StopAt, Vary and Constraint, takes the subject from
# the first argument. `element` takes the spacecraft first so it works in all of them, and a Calc
# of it reports by representation and field, with that field's unit and wrap. These tests hold
# that contract; AstroProp and AstroSolve test the propagation and solver ends.
# =============================================================================

using Test
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: label, unit, is_cyclic, cycle, is_settable, set_quantity!

const _EQ_EPOCH = Time("2024-01-01T12:00:00", UTC(), ISOT())
_eq_sat() = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]),
                       time = _EQ_EPOCH)

@testset "element reads with the subject first, in any frame" begin
    sat = _eq_sat()
    @test element(sat, KeplerianState, :sma) == semi_major_axis(sat)
    @test element(sat, KeplerianState, EarthMJ2000Ec, :inc) ≈ inclination(sat, EarthMJ2000Ec)
    @test element(sat, ModifiedKeplerianState, :rp) ≈ periapsis_radius(sat)
    @test element(sat, DelaunayState, :G) ≈ angular_momentum_magnitude(sat) rtol = 1e-13
    @test element(sat, CartesianState, :position) ≈ position_vector(sat)
    # The same read through a deferred Calc
    @test Calc(element, sat, KeplerianState, :sma)() == semi_major_axis(sat)
    @test Calc(element, sat, KeplerianState, EarthMJ2000Ec, :inc)() ≈ inclination(sat, EarthMJ2000Ec)
end

@testset "a Calc of element reports its representation, field, unit and wrap" begin
    sat = _eq_sat()
    L  = Calc(element, sat, DelaunayState, :L)
    l  = Calc(element, sat, DelaunayState, EarthMJ2000Eq, :l)
    h  = Calc(element, sat, EquinoctialState, :h)
    ML = Calc(element, sat, ModifiedEquinoctialState, :L)
    @test label(L) == "Delaunay L"   && unit(L) == "km²/s" && !is_cyclic(L)
    @test label(l) == "Delaunay l"   && unit(l) == "rad"   && is_cyclic(l) && cycle(l) == 2π
    @test label(h) == "Equinoctial h" && unit(h) == ""     && !is_cyclic(h)
    # The same field name means something else in another representation
    @test label(ML) == "ModifiedEquinoctial L" && unit(ML) == "rad" && is_cyclic(ML)
    @test unit(Calc(element, sat, OutGoingAsymptoteState, :c3)) == "km²/s²"
    @test unit(Calc(element, sat, CartesianState, :velocity)) == "km/s"
    @test unit(Calc(element, sat, BrouwerMeanShortState, :ma)) == "rad"
    # Every other quantity's Calc passes its traits straight through
    @test label(Calc(semi_major_axis, sat)) == "Semi-major axis"
    @test unit(Calc(semi_major_axis, sat)) == "km"
    @test label(element) == "Orbit element"
    @test element in quantities(Spacecraft)
end

@testset "every representation's fields have a unit entry" begin
    reps = (CartesianState, KeplerianState, ModifiedKeplerianState, SphericalRADECState,
            SphericalAZIFPAState, ModifiedEquinoctialState, EquinoctialState,
            AlternateEquinoctialState, OutGoingAsymptoteState, IncomingAsymptoteState,
            BrouwerMeanLongState, BrouwerMeanShortState, DelaunayState)
    for S in reps, f in fieldnames(S)
        @test haskey(AstroCallbacks._ELEMENT_TRAITS, (nameof(S), f))
    end
end

@testset "element! writes, holding the rest of the representation" begin
    sat = _eq_sat()
    @test is_settable(sat, element, KeplerianState, :sma)
    @test is_settable(sat, element, KeplerianState, EarthMJ2000Ec, :inc)
    @test !is_settable(Coordinate(sat, EarthMJ2000Eq), element, KeplerianState, :sma)

    ecc0 = eccentricity(sat)
    element!(sat, KeplerianState, :sma; to = 7300.0)
    @test semi_major_axis(sat) ≈ 7300.0 rtol = 1e-12
    @test eccentricity(sat) ≈ ecc0 rtol = 1e-10

    element!(sat, EquinoctialState, EarthMJ2000Eq, :h; to = 0.0)
    @test element(sat, EquinoctialState, EarthMJ2000Eq, :h) ≈ 0.0 atol = 1e-12

    # A solver writes through the Calc
    q = Calc(element, sat, DelaunayState, :l)
    set_calc!(q, 1.25)
    @test q() ≈ 1.25 atol = 1e-10

    # A Cartesian field takes a vector
    set_quantity!(sat, element, CartesianState, :position; to = [7100.0, 0.0, 100.0])
    @test position_vector(sat) ≈ [7100.0, 0.0, 100.0]

    before = position_vector(sat)
    @test_throws ArgumentError element!(sat, KeplerianState, :ecc; to = -0.5)
    @test position_vector(sat) ≈ before
end

nothing
