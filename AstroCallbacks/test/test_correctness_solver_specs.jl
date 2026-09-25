# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Writing through a quantity, and the constraint spec.
#
# The property under test is that reading and writing are one interface. A
# settability is decided by the subject and the quantity together, so naming the
# quantity once gets
# both directions and a script never says "set". If that breaks, a solver
# silently stops moving the thing it was asked to move.
# =============================================================================

using Test
using AstroCallbacks
using AstroManeuvers
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: label, is_settable, set_quantity!

const _S_EPOCH = Time(2458849.5, 0.0, :tdb, :jd)

_burn() = ImpulsiveManeuver(axes = VNB(), element1 = 1.0, element2 = 2.0, element3 = 3.0)
_sat()  = Spacecraft(CartesianState([7000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]), _S_EPOCH;
                     coord_sys = CoordinateSystem(earth, ICRF()))

@testset "delta_v reads and writes a maneuver" begin
    man = _burn()
    @test delta_v(man) == [1.0, 2.0, 3.0]

    # Agrees with the shipped tag, which took a spacecraft it never used.
    @test delta_v(man) == get_calc(ManeuverCalc(man, _sat(), DeltaVVector()))

    delta_v!(man; to = [4.0, 5.0, 6.0])
    @test delta_v(man) == [4.0, 5.0, 6.0]
    @test man.element1 == 4.0 && man.element2 == 5.0 && man.element3 == 6.0

    # A wrong length is refused before anything is written.
    @test_throws ArgumentError delta_v!(man; to = [1.0, 2.0])
    @test delta_v(man) == [4.0, 5.0, 6.0]
end

@testset "the subject and the quantity decide settability together" begin
    # This is the whole claim: one name, both directions. Settability is not a
    # property of the quantity alone — `set_quantity!` dispatches on the pair,
    # and `is_settable` reports whether that pair has a method.
    man = _burn()
    @test is_settable(man, delta_v)
    @test !is_settable(man, position_magnitude)

    set_quantity!(man, delta_v; to = [7.0, 8.0, 9.0])
    @test delta_v(man) == [7.0, 8.0, 9.0]

    # Typed on the maneuver: a spacecraft has no delta_v to write.
    @test !is_settable(_sat(), delta_v)
end

@testset "set_calc! writes through a Calc" begin
    man = _burn()
    c   = Calc(delta_v, man)

    @test c() == [1.0, 2.0, 3.0]
    @test calc_numvars(c) == 3
    @test calc_is_settable(c)

    set_calc!(c, [7.0, 8.0, 9.0])
    @test c() == [7.0, 8.0, 9.0]
    @test man.element1 == 7.0          # it really reached the maneuver

    # A pair with no set_quantity! method refuses rather than failing later.
    e = try; set_calc!(Calc(position_dot_velocity, _sat()), 1.0); nothing; catch e; e; end
    @test e isa ArgumentError
    @test occursin("cannot be set", e.msg)
end

@testset "a scalar quantity is one variable" begin
    sat = _sat()
    @test calc_numvars(Calc(position_magnitude, sat)) == 1
    @test calc_numvars(Calc(position_vector, sat))    == 3
end

