# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# `Vary` — the spec that lets a solver change a quantity.
#
# Two properties matter. A `Vary` must round-trip through the same path the
# solver uses, `get_sol_var` and `set_sol_var`, or the optimizer will read one
# thing and write another. And a quantity with no setter must be refused when
# the spec is built, not when the solve is halfway through.
# =============================================================================

using Test
using AstroSolve
using AstroSolve: get_sol_var, set_sol_var
using AstroCallbacks
using AstroManeuvers
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse

_v_burn() = ImpulsiveManeuver(axes = VNB(), element1 = 1.0, element2 = 2.0, element3 = 3.0)

@testset "Vary round-trips through the solver's own path" begin
    man = _v_burn()
    v   = Vary(delta_v, man; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 8.0, 8.0],
                             name = "dv")

    @test v.numvars == 3
    @test v.name == "dv"
    @test v.lower_bound == [0.0, 0.0, 0.0]
    @test v.upper_bound == [8.0, 8.0, 8.0]

    @test get_sol_var(v) == [1.0, 2.0, 3.0]

    set_sol_var(v, [4.0, 5.0, 6.0])
    @test get_sol_var(v) == [4.0, 5.0, 6.0]
    @test delta_v(man)   == [4.0, 5.0, 6.0]     # it reached the maneuver itself
end

@testset "the initial guess lives in the spec" begin
    man = _v_burn()
    Vary(delta_v, man; guess = [0.5, 0.6, 0.7])
    @test delta_v(man) == [0.5, 0.6, 0.7]       # written when the spec was built
end

@testset "the subject is the maneuver, not the spacecraft" begin
    # The old ManeuverCalc took both; the spacecraft was never used to evaluate
    # a delta-V. Vary names one subject, in the slot every spec uses.
    man = _v_burn()
    v   = Vary(delta_v, man)
    @test AstroCallbacks._subjects_from_calc(v.calc) === (man,)
end

@testset "a quantity with no setter is refused at build time" begin
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                     Time(2458849.5, 0.0, :tdb, :jd);
                     coord_sys = CoordinateSystem(earth, ICRF()))
    e = try; Vary(position_dot_velocity, sat); nothing; catch e; e; end
    @test e isa ArgumentError
    # Assert the meaning, not the mechanism: this used to check for the word
    # "setter", which tied the test to the trait table that set_quantity!
    # replaced.
    @test occursin("cannot be set", e.msg)
    @test occursin("Spacecraft", e.msg)      # the error names the pair
end

@testset "Vary and SolverVariable describe the same thing" begin
    # Vary is the quantity spelling of what SolverVariable already was, so the
    # two must produce the same numbers. When the tag stack goes, this goes.
    man_a, man_b = _v_burn(), _v_burn()
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                     Time(2458849.5, 0.0, :tdb, :jd);
                     coord_sys = CoordinateSystem(earth, ICRF()))

    new_way = Vary(delta_v, man_a; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 8.0, 8.0])
    old_way = SolverVariable(calc = ManeuverCalc(man_b, sat, DeltaVVector()),
                             lower_bound = [0.0, 0.0, 0.0],
                             upper_bound = [8.0, 8.0, 8.0])

    @test new_way.numvars     == old_way.numvars
    @test get_sol_var(new_way) == get_sol_var(old_way)

    set_sol_var(new_way, [1.5, 2.5, 3.5])
    set_sol_var(old_way, [1.5, 2.5, 3.5])
    @test get_sol_var(new_way) == get_sol_var(old_way)
end
