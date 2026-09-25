# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A solver varying an orbital element through the quantity interface.
#
# `Vary(semi_major_axis, sat)` writes through `semi_major_axis!`, which holds the other five
# Keplerian elements. Before the setters existed only the legacy `OrbitCalc(sat, SMA())` could do
# this, so the two must now reach the same answer.
#
# Truth: **analytic**. Two-body motion, starting just past periapsis, stopped at the first
# apoapsis. With the eccentricity held, the apoapsis radius is a(1 + e), so the semi-major axis
# that puts apoapsis at 10000 km is 10000 / (1 + e). Both formulations must land there.

using Test
using Epicycle

const _VE_MU   = 398600.4418
const _VE_ECC  = 0.1
const _VE_RA   = 10000.0
const _VE_T0   = Time("2024-01-01T12:00:00", UTC(), ISOT())

_ve_sat() = Spacecraft(
    state = CartesianState(KeplerianState(7500.0, _VE_ECC, 0.5, 0.3, 0.2, 0.1), _VE_MU),
    time  = _VE_T0, coord_sys = EarthMJ2000Eq)

"""Target apoapsis radius by varying the initial semi-major axis, given as `vary(sat)`."""
function _ve_solve(vary)
    sat  = _ve_sat()
    prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                           IntegratorConfig(Tsit5(); dt = 10.0, reltol = 1e-11, abstol = 1e-11))

    var   = vary(sat)
    set   = Event(name = "Set the orbit", event = () -> nothing, vars = [var], funcs = [])
    apo   = Event(name = "Propagate to apoapsis",
                  event = () -> propagate!(prop, sat, StopAt(position_dot_velocity, sat;
                                                             equals = 0.0, direction = -1)),
                  funcs = [Constraint(position_magnitude, sat; equals = _VE_RA)])
    seq = Sequence()
    add_sequence!(seq, set, apo)
    result = solve_trajectory!(seq)
    return result, first(result.variables)
end

@testset "Vary(semi_major_axis, sat) targets apoapsis, as the legacy OrbitCalc does" begin
    a_exact = _VE_RA / (1 + _VE_ECC)

    r_new, a_new = _ve_solve(sat -> Vary(semi_major_axis, sat; guess = 7500.0,
                                         lower_bound = 7000.0, upper_bound = 12000.0))
    @test r_new.info === :Solve_Succeeded
    @test a_new ≈ a_exact rtol = 1e-7

    r_old, a_old = _ve_solve(sat -> SolverVariable(calc = OrbitCalc(sat, SMA()), name = "sma",
                                                   lower_bound = 7000.0, upper_bound = 12000.0))
    @test r_old.info === :Solve_Succeeded
    @test a_new ≈ a_old rtol = 1e-7
end
