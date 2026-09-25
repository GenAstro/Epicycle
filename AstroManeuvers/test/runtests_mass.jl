# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A burn draws propellant from the spacecraft, and with no tanks modeled that is the
# whole mass pool. The GMAT-referenced mass values in runtests.jl are the check that this change moved no
# numbers; these cover the behaviour around the edges.

using Test
using Logging

using EpicycleBase
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using AstroModels
using AstroManeuvers

_burn_sc(; mass = 1000.0) = Spacecraft(
    state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
    time  = Time("2015-09-21T00:00:00", TAI(), ISOT()),
    mass  = mass,
    name  = "SC-burn")

@testset "a lumped-mass burn says so once" begin
    sc = _burn_sc()
    small = ImpulsiveManeuver(axes = Inertial(), Isp = 300.0, element1 = 0.01)

    # First burn explains that the pool includes payload; the second is silent.
    logs = @test_logs (:warn,) match_mode = :any maneuver!(sc, small)
    @test_logs min_level = Logging.Warn maneuver!(sc, small)

    @test total_mass(sc) < 1000.0
end

@testset "an over-sized burn returns rather than throwing" begin
    # A solver varying Δv proposes infeasible burns as a matter of course. Throwing
    # would end the run instead of letting it step back, so this asserts the
    # absence of an exception rather than a particular mass.
    sc = _burn_sc(mass = 10.0)
    huge = ImpulsiveManeuver(axes = Inertial(), Isp = 300.0, element1 = 50.0)

    @test maneuver!(sc, huge) isa Spacecraft
    @test total_mass(sc) > 0.0
end

@testset "an impulsive burn approaches zero mass without crossing it" begin
    # Tsiolkovsky gives m₁ = m₀·exp(−Δv/(g₀·Isp)), so mass decays toward zero and never
    # goes negative. The hazard is not a sign change — it is that drag and SRP divide by
    # a mass heading for zero, so the acceleration diverges while every value stays
    # positive and nothing raises. Positivity belongs in the solve as a constraint.
    m0, dv, isp, g0 = 1000.0, 20.0, 300.0, 9.81
    sc  = _burn_sc(mass = m0)
    big = ImpulsiveManeuver(axes = Inertial(), g0 = g0, Isp = isp, element1 = dv)

    maneuver!(sc, big)
    m = total_mass(sc)

    # Truth source: the rocket equation itself. m₁ = m₀·exp(−Δv/(g₀·Isp)), which is
    # strictly positive for every finite Δv.
    @test m ≈ m0 * exp(-dv / (g0 / 1000 * isp)) rtol = 1e-12
    @test m > 0.0
    @test m < m0 / 100                        # three orders down and still not negative
end

@testset "the non-positive notice fires once and is separate from the lumped one" begin
    # Reaching zero exactly needs the exponential to underflow, which takes a Δv near
    # 750·g₀·Isp. Absurd as a trajectory, but it is the only way an impulsive burn
    # reaches the guard, and the guard has to hold for the finite-thrust mass ODE later.
    sc = _burn_sc(mass = 10.0)
    absurd = ImpulsiveManeuver(axes = Inertial(), Isp = 300.0, element1 = 5000.0)

    @test_logs (:warn,) (:warn,) match_mode = :any maneuver!(sc, absurd)
    @test total_mass(sc) == 0.0
    @test_logs min_level = Logging.Warn maneuver!(sc, absurd)

    @test :lumped_mass_burn in getfield(sc, :notified)
    @test :non_positive_mass in getfield(sc, :notified)
end

@testset "a feasible burn raises no mass notice" begin
    sc = _burn_sc()
    small = ImpulsiveManeuver(axes = Inertial(), Isp = 300.0, element1 = 0.01)
    maneuver!(sc, small)

    @test :lumped_mass_burn in getfield(sc, :notified)
    @test !(:non_positive_mass in getfield(sc, :notified))
end
