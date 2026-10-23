# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Robustness / input-validation dimension (testing standards §2; coding standards §9) for the
# two-body quantities. Satisfies FR-TWOB-9: a quantity undefined for the orbit fails loudly, and
# its message names the energy or eccentricity.
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testitem "quantities undefined for the orbit are rejected loudly" tags=[:FR_TWOB_9, :Robustness] begin
    using AstroRoutines
    μ = 398600.4418
    closed   = ([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], μ)
    open     = ([7000.0, 300.0, 0.0], [0.0, 12.5, 1.0], μ)
    parabola = ([2.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0)              # zero energy, exactly
    radial   = ([7000.0, 0.0, 0.0], [3.0, 0.0, 0.0])

    message(f, args) = try
        f(args...)
        ""
    catch err
        sprint(showerror, err)
    end

    @test_throws ArgumentError orbital_period(open...)
    @test occursin("specific energy", message(orbital_period, open))
    @test_throws ArgumentError apoapsis_speed(open...)
    @test occursin("eccentricity", message(apoapsis_speed, open))
    @test_throws ArgumentError hyperbolic_excess_speed(closed...)
    @test occursin("specific energy", message(hyperbolic_excess_speed, closed))
    @test_throws ArgumentError mean_motion(parabola...)
    @test_throws ArgumentError orbit_normal(radial...)
    @test occursin("radial", message(orbit_normal, radial))
end

@testitem "state vectors of the wrong length are rejected" tags=[:FR_TWOB_9, :Robustness] begin
    using AstroRoutines
    @test_throws ArgumentError specific_energy([7000.0, 0.0], [0.0, 7.5, 0.0], 398600.4418)
    @test_throws ArgumentError angular_momentum([7000.0, 0.0, 0.0], [0.0, 7.5])
end
