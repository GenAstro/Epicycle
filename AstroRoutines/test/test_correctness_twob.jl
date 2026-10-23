# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Correctness dimension (testing standards §2) for the two-body quantities (FR-TWOB-1…8).
#
# Truth source: tier-1 analytic. Each state is built from classical elements (a, e, i, Ω, ω, ν)
# with the perifocal closed form, and each quantity is checked against its expression in those
# elements, which the routines never see. The same quantities are checked against GMAT through
# the Spacecraft quantities in AstroCallbacks (test_correctness_orbit_parameters.jl).
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testsnippet TwoBodyStates begin
    using AstroRoutines
    const μ = 398600.4418

    cross_(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]

    # Position and velocity from classical elements, through the perifocal frame; also the unit
    # vectors to periapsis (P̂) and along the orbit normal (Ŵ).
    function state_from_elements(a, e, i, Ω, ω, ν)
        p = a * (1 - e^2)
        P̂ = [cos(Ω) * cos(ω) - sin(Ω) * sin(ω) * cos(i),
             sin(Ω) * cos(ω) + cos(Ω) * sin(ω) * cos(i),
             sin(ω) * sin(i)]
        Q̂ = [-cos(Ω) * sin(ω) - sin(Ω) * cos(ω) * cos(i),
             -sin(Ω) * sin(ω) + cos(Ω) * cos(ω) * cos(i),
             cos(ω) * sin(i)]
        r = p / (1 + e * cos(ν)) .* (cos(ν) .* P̂ .+ sin(ν) .* Q̂)
        v = sqrt(μ / p) .* (-sin(ν) .* P̂ .+ (e + cos(ν)) .* Q̂)
        return r, v, P̂, cross_(P̂, Q̂)
    end

    # Near-circular, eccentric, polar, retrograde; and two hyperbolas, inbound and outbound.
    const CLOSED = ((7000.0, 0.001, 0.5, 1.0, 2.0, 0.3), (26600.0, 0.74, 1.1, 0.7, 4.7, 2.5),
                    (7200.0, 0.05, π / 2, 2.1, 0.8, 4.0), (9000.0, 0.2, 2.8, 5.2, 1.4, 1.3))
    const OPEN   = ((-20000.0, 1.4, 0.5, 1.0, 0.3, -0.9), (-8000.0, 2.5, 1.2, 4.0, 2.0, 1.1))
end

@testitem "energy, angular momentum and eccentricity vector follow the elements" tags=[:FR_TWOB_1, :FR_TWOB_2, :FR_TWOB_5, :Correctness] setup=[TwoBodyStates] begin
    for (a, e, i, Ω, ω, ν) in (CLOSED..., OPEN...)
        r, v, P̂, Ŵ = state_from_elements(a, e, i, Ω, ω, ν)
        p = a * (1 - e^2)
        @test specific_energy(r, v, μ) ≈ -μ / (2a) rtol = 1e-12
        @test angular_momentum_magnitude(r, v) ≈ sqrt(μ * p) rtol = 1e-12
        @test angular_momentum(r, v) ≈ sqrt(μ * p) .* Ŵ rtol = 1e-12
        @test eccentricity_vector(r, v, μ) ≈ e .* P̂ rtol = 1e-10 atol = 1e-13
    end
end

@testitem "orbit normal and angular velocity follow the plane and the radius" tags=[:FR_TWOB_3, :FR_TWOB_4, :Correctness] setup=[TwoBodyStates] begin
    for (a, e, i, Ω, ω, ν) in (CLOSED..., OPEN...)
        r, v, _, Ŵ = state_from_elements(a, e, i, Ω, ω, ν)
        n̂ = orbit_normal(r, v)
        @test n̂ ≈ Ŵ rtol = 1e-12
        @test abs(sum(n̂ .* r)) < 1e-9
        @test abs(sum(n̂ .* v)) < 1e-12
        rate = orbital_angular_velocity(r, v)
        @test rate ≈ sqrt(μ * a * (1 - e^2)) / sum(abs2, r) .* Ŵ rtol = 1e-12     # ν̇ = h/r²
    end
end

@testitem "period, mean motion and characteristic speeds follow the elements" tags=[:FR_TWOB_6, :FR_TWOB_7, :FR_TWOB_8, :Correctness] setup=[TwoBodyStates] begin
    for (a, e, i, Ω, ω, ν) in CLOSED
        r, v = state_from_elements(a, e, i, Ω, ω, ν)
        p = a * (1 - e^2)
        @test orbital_period(r, v, μ) ≈ 2π * sqrt(a^3 / μ) rtol = 1e-12
        @test mean_motion(r, v, μ) ≈ sqrt(μ / a^3) rtol = 1e-12
        @test apoapsis_speed(r, v, μ) ≈ sqrt(μ / p) * (1 - e) rtol = 1e-10
        @test periapsis_speed(r, v, μ) ≈ sqrt(μ / p) * (1 + e) rtol = 1e-12
        # vis-viva at the apsides, an independent check of both speeds
        @test apoapsis_speed(r, v, μ)^2 ≈ μ * (2 / (a * (1 + e)) - 1 / a) rtol = 1e-10
        @test periapsis_speed(r, v, μ)^2 ≈ μ * (2 / (a * (1 - e)) - 1 / a) rtol = 1e-12
    end
    for (a, e, i, Ω, ω, ν) in OPEN
        r, v = state_from_elements(a, e, i, Ω, ω, ν)
        @test hyperbolic_excess_speed(r, v, μ) ≈ sqrt(-μ / a) rtol = 1e-12
        @test mean_motion(r, v, μ) ≈ sqrt(μ / (-a)^3) rtol = 1e-12
        @test periapsis_speed(r, v, μ)^2 ≈ μ * (2 / (a * (1 - e)) - 1 / a) rtol = 1e-12
    end
end

@testitem "circular orbit: e = 0, and both apsis speeds are the circular speed" tags=[:FR_TWOB_5, :FR_TWOB_8, :Correctness] setup=[TwoBodyStates] begin
    r, v = state_from_elements(7000.0, 0.0, 0.4, 1.0, 0.0, 2.0)
    vc = sqrt(μ / 7000.0)
    @test sqrt(sum(abs2, eccentricity_vector(r, v, μ))) < 1e-15
    @test apoapsis_speed(r, v, μ) ≈ vc rtol = 1e-13
    @test periapsis_speed(r, v, μ) ≈ vc rtol = 1e-13
end
