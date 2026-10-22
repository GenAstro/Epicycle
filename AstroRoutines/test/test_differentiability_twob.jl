# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Differentiability dimension (TestingStandards.md §3) for the two-body quantities. Validation
# tier "AD vs finite difference": ForwardDiff gradients and Jacobians with respect to the state and
# μ are checked against FiniteDiff. Satisfies SR-AD-1 for this area.
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testitem "two-body quantities are differentiable (AD vs finite difference)" tags=[:SR_AD_1, :Differentiability] begin
    using AstroRoutines, ForwardDiff, FiniteDiff
    μ = 398600.4418
    for x in ([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9],          # closed
              [7000.0, 300.0, 0.0, 0.0, 12.5, 1.0])              # hyperbolic
        y = [x; μ]
        closed = specific_energy(x[1:3], x[4:6], μ) < 0
        scalars = Function[z -> specific_energy(z[1:3], z[4:6], z[7]),
                           z -> angular_momentum_magnitude(z[1:3], z[4:6]),
                           z -> mean_motion(z[1:3], z[4:6], z[7]),
                           z -> periapsis_speed(z[1:3], z[4:6], z[7])]
        if closed
            push!(scalars, z -> orbital_period(z[1:3], z[4:6], z[7]),
                           z -> apoapsis_speed(z[1:3], z[4:6], z[7]))
        else
            push!(scalars, z -> hyperbolic_excess_speed(z[1:3], z[4:6], z[7]))
        end
        for f in scalars
            @test isapprox(ForwardDiff.gradient(f, y), FiniteDiff.finite_difference_gradient(f, y);
                           rtol = 1e-6, atol = 1e-10)
        end
        for f in (z -> angular_momentum(z[1:3], z[4:6]),
                  z -> orbit_normal(z[1:3], z[4:6]),
                  z -> orbital_angular_velocity(z[1:3], z[4:6]),
                  z -> eccentricity_vector(z[1:3], z[4:6], z[7]))
            @test isapprox(ForwardDiff.jacobian(f, y), FiniteDiff.finite_difference_jacobian(f, y);
                           rtol = 1e-6, atol = 1e-10)
        end
    end
end
