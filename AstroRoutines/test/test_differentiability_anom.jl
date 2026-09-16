# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Differentiability dimension (TestingStandards.md §3) for the anomaly-conversion functional area.
# Validation tier "AD vs analytic": ForwardDiff partials are checked against the hand-derived
# analytic derivatives, and against finite differences for the closed forms. Satisfies SR-AD-1.
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testitem "Kepler solve is differentiable (dE/dM, dE/de vs analytic)" tags=[:SR_AD_1, :Differentiability] begin
    using AstroRoutines, ForwardDiff
    for e in (0.0, 0.2, 0.7, 0.95)
        for M in range(-3.0, 3.0; length = 9)
            E = mean_to_eccentric_anomaly(M, e)
            dEdM = ForwardDiff.derivative(m -> mean_to_eccentric_anomaly(m, e), M)
            dEde = ForwardDiff.derivative(ee -> mean_to_eccentric_anomaly(M, ee), e)
            @test isapprox(dEdM, 1 / (1 - e * cos(E));            rtol = 1e-10)
            @test isapprox(dEde, sin(E) / (1 - e * cos(E));       rtol = 1e-10, atol = 1e-12)
        end
    end
end

@testitem "closed-form conversions are differentiable (AD vs finite difference)" tags=[:SR_AD_1, :Differentiability] begin
    using AstroRoutines, ForwardDiff, FiniteDiff
    e = 0.3
    for x in range(-2.5, 2.5; length = 9)
        for f in (E -> eccentric_to_true_anomaly(E, e),
                  ν -> true_to_eccentric_anomaly(ν, e),
                  E -> eccentric_to_mean_anomaly(E, e))
            ad = ForwardDiff.derivative(f, x)
            fd = FiniteDiff.finite_difference_derivative(f, x)
            @test isapprox(ad, fd; rtol = 1e-6, atol = 1e-8)
        end
    end
end

@testitem "AD propagates through the composed mean→true path" tags=[:SR_AD_1, :Differentiability] begin
    using AstroRoutines, ForwardDiff
    g = m -> mean_to_true_anomaly(m, 0.4)
    d = ForwardDiff.derivative(g, 1.0)
    @test isfinite(d)
    @test d > 0                                             # ν increases with M on (−π, π)
end
