# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Correctness dimension (TestingStandards.md §2) for the anomaly-conversion functional area.
#
# Truth source: tier-1 analytic. `eccentric_to_mean_anomaly` (M = E − e·sin E) is exact, so it is
# the reference the iterative `mean_to_eccentric_anomaly` is checked against; E↔ν closed forms are
# checked by round-trip (Principles.md — "name the truth source").
#
# Each @testitem runs in its own isolated module, so every block is self-contained.

@testitem "mean ↔ eccentric round trip and Kepler residual" tags=[:FR_ANOM_1, :FR_ANOM_2, :PR_ANOM_1, :Correctness] begin
    using AstroRoutines
    for e in (0.0, 0.01, 0.2, 0.5, 0.7, 0.9, 0.99)
        for E in range(-π, π; length = 33)
            M  = eccentric_to_mean_anomaly(E, e)             # exact truth
            E2 = mean_to_eccentric_anomaly(M, e)             # iterative solve
            @test isapprox(E2, E; atol = 1e-12, rtol = 0)    # round trip (FR-ANOM-1/2)
            @test abs(E2 - e * sin(E2) - M) < 1e-12          # Kepler residual (PR-ANOM-1)
        end
    end
end

@testitem "eccentric ↔ true round trip" tags=[:FR_ANOM_3, :FR_ANOM_4, :Correctness] begin
    using AstroRoutines
    for e in (0.0, 0.2, 0.5, 0.9, 0.99)
        for E in range(-3.1, 3.1; length = 33)              # within (−π, π]
            ν  = eccentric_to_true_anomaly(E, e)
            E2 = true_to_eccentric_anomaly(ν, e)
            @test isapprox(E2, E; atol = 1e-12, rtol = 0)
        end
    end
end

@testitem "mean ↔ true round trip" tags=[:FR_ANOM_5, :FR_ANOM_6, :Correctness] begin
    using AstroRoutines
    for e in (0.0, 0.2, 0.5, 0.9, 0.99)
        for M in range(-3.1, 3.1; length = 33)
            ν  = mean_to_true_anomaly(M, e)
            M2 = true_to_mean_anomaly(ν, e)
            @test isapprox(M2, M; atol = 1e-11, rtol = 0)
        end
    end
end

@testitem "circular limit e = 0: M = E = ν" tags=[:FR_ANOM_1, :FR_ANOM_3, :Correctness] begin
    using AstroRoutines
    for θ in range(-3.0, 3.0; length = 13)
        @test isapprox(mean_to_eccentric_anomaly(θ, 0.0), θ; atol = 1e-13, rtol = 0)
        @test isapprox(eccentric_to_true_anomaly(θ, 0.0), θ; atol = 1e-13, rtol = 0)
        @test isapprox(mean_to_true_anomaly(θ, 0.0), θ; atol = 1e-13, rtol = 0)
    end
end
