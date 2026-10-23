# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Robustness / input-validation dimension (testing standards §2; coding standards §9) for the
# anomaly-conversion functional area. Satisfies FR-ANOM-7: non-elliptic eccentricity fails loudly.
#
# The @testitem runs in its own isolated module and is self-contained.

@testitem "non-elliptic eccentricity is rejected loudly" tags=[:FR_ANOM_7, :Robustness] begin
    using AstroRoutines
    for f in (mean_to_eccentric_anomaly, eccentric_to_mean_anomaly,
              eccentric_to_true_anomaly, true_to_eccentric_anomaly,
              mean_to_true_anomaly, true_to_mean_anomaly)
        @test_throws ArgumentError f(1.0, 1.0)             # e = 1 (parabolic boundary)
        @test_throws ArgumentError f(1.0, 1.5)             # e > 1 (hyperbolic)
        @test_throws ArgumentError f(1.0, -0.1)            # e < 0
    end

    # The message names the offending argument (coding standards §9.4).
    msg = try
        mean_to_eccentric_anomaly(1.0, 1.0)
    catch err
        sprint(showerror, err)
    end
    @test occursin("e", msg)
    @test occursin("[0, 1)", msg)
end

@testitem "non-hyperbolic eccentricity and anomalies beyond the asymptotes are rejected loudly" tags=[:FR_ANOM_12, :Robustness] begin
    using AstroRoutines
    for f in (true_to_hyperbolic_anomaly, hyperbolic_to_true_anomaly,
              hyperbolic_to_mean_anomaly, mean_to_hyperbolic_anomaly)
        @test_throws ArgumentError f(0.5, 1.0)             # e = 1 (parabolic boundary)
        @test_throws ArgumentError f(0.5, 0.5)             # e < 1 (elliptic)
    end
    msg = try
        hyperbolic_to_mean_anomaly(1.0, 0.5)
    catch err
        sprint(showerror, err)
    end
    @test occursin("greater than 1", msg)

    # cos ν∞ = −1/e: at and beyond the asymptote the true anomaly has no hyperbolic anomaly.
    e = 2.0
    ν∞ = acos(-1 / e)
    @test_throws ArgumentError true_to_hyperbolic_anomaly(ν∞ + 1e-6, e)
    @test_throws ArgumentError true_to_hyperbolic_anomaly(π, e)
    msg = try
        true_to_hyperbolic_anomaly(π, e)
    catch err
        sprint(showerror, err)
    end
    @test occursin("asymptotes", msg)
end

@testitem "the hyperbolic Kepler solve converges for large mean anomaly and near-parabolic e" tags=[:PR_ANOM_2, :Robustness] begin
    using AstroRoutines
    for e in (1.000001, 1.0001, 1.01, 2.0, 1e4), M in (-1e5, -10.0, -0.02, 1e-9, 0.1, 10.0, 1e5)
        H = mean_to_hyperbolic_anomaly(M, e)
        @test abs(e * sinh(H) - H - M) < 1e-10 * max(1, abs(M))
    end
end
