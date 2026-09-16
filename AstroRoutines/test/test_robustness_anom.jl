# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Robustness / input-validation dimension (TestingStandards.md §2; CodingStandards.md §9) for the
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

    # The message names the offending argument (CodingStandards §9.4).
    msg = try
        mean_to_eccentric_anomaly(1.0, 1.0)
    catch err
        sprint(showerror, err)
    end
    @test occursin("e", msg)
    @test occursin("[0, 1)", msg)
end
