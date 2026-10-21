# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Robustness / input-validation dimension (TestingStandards.md §2; CodingStandards.md §9) for
# Lagrange and Hermite interpolation (ticket #9): bad input fails loudly, naming what is wrong.

@testitem "invalid methods are rejected" tags=[:Interpolation, :Robustness] begin
    using AstroRoutines
    @test_throws ArgumentError Lagrange(0)
    @test_throws ArgumentError Hermite(0)
    @test_throws ArgumentError Hermite(4)                       # even: not reachable from pairs
    msg = try Hermite(4) catch e sprint(showerror, e) end
    @test occursin("odd", msg) && occursin("4", msg)
end

@testitem "no extrapolation" tags=[:Interpolation, :Robustness] begin
    using AstroRoutines
    t = collect(0.0:1.0:10.0)
    y = t .^ 2
    @test_throws ArgumentError interpolate(Lagrange(3), t, y, -0.001)
    @test_throws ArgumentError interpolate(Lagrange(3), t, y, 10.001)
    @test_throws ArgumentError interpolate(Hermite(3), t, y, 2t, 11.0)
    msg = try interpolate(Lagrange(3), t, y, 12.0) catch e sprint(showerror, e) end
    @test occursin("outside", msg) && occursin("10.0", msg)
    @test interpolate(Lagrange(3), t, y, 10.0) == 100.0          # the ends themselves are in
end

@testitem "too few points for the degree" tags=[:Interpolation, :Robustness] begin
    using AstroRoutines
    t = [0.0, 1.0, 2.0]
    @test_throws ArgumentError interpolate(Lagrange(3), t, t, 0.5)          # needs 4
    @test_throws ArgumentError interpolate(Hermite(7), t, t, t, 0.5)        # needs 4
    @test interpolate(Hermite(5), t, t, ones(3), 0.5) ≈ 0.5                 # 3 is enough
    msg = try interpolate(Lagrange(3), t, t, 0.5) catch e sprint(showerror, e) end
    @test occursin("needs 4 points", msg) && occursin("has 3", msg)
end

@testitem "mismatched lengths and unordered times" tags=[:Interpolation, :Robustness] begin
    using AstroRoutines
    t = collect(0.0:1.0:5.0)
    @test_throws ArgumentError interpolate(Lagrange(3), t, ones(5), 1.0)
    @test_throws ArgumentError interpolate(Hermite(3), t, ones(6), ones(5), 1.0)
    bad = [0.0, 1.0, 1.0, 3.0, 4.0, 5.0]                                    # a repeated time
    @test_throws ArgumentError interpolate(Lagrange(3), bad, ones(6), 1.5)
    msg = try interpolate(Lagrange(3), bad, ones(6), 1.5) catch e sprint(showerror, e) end
    @test occursin("increase", msg)
end
