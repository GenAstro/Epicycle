# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Correctness dimension (TestingStandards.md §2) for Lagrange and Hermite interpolation
# (ticket #9). Truth: polynomials, which an interpolant of sufficient degree reproduces exactly,
# and DataInterpolations as an independent implementation where the two define the same window.

@testitem "Lagrange reproduces polynomials up to its degree" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    t = [0.0, 0.7, 1.9, 3.0, 4.2, 5.0, 6.5, 7.1, 8.0, 9.4, 10.0]       # uneven spacing
    for d in 1:9
        p(x) = sum(((-1)^k * 0.3^k / (k + 1)) * x^k for k in 0:d)
        y = p.(t)
        err = maximum(abs(interpolate(Lagrange(d), t, y, s) - p(s)) for s in 0:0.05:10)
        @test err < 1e-10 * max(1, maximum(abs, y))
    end
    # One degree short is not exact: the test above can fail.
    q(x) = x^4
    @test abs(interpolate(Lagrange(3), t, q.(t), 2.5) - q(2.5)) > 1e-3
end

@testitem "Hermite reproduces polynomials up to its degree" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    t = [0.0, 0.7, 1.9, 3.0, 4.2, 5.0, 6.5, 7.1, 8.0]
    for d in (1, 3, 5, 7, 9)
        p(x) = sum((0.5^k / (k + 1)) * x^k for k in 0:d)
        dp(x) = sum((k * 0.5^k / (k + 1)) * x^(k - 1) for k in 1:d)
        y, v = p.(t), dp.(t)
        for s in 0:0.05:8
            val, rate = interpolate_with_derivative(Hermite(d), t, y, v, s)
            @test isapprox(val, p(s); rtol = 1e-11, atol = 1e-11)
            @test isapprox(rate, dp(s); rtol = 1e-10, atol = 1e-10)
        end
    end
end

@testitem "Hermite degree 7 needs only four points" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    @test required_points(Hermite(7)) == 4
    @test required_points(Lagrange(7)) == 8
    t = [0.0, 1.0, 2.5, 4.0]
    p(x) = 1 + x - x^3 / 7 + x^7 / 5000
    dp(x) = 1 - 3x^2 / 7 + 7x^6 / 5000
    @test isapprox(interpolate(Hermite(7), t, p.(t), dp.(t), 3.3), p(3.3); rtol = 1e-13)
end

@testitem "tabulated points are returned exactly" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    t = collect(0.0:1.0:10.0)
    y = sin.(t) .* 1e4
    v = cos.(t) .* 1e4
    for i in eachindex(t)
        @test interpolate(Lagrange(7), t, y, t[i]) === y[i]
        @test interpolate(Hermite(5), t, y, v, t[i]) === y[i]
        @test interpolate_with_derivative(Hermite(5), t, y, v, t[i])[2] === v[i]
    end
end

@testitem "window is centred, and shifts inward at the ends" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    # A step function shows which points the window uses: the value is exact only if all
    # used points lie on the same side of the step.
    t = collect(0.0:1.0:20.0)
    y = [x < 10.5 ? 0.0 : 1.0 for x in t]
    @test interpolate(Lagrange(3), t, y, 2.5) == 0.0          # window 1:4, all zero
    @test interpolate(Lagrange(3), t, y, 0.2) == 0.0          # shifted inward at the start
    @test interpolate(Lagrange(3), t, y, 19.8) == 1.0         # and at the end
    @test interpolate(Lagrange(3), t, y, 8.5) == 0.0          # points 8..11 = t 7..10: all zero
    @test interpolate(Lagrange(3), t, y, 9.5) != 0.0          # points 9..12 cross the step
    # Odd count: the extra point is taken on the side nearer t.
    @test interpolate(Lagrange(2), t, y, 9.2) == 0.0          # nearer 9: points 8,9,10
    @test interpolate(Lagrange(2), t, y, 9.8) != 0.0          # nearer 10: points 9,10,11 straddle
end

@testitem "vector values: positions and velocities" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines, StaticArrays
    t = collect(0.0:60.0:3600.0)
    n = 2π / 5400
    r = [SVector(7000cos(n * s), 7000sin(n * s), 0.0) for s in t]
    v = [SVector(-7000n * sin(n * s), 7000n * cos(n * s), 0.0) for s in t]
    s = 1234.5
    truth_r = SVector(7000cos(n * s), 7000sin(n * s), 0.0)
    truth_v = SVector(-7000n * sin(n * s), 7000n * cos(n * s), 0.0)
    rl = interpolate(Lagrange(7), t, r, s)
    rh, vh = interpolate_with_derivative(Hermite(7), t, r, v, s)
    @test rl isa SVector{3, Float64}
    # Lagrange 7 at 60 s steps: the error bound h⁸·max|r⁽⁸⁾|·Π|s − tₖ|/8!/h⁸ is about 4e-9 km.
    @test maximum(abs.(rl - truth_r)) < 1e-8          # km
    @test maximum(abs.(rh - truth_r)) < 1e-9
    @test maximum(abs.(vh - truth_v)) < 1e-12         # km/s
    # Plain Vector values work too.
    @test interpolate(Lagrange(7), t, Vector.(r), s) ≈ Vector(rl)
end

@testitem "agrees with DataInterpolations" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines, DataInterpolations
    t = [0.0, 0.4, 1.1, 1.5, 2.6, 3.0, 3.9, 4.4, 5.0]
    y = exp.(sin.(t))
    v = cos.(t) .* y
    # Lagrange: with exactly degree + 1 points, both use the whole table.
    li = LagrangeInterpolation(y, t, length(t) - 1)
    for s in 0:0.1:5
        @test isapprox(interpolate(Lagrange(length(t) - 1), t, y, s), li(s); rtol = 1e-12)
    end
    # Hermite degree 3 is the cubic Hermite spline: both use the two points bracketing s.
    hi = CubicHermiteSpline(v, y, t)
    for s in 0:0.1:5
        @test isapprox(interpolate(Hermite(3), t, y, v, s), hi(s); rtol = 1e-12)
    end
end

@testitem "safe to evaluate from many tasks at once" tags=[:Interpolation, :Correctness] begin
    using AstroRoutines
    t = collect(0.0:30.0:6000.0)
    y = sin.(t ./ 900)
    v = cos.(t ./ 900) ./ 900
    s = collect(range(0.0, 6000.0; length = 4001))
    serial = [interpolate(Hermite(7), t, y, v, x) for x in s]
    tasks = [Threads.@spawn([interpolate(Hermite(7), t, y, v, x) for x in s]) for _ in 1:16]
    @test all(fetch(task) == serial for task in tasks)
    lag = [interpolate(Lagrange(9), t, y, x) for x in s]
    tasks = [Threads.@spawn([interpolate(Lagrange(9), t, y, x) for x in s]) for _ in 1:16]
    @test all(fetch(task) == lag for task in tasks)
end
