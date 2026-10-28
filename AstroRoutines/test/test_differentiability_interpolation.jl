# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Differentiability dimension (TestingStandards.md §2) for Lagrange and Hermite interpolation
# (ticket #9): ForwardDiff through the time and through the tabulated values, checked against
# the analytic derivative of the polynomial and against finite differences.

@testitem "derivative with respect to time" tags=[:Interpolation, :Differentiability] begin
    using AstroRoutines, ForwardDiff, FiniteDiff
    t = [0.0, 0.6, 1.3, 2.0, 2.4, 3.1, 4.0, 4.8, 5.5, 6.0]
    y = sin.(t)
    v = cos.(t)
    for s in (0.0, 0.31, 1.3, 2.77, 6.0)                         # includes nodes and the ends
        for m in (Lagrange(5), Lagrange(7))
            ad = ForwardDiff.derivative(x -> interpolate(m, t, y, x), s)
            @test ad ≈ interpolate_with_derivative(m, t, y, s)[2] rtol = 1e-10
            # A finite difference would step outside the table at its ends; there the check
            # above, against the polynomial's own derivative, is the reference.
            if 0.0 < s < 6.0
                fd = FiniteDiff.finite_difference_derivative(x -> interpolate(m, t, y, x), s)
                @test ad ≈ fd rtol = 1e-5
            end
        end
        for m in (Hermite(3), Hermite(7))
            ad = ForwardDiff.derivative(x -> interpolate(m, t, y, v, x), s)
            @test ad ≈ interpolate_with_derivative(m, t, y, v, s)[2] rtol = 1e-10
        end
    end
    # At a node, the derivative of the Hermite interpolant is the tabulated derivative.
    @test ForwardDiff.derivative(x -> interpolate(Hermite(5), t, y, v, x), t[4]) ≈ v[4] rtol = 1e-12
end

@testitem "derivative with respect to the tabulated values" tags=[:Interpolation, :Differentiability] begin
    using AstroRoutines, ForwardDiff, FiniteDiff
    t = collect(0.0:0.5:5.0)
    y0 = cos.(t)
    v0 = -sin.(t)
    s = 2.2
    f(y) = interpolate(Lagrange(7), t, y, s)
    @test ForwardDiff.gradient(f, y0) ≈ FiniteDiff.finite_difference_gradient(f, y0) rtol = 1e-6
    # The interpolant is linear in the values, so the gradient sums to one.
    @test sum(ForwardDiff.gradient(f, y0)) ≈ 1.0 rtol = 1e-12
    g(x) = interpolate(Hermite(5), t, x[1:end÷2], x[end÷2+1:end], s)
    x0 = vcat(y0, v0)
    @test ForwardDiff.gradient(g, x0) ≈ FiniteDiff.finite_difference_gradient(g, x0) rtol = 1e-6
end
