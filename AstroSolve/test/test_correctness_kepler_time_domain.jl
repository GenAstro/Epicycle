# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Time-domain Kepler propagation over zero and short steps (issue #26).
#
# A zero step returned NaN on an ellipse: the closed-form Stumpff functions are 0/0 at y = 0, and
# lose digits by cancellation for small y on both conic branches. Truth here is:
#
#   · **analytic** for a zero step: the initial state, the identity STM and [v; a];
#   · a **Taylor step** in position and velocity for steps of 1e-9 and 1e-6 s;
#   · **central differences** of the propagated state for the STM and the time partial;
#   · **composition** across the series limit: one step against two half steps, with the half
#     steps on the series side and the full step on the closed-form side.
#
# Covers the series and zero-step paths of kepler_time_domain.jl.

using LinearAlgebra
using AstroSolve
using ForwardDiff
using Test

const _KTD_MU = 398600.4418           # Earth, km³/s²

# Ellipse, hyperbola, and orbits either side of the parabola
const _KTD_CASES = (
    ellipse            = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0],
    inclined_ellipse   = [6500.0, 1200.0, -800.0, -1.1, 7.2, 2.5],
    hyperbola          = [7000.0, 0.0, 0.0, 0.0, 12.0, 0.0],
    parabola           = [7000.0, 0.0, 0.0, 0.0, sqrt(2 * 398600.4418 / 7000.0), 0.0],
    near_parabola_hyp  = [7000.0, 0.0, 0.0, 0.0, sqrt(2 * 398600.4418 / 7000.0) * (1 + 1e-6), 0.0],
    near_parabola_ell  = [7000.0, 0.0, 0.0, 0.0, sqrt(2 * 398600.4418 / 7000.0) * (1 - 1e-6), 0.0],
)

_ktd_accel(x) = -_KTD_MU / norm(x[1:3])^3 * x[1:3]

# Two-body state after a short step: position to third order, velocity to second
function _ktd_taylor(x, dt)
    r = x[1:3]; v = x[4:6]; rn = norm(r)
    a = _ktd_accel(x)
    j = -_KTD_MU / rn^3 * (v - 3 * dot(r, v) / rn^2 * r)
    [r + v * dt + a * dt^2 / 2 + j * dt^3 / 6; v + a * dt + j * dt^2 / 2]
end

_ktd_prop(x, dt) = AstroSolve.kepler_propagate_time_domain(x, _KTD_MU, dt).state

# STM by central differences. Over a microsecond the map is linear to far below roundoff, so a
# large step keeps the rounding of a 7000 km position out of the quotient.
function _ktd_stm_fd(x, dt)
    Φ = zeros(6, 6)
    for j in 1:6
        h = 0.1
        e = zeros(6); e[j] = h
        Φ[:, j] = (_ktd_prop(x + e, dt) - _ktd_prop(x - e, dt)) / 2h
    end
    Φ
end

@testset "kepler_time_domain: zero and short steps" begin

    @testset "zero step: $name" for (name, x0) in pairs(_KTD_CASES)
        res = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, 0.0)
        @test res.converged
        @test res.state == x0
        @test res.stm == Matrix(1.0I, 6, 6)
        @test res.dstate_dt ≈ [x0[4:6]; _ktd_accel(x0)] rtol = 1e-14
        @test (res.F, res.G, res.Ft, res.Gt) == (1.0, 0.0, 0.0, 1.0)

        bare = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, 0.0;
                                                       need_stm = false, need_time_partials = false)
        @test bare.state == x0
        @test bare.stm === nothing && bare.dstate_dt === nothing
    end

    @testset "short step against Taylor: $name, dt = $dt" for (name, x0) in pairs(_KTD_CASES),
                                                              dt in (1e-9, -1e-9, 1e-6, -1e-6)
        res = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, dt)
        @test res.converged
        @test all(isfinite, res.state) && all(isfinite, res.stm) && all(isfinite, res.dstate_dt)
        # The solve stops on an absolute universal-anomaly tolerance (Xtol = 1e-12), which is a
        # timing error of up to r0·Xtol/√μ ≈ 1e-11 s, so about 1e-10 km at orbital speed.
        xt = _ktd_taylor(x0, dt)
        @test norm(res.state[1:3] - xt[1:3]) < 2e-10          # km
        @test norm(res.state[4:6] - xt[4:6]) < 2e-13          # km/s
        @test res.dstate_dt ≈ [res.state[4:6]; _ktd_accel(res.state)] rtol = 1e-10
        @test res.stm ≈ _ktd_stm_fd(x0, dt) atol = 1e-9
    end

    @testset "time partial against central differences: $name" for (name, x0) in pairs(_KTD_CASES)
        dt = 1e-6; h = 1e-7
        res = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, dt)
        fd = (_ktd_prop(x0, dt + h) - _ktd_prop(x0, dt - h)) / 2h
        @test res.dstate_dt ≈ fd rtol = 1e-6
    end

    # Elliptic and hyperbolic steps whose y = αX² falls either side of the series limit
    @testset "continuity across the series limit: $name" for name in (:ellipse, :inclined_ellipse, :hyperbola)
        x0 = _KTD_CASES[name]
        for dt in (150.0, 300.0, 600.0, -300.0)
            full = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, dt)
            half = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, dt / 2)
            two  = AstroSolve.kepler_propagate_time_domain(half.state, _KTD_MU, dt / 2)
            @test norm(full.state[1:3] - two.state[1:3]) < 1e-8      # km
            @test norm(full.state[4:6] - two.state[4:6]) < 1e-11     # km/s
            @test full.stm ≈ two.stm * half.stm rtol = 1e-9
        end
    end

    # A dual dt at zero carries a time partial. The rates of the state, the STM and the time
    # partial are checked against [v; a] and central differences of the outputs over ±1 ms.
    @testset "AD through a zero step: $name" for (name, x0) in pairs(_KTD_CASES)
        prop(t) = AstroSolve.kepler_propagate_time_domain(x0, _KTD_MU, t)
        @test ForwardDiff.derivative(t -> prop(t).state, 0.0) ≈ [x0[4:6]; _ktd_accel(x0)] rtol = 1e-12
        h = 1e-3
        @test ForwardDiff.derivative(t -> prop(t).stm, 0.0) ≈ (prop(h).stm - prop(-h).stm) / 2h atol = 1e-9
        @test ForwardDiff.derivative(t -> prop(t).dstate_dt, 0.0) ≈
              (prop(h).dstate_dt - prop(-h).dstate_dt) / 2h rtol = 1e-6
    end
end
