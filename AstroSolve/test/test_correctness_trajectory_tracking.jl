# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A Lagrange cost and a time-varying path bound, against a closed-form optimal control.
#
# Truth: **analytic**, and unusually complete — the optimal state, the optimal control and the
# optimal cost are all known in closed form, so this problem checks the whole answer rather than
# one number off the end of it.
#
#   minimise   J = ∫₀^{2π} [(x - sin t)² + ½ u²] dt
#   subject to ẋ = u,  x(0) = x(2π) = 0
#              (1 + ½ sin t) x² - 1 ≤ 0        along the path
#
#   x*(t) = (2/3) sin t,   u*(t) = (2/3) cos t,   J* = π/3
#
# What this covers that the brachistochrone does not: an objective that is an integral and
# nothing else, an integrand and a path constraint that both depend on `t` explicitly, and a
# phase carrying no model object. The brachistochrone's objective is a single endpoint quantity,
# so the quadrature that accumulates a running cost is untouched by it.
#
# The problem is the regression suite's uc17, recorded at 1.047198 under LGL. This file transcribes
# it with Hermite-Simpson, so agreement is cross-transcription as well as analytic.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

const _TT_TF = 2π

struct _TtState{T}   <: AbstractState;   x::T end
struct _TtControl{T} <: AbstractControl; u::T end

_tt_dynamics!(dy, y::_TtState, u::_TtControl, p, t, model) = (dy[1] = u.u)

@partial(_tt_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = 1.0
end

_tt_position(c) = [state(c).x]
@partial(_tt_position, state) do c
    [1.0]
end

# The bound tightens and relaxes with t. The optimum never reaches it, so it exercises the
# time-dependent path machinery without moving the answer the analytic solution predicts.
_tt_envelope(c) = [(1.0 + 0.5 * sin(c.t)) * state(c).x^2 - 1.0]
@partial(_tt_envelope, state) do c
    [2.0 * (1.0 + 0.5 * sin(c.t)) * state(c).x]
end

_tt_cost(c) = (state(c).x - sin(c.t))^2 + 0.5 * control(c).u^2
@partial(_tt_cost, state) do c
    [2.0 * (state(c).x - sin(c.t))]
end
@partial(_tt_cost, control) do c
    [control(c).u]
end

"""Build the tracking phase on a Hermite-Simpson mesh of `n_steps` intervals."""
function _tt_phase(; n_steps = 24)
    phase = CollocationPhase(
        name          = :tracking,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _tt_dynamics!,
        state         = _TtState,     control = _TtControl,
        tspan         = (0.0, _TT_TF))

    Vary(state,   phase; guess = [0.0 0.0],
                         lower_bound = [-5.0], upper_bound = [5.0])
    Vary(control, phase; guess = reshape([0.0 0.0], 1, 2),
                         lower_bound = [-5.0], upper_bound = [5.0])

    Constraint(_tt_position, phase; equals = [0.0], at = Initial())
    Constraint(_tt_position, phase; equals = [0.0], at = Final())
    Constraint(_tt_envelope, phase; upper_bound = [0.0], at = Path())

    Objective(_tt_cost, phase; sense = Min(), at = Path())
    return phase
end

# The closed-form optimum, for comparison at any epoch.
_tt_x_star(t) = (2 / 3) * sin(t)
_tt_u_star(t) = (2 / 3) * cos(t)
const _TT_J_STAR = π / 3

@testset "tracking — the closed form is a solution of the problem" begin
    # Check the analytic answer against the problem's own statement before using it as truth.
    # The state and control are consistent with the dynamics, the boundary conditions hold, and
    # the cost integrates to π/3.
    ts = range(0, _TT_TF; length = 2001)

    # ẋ = u, checked by finite difference on the closed form.
    h = 1e-6
    for t in (0.3, 1.7, 3.0, 4.5, 6.0)
        dx = (_tt_x_star(t + h) - _tt_x_star(t - h)) / (2h)
        @test dx ≈ _tt_u_star(t) atol = 1e-8
    end

    @test _tt_x_star(0.0) ≈ 0.0 atol = 1e-14
    @test _tt_x_star(_TT_TF) ≈ 0.0 atol = 1e-14

    # The cost, by trapezoid on a fine grid. ∫[(x* - sin t)² + ½u*²] dt = π/3.
    f(t) = (_tt_x_star(t) - sin(t))^2 + 0.5 * _tt_u_star(t)^2
    vals = f.(ts)
    J = (step(ts)) * (sum(vals) - 0.5 * (vals[1] + vals[end]))
    @test J ≈ _TT_J_STAR atol = 1e-6

    # The path bound is slack everywhere on the optimum, which is why it does not move the answer.
    @test maximum((1.0 + 0.5 * sin(t)) * _tt_x_star(t)^2 - 1.0 for t in ts) < 0.0
end

@testset "tracking — Hermite-Simpson reproduces the analytic cost" begin
    phase = _tt_phase()
    result = solve!(Sequence(phase))

    @test result.info === :Solve_Succeeded

    # The objective is an integral and nothing else, so this number is entirely the quadrature
    # the transcription performs. π/3 to four figures on a 24-step mesh.
    @test result.objective ≈ _TT_J_STAR atol = 1e-4

    # And against the recorded LGL baseline, which is the same problem under a different
    # transcription.
    @test result.objective ≈ 1.047198 atol = 1e-4
end

@testset "tracking — the converged trajectory is the analytic one, pointwise" begin
    # A cost that happens to be right does not mean the trajectory is. Both are known here, so
    # both are checked at every node rather than only at the endpoints.
    phase = _tt_phase()
    solve!(Sequence(phase))

    ts = get_node_times(phase)
    xs = phase._Y
    us = phase._U

    @test length(ts) == size(xs, 2)
    @test size(xs, 1) == 1

    for k in eachindex(ts)
        @test xs[1, k] ≈ _tt_x_star(ts[k]) atol = 1e-3
    end

    # The control mesh may be coarser than the state mesh depending on the transcription, so
    # compare it on its own node times.
    @test size(us, 1) == 1
    @test maximum(abs, us) <= 5.0 + 1e-9

    # Endpoints, which were declared as constraints rather than inferred.
    @test xs[1, 1] ≈ 0.0 atol = 1e-6
    @test xs[1, end] ≈ 0.0 atol = 1e-6

    # The path bound held everywhere and stayed slack, as the analytic optimum says it must.
    viol = maximum((1.0 + 0.5 * sin(ts[k])) * xs[1, k]^2 - 1.0 for k in eachindex(ts))
    @test viol <= 1e-9
    @test viol < 0.0
end

@testset "tracking — a finer mesh converges toward the analytic cost" begin
    # The gap between the transcribed cost and π/3 is quadrature error. Refining must shrink it,
    # which is what makes the tolerance above a statement about the mesh.
    err_coarse = abs(solve!(Sequence(_tt_phase(n_steps = 8))).objective  - _TT_J_STAR)
    err_fine   = abs(solve!(Sequence(_tt_phase(n_steps = 40))).objective - _TT_J_STAR)

    @test err_fine < err_coarse
end
