# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A Bolza cost: an endpoint term and an integral, declared as two objectives on one phase.
#
# Truth: **analytic**, and derived here rather than quoted. The Hull problem is
#
#   minimise   J = 2.5 x_f² + ½ ∫₀¹ u² dt
#   subject to ẋ = u,  x(0) = 1.5
#
# With ẋ = u and no bound active, the optimal control is constant: any variation that keeps x_f
# fixed can only raise ∫u², and the integral of a fixed area is smallest when spread evenly. So
# writing u ≡ c gives x_f = 1.5 + c and J(c) = 2.5(1.5 + c)² + ½c², a scalar quadratic whose
# minimum is c = -1.25, x_f = 0.25, J = 0.9375. The test solves that quadratic rather than
# hard-coding its root.
#
# What this covers that the other collocation cases do not: **two objectives on one phase**, a
# Mayer term and a Lagrange term, which the assembly has to add rather than choose between. A
# transcription that silently dropped one would still converge, to 0.78125 or to 0.15625, so the
# total is the number that catches it.
#
# The problem is the regression suite's uc14, recorded at 0.9375 under LGL. This file transcribes
# it with Hermite-Simpson.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

struct _HuState{T}   <: AbstractState;   x::T end
struct _HuControl{T} <: AbstractControl; u::T end

const _HU_X0 = 1.5

_hu_dynamics!(dy, y::_HuState, u::_HuControl, p, t, model) = (dy[1] = u.u)

@partial(_hu_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = 1.0
end
# ∂ẋ/∂x is structurally zero and is left undeclared, so the state Jacobian falls back to AD.
# That mix is the ordinary case and is what check_partials reports on.

_hu_initial(c) = [state(c).x]
@partial(_hu_initial, state) do c
    [1.0]
end

_hu_terminal_cost(c) = 2.5 * state(c).x^2
@partial(_hu_terminal_cost, state) do c
    [5.0 * state(c).x]
end

_hu_running_cost(c) = 0.5 * control(c).u^2
@partial(_hu_running_cost, control) do c
    [control(c).u]
end

"""Build the Hull phase on a Hermite-Simpson mesh of `n_steps` intervals."""
function _hu_phase(; n_steps = 20)
    phase = CollocationPhase(
        name          = :hull,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _hu_dynamics!,
        state         = _HuState,      control = _HuControl,
        tspan         = (0.0, 1.0))

    Vary(state,   phase; guess = [1.5 0.25],
                         lower_bound = [-10.0], upper_bound = [10.0])
    Vary(control, phase; guess = reshape([-1.25 -1.25], 1, 2),
                         lower_bound = [-10.0], upper_bound = [10.0])

    Constraint(_hu_initial, phase; equals = [_HU_X0], at = Initial())

    Objective(_hu_terminal_cost, phase; sense = Min())
    Objective(_hu_running_cost,  phase; sense = Min(), at = Path())
    return phase
end

"""The optimal constant control, final state and cost, from the scalar quadratic."""
function _hu_analytic(; x0 = _HU_X0)
    # J(c) = 2.5 (x0 + c)^2 + 0.5 c^2 ;  dJ/dc = 5(x0 + c) + c = 0
    c   = -5x0 / 6
    xf  = x0 + c
    J   = 2.5 * xf^2 + 0.5 * c^2
    return c, xf, J
end

@testset "Bolza — the analytic optimum solves the problem it claims to" begin
    c, xf, J = _hu_analytic()

    # The stationarity condition the derivation used, checked rather than assumed.
    @test 5 * (_HU_X0 + c) + c ≈ 0.0 atol = 1e-14

    # And it really is a minimum: perturbing the constant control raises the cost either way.
    Jof(cc) = 2.5 * (_HU_X0 + cc)^2 + 0.5 * cc^2
    @test Jof(c + 1e-3) > J
    @test Jof(c - 1e-3) > J

    # The figures this problem is quoted with.
    @test c  ≈ -1.25 atol = 1e-12
    @test xf ≈ 0.25  atol = 1e-12
    @test J  ≈ 0.9375 atol = 1e-12
end

@testset "Bolza — the transcription adds both objective terms" begin
    phase = _hu_phase()
    result = solve!(Sequence(phase))
    _, xf_exact, J_exact = _hu_analytic()

    @test result.info === :Solve_Succeeded

    # The total. A transcription that dropped the Mayer term would report 0.78125 and one that
    # dropped the Lagrange term 0.15625, so this single number separates all three cases.
    @test result.objective ≈ J_exact atol = 1e-6
    @test result.objective ≈ 0.9375 atol = 1e-6

    # The state the two terms trade off against each other to reach.
    @test get_final_state(phase).x ≈ xf_exact atol = 1e-6

    # Both halves, read off the converged answer, so a compensating pair of errors that happens
    # to sum correctly is caught as well.
    xs = phase._Y
    us = phase._U
    @test 2.5 * xs[1, end]^2 ≈ 2.5 * xf_exact^2 atol = 1e-6
    @test all(abs.(us .- (-1.25)) .< 1e-5)
end

@testset "Bolza — the converged trajectory is the straight line it should be" begin
    # A constant control on ẋ = u makes x affine in t, so the whole trajectory is known, not just
    # its endpoints. This is what says the mesh integrated the dynamics rather than only
    # satisfying the boundary condition and the cost.
    phase = _hu_phase()
    solve!(Sequence(phase))

    ts = get_node_times(phase)
    xs = phase._Y
    c, _, _ = _hu_analytic()

    @test length(ts) == size(xs, 2)
    for k in eachindex(ts)
        @test xs[1, k] ≈ _HU_X0 + c * ts[k] atol = 1e-5
    end

    # The declared initial condition, which is the one thing here that was constrained.
    @test xs[1, 1] ≈ _HU_X0 atol = 1e-8
    @test get_initial_state(phase).x ≈ _HU_X0 atol = 1e-8
end

@testset "Bolza — the answer does not depend on the mesh" begin
    # The optimum is a straight line and a constant control, both of which Hermite-Simpson
    # represents exactly, so refining must not move the answer at all. That is a stronger claim
    # than the convergence checks in the other collocation files and is available only because
    # this solution lies in the transcription's own function space.
    _, _, J_exact = _hu_analytic()
    for n in (6, 20, 40)
        r = solve!(Sequence(_hu_phase(n_steps = n)))
        @test r.objective ≈ J_exact atol = 1e-6
    end
end
