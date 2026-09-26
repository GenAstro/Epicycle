# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Two derivative paths that nothing else reaches.
#
# `test_correctness_ad_fallback.jl` covers a phase with no declared partials, but its problem has
# no path constraint and a fixed time span, so two of the AD routines never run:
# `_path_jacobian_ad_chunk`, which differentiates a path constraint node by node when no analytic
# Jacobian is registered, and `_add_nonauto_time_correction!`, which corrects the time columns of
# the defect Jacobian when the dynamics depend on `t` explicitly.
#
# Truth: **finite differences of the constraint vector**. Declared-versus-AD cannot be used here,
# because for the time correction there is no declared version to compare against — the
# correction is what the analytic path itself relies on. Central differences of `get_functions`
# with respect to the decision vector are independent of both.
#
# The failure these guard against is silent. A wrong derivative does not raise; the solver takes
# small steps, or converges somewhere plausible and wrong, which is how the dense-sparsity
# defect behaved for months before anyone measured it.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: get_decision_vector, get_functions, get_jacobian, initialize!,
    set_decision_vector!
using Test

# ── The problem ──────────────────────────────────────────────────────────────
#
# Dynamics that depend on `t` explicitly, so ∂F/∂t is non-zero and the time correction has
# something to correct. The autonomous twin differs only in that term.

struct _TdState{T}   <: AbstractState;   x::T end
struct _TdControl{T} <: AbstractControl; u::T end

_td_dyn_nonauto!(dy, y::_TdState, u::_TdControl, p, t, model) = (dy[1] = u.u + cos(t))
_td_dyn_auto!(dy, y::_TdState, u::_TdControl, p, t, model)    = (dy[1] = u.u)

_td_pos(c)   = [state(c).x]
_td_speed(c) = [control(c).u]                       # no @partial: differentiated
_td_cost(c)  = state(c).x^2 + 0.5 * control(c).u^2

"""
Build a small phase, optionally with a varying final time and a path constraint.

The mesh is deliberately coarse. These tests compare derivative entries rather than converged
answers, so accuracy of the transcription is not what is under test and a fine mesh only makes
the finite-difference loop slower.
"""
function _td_phase(; nonauto::Bool = true, vary_tf::Bool = false, with_path::Bool = false,
                     n_steps::Int = 4)
    phase = CollocationPhase(
        name          = Symbol(:td_, nonauto, :_, vary_tf, :_, with_path),
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = nonauto ? _td_dyn_nonauto! : _td_dyn_auto!,
        state         = _TdState, control = _TdControl,
        tspan         = (0.0, 1.0))

    Vary(state,   phase; guess = [0.2 0.4],
                         lower_bound = [-5.0], upper_bound = [5.0])
    Vary(control, phase; guess = reshape([0.3 0.1], 1, 2),
                         lower_bound = [-5.0], upper_bound = [5.0])
    vary_tf && Vary(final_time, phase; guess = 1.0, lower_bound = 0.5, upper_bound = 2.0)
    with_path && Constraint(_td_speed, phase; upper_bound = 3.0, at = Path())

    Constraint(_td_pos, phase; equals = [0.0], at = Initial())
    Objective(_td_cost, phase; sense = Min(), at = Path())
    return phase
end

"""
Central-difference the constraint vector with respect to the decision vector.

The sequence is left at `x` afterwards, so the caller can ask for the analytic Jacobian at the
same point.
"""
function _fd_jacobian(seq, x; h = 1e-6)
    n = length(x)
    set_decision_vector!(seq, copy(x))
    m = length(get_functions(seq))
    J = zeros(m, n)
    for j in 1:n
        xp = copy(x); xp[j] += h
        set_decision_vector!(seq, xp)
        fp = copy(get_functions(seq))

        xm = copy(x); xm[j] -= h
        set_decision_vector!(seq, xm)
        fm = copy(get_functions(seq))

        J[:, j] = (fp .- fm) ./ (2h)
    end
    set_decision_vector!(seq, copy(x))
    return J
end

"Build, initialise, and return the sequence together with its decision vector."
function _td_ready(phase)
    seq = Sequence(phase)
    initialize!(seq)
    return seq, copy(get_decision_vector(seq))
end

@testset "AD path Jacobian — a path constraint with no declared partial" begin
    # `_path_jacobian_ad_chunk`. The constraint is on the control, so its rows pick out the
    # control columns and nothing else; getting that block wrong is a constraint the solver
    # cannot see itself violating.
    seq, x = _td_ready(_td_phase(with_path = true))

    J  = Matrix(get_jacobian(seq))
    Jd = _fd_jacobian(seq, x)

    @test size(J) == size(Jd)
    @test J ≈ Jd atol = 1e-5

    # The path rows are not all zero, so the block really was assembled rather than skipped.
    @test any(!iszero, J)
end

@testset "AD path Jacobian — non-autonomous dynamics with a varying final time" begin
    # `_add_nonauto_time_correction!`. Without the correction the time column of the defect
    # Jacobian is short by ∂F/∂t · τ, which for these dynamics is order one.
    seq, x = _td_ready(_td_phase(nonauto = true, vary_tf = true))

    J  = Matrix(get_jacobian(seq))
    Jd = _fd_jacobian(seq, x)

    @test size(J) == size(Jd)
    @test J ≈ Jd atol = 1e-5
end

@testset "AD path Jacobian — the time correction vanishes for autonomous dynamics" begin
    # The comment on `_add_nonauto_time_correction!` says it is safe to apply always because
    # ∂F/∂t is zero when the dynamics do not mention `t`. That is a claim about the code, so it
    # is worth running rather than reading: the same problem with `t` removed must still match
    # finite differences.
    seq, x = _td_ready(_td_phase(nonauto = false, vary_tf = true))

    J  = Matrix(get_jacobian(seq))
    Jd = _fd_jacobian(seq, x)

    @test J ≈ Jd atol = 1e-5
end

@testset "AD path Jacobian — a path constraint and a varying final time together" begin
    # Both paths in one problem, which is the combination a real low-thrust case presents.
    seq, x = _td_ready(_td_phase(nonauto = true, vary_tf = true, with_path = true))

    J  = Matrix(get_jacobian(seq))
    Jd = _fd_jacobian(seq, x)

    @test size(J) == size(Jd)
    @test J ≈ Jd atol = 1e-5
end
