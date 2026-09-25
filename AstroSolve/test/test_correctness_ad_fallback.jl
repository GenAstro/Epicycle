# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What happens when a user declares no partials at all.
#
# Every derivative a collocation phase needs can be declared with `@partial` or left to automatic
# differentiation. The other collocation files in this suite declare nearly everything, which is
# the careful case; this one declares nothing, which is the case a user reaches first and the one
# that decides whether the library is usable before anyone has written a Jacobian by hand.
#
# Truth: **AD-versus-analytic**, the §2 tier, applied at three levels — the converged answer, the
# constraint Jacobian, and the objective gradient. The same problem is transcribed twice from
# identical data, once with hand-written partials and once with none, and the two must agree.
#
# The problem is the tracking case from test_correctness_trajectory_tracking.jl, whose optimal
# state, control and cost are all in closed form, so "agree" can be checked against the
# mathematics rather than only against each other.
#
# The last two testsets cover a Bolza objective with no declared partials. Its Mayer half once
# contributed nothing to the gradient, and declaring the running cost before the terminal cost
# once dropped the running cost; both solved "successfully" to the wrong answer.
#
# Covers autodiff.jl, which every other file in this suite bypasses by declaring its partials.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: get_decision_vector, get_functions, get_jacobian, get_objective,
    get_objective_gradient, initialize!, set_decision_vector!
using Test

# ── The same problem, declared twice ─────────────────────────────────────────
#
# `@partial` registers against the function object, so the two versions need distinct functions.
# Everything else about them is identical.

const _AD_TF = 2π

struct _AdState{T}   <: AbstractState;   x::T end
struct _AdControl{T} <: AbstractControl; u::T end

_ad_dyn_bare!(dy, y::_AdState, u::_AdControl, p, t, model) = (dy[1] = u.u)
_ad_dyn_full!(dy, y::_AdState, u::_AdControl, p, t, model) = (dy[1] = u.u)
@partial(_ad_dyn_full!, control) do dF, y, u, p, t, model
    dF[1, 1] = 1.0
end

_ad_pos_bare(c) = [state(c).x]
_ad_pos_full(c) = [state(c).x]
@partial(_ad_pos_full, state) do c
    [1.0]
end

_ad_cost_bare(c) = (state(c).x - sin(c.t))^2 + 0.5 * control(c).u^2
_ad_cost_full(c) = (state(c).x - sin(c.t))^2 + 0.5 * control(c).u^2
@partial(_ad_cost_full, state) do c
    [2.0 * (state(c).x - sin(c.t))]
end
@partial(_ad_cost_full, control) do c
    [control(c).u]
end

_ad_x_star(t) = (2 / 3) * sin(t)
const _AD_J_STAR = π / 3

"""Build the tracking phase from the `bare` (no partials) or `full` (declared) function set."""
function _ad_phase(which::Symbol; n_steps = 24)
    dyn, pos, cost = which === :bare ?
        (_ad_dyn_bare!, _ad_pos_bare, _ad_cost_bare) :
        (_ad_dyn_full!, _ad_pos_full, _ad_cost_full)

    phase = CollocationPhase(
        name          = Symbol(:tracking_, which),
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = dyn,
        state         = _AdState,      control = _AdControl,
        tspan         = (0.0, _AD_TF))

    Vary(state,   phase; guess = [0.0 0.0],
                         lower_bound = [-5.0], upper_bound = [5.0])
    Vary(control, phase; guess = reshape([0.0 0.0], 1, 2),
                         lower_bound = [-5.0], upper_bound = [5.0])

    Constraint(pos, phase; equals = [0.0], at = Initial())
    Constraint(pos, phase; equals = [0.0], at = Final())
    Objective(cost, phase; sense = Min(), at = Path())
    return phase
end

@testset "AD fallback — a phase with no declared partials solves the problem" begin
    # The case a user reaches before writing any derivative by hand. If this does not work, the
    # library is unusable until someone has done the calculus.
    phase = _ad_phase(:bare)
    result = solve!(Sequence(phase))

    @test result.info === :Solve_Succeeded
    @test result.objective ≈ _AD_J_STAR atol = 1e-4

    # And the trajectory is the analytic one, not merely a cost that happens to land.
    ts = get_node_times(phase)
    xs = phase._Y
    for k in eachindex(ts)
        @test xs[1, k] ≈ _ad_x_star(ts[k]) atol = 1e-3
    end
end

@testset "AD fallback — declared and differentiated give the same answer" begin
    # The comparison this file exists for. Two transcriptions of one problem, differing only in
    # where their derivatives come from, must land in the same place.
    bare = _ad_phase(:bare)
    full = _ad_phase(:full)

    r_bare = solve!(Sequence(bare))
    r_full = solve!(Sequence(full))

    @test r_bare.info === :Solve_Succeeded
    @test r_full.info === :Solve_Succeeded
    @test r_bare.objective ≈ r_full.objective atol = 1e-6

    # Not only the cost: the whole converged mesh.
    @test size(bare._Y) == size(full._Y)
    @test bare._Y ≈ full._Y atol = 1e-4
    @test bare._U ≈ full._U atol = 1e-4
end

@testset "AD fallback — the differentiated constraint Jacobian is the declared one" begin
    # One level below the answer. Both phases are put at the same point in the decision space and
    # their constraint Jacobians compared entry by entry: same constraints, same state, different
    # provenance for every number.
    bare = _ad_phase(:bare); sb = Sequence(bare); initialize!(sb)
    full = _ad_phase(:full); sf = Sequence(full); initialize!(sf)

    x = get_decision_vector(sf)
    set_decision_vector!(sb, copy(x))
    set_decision_vector!(sf, copy(x))

    @test get_functions(sb) ≈ get_functions(sf) atol = 1e-12   # same constraints

    Jb = get_jacobian(sb)
    Jf = get_jacobian(sf)
    @test size(Jb) == size(Jf)
    @test maximum(abs, Jb .- Jf) < 1e-8

    # And the shared value is a derivative rather than a shared mistake.
    J_fd = zeros(size(Jf))
    h = 1e-6
    for j in eachindex(x)
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(sf, xp); Fp = copy(get_functions(sf))
        set_decision_vector!(sf, xm); Fm = copy(get_functions(sf))
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(sf, x); get_functions(sf)

    @test maximum(abs, Jf .- J_fd) < 1e-5
    @test maximum(abs, Jb .- J_fd) < 1e-5
end

@testset "AD fallback — a Lagrange objective's gradient is differentiated" begin
    # The running cost is the part a user writes last and declares least often, so its fallback is
    # the one that matters most. Here it is differentiated correctly: the two gradients agree with
    # each other and with a difference of the objective.
    bare = _ad_phase(:bare); sb = Sequence(bare); initialize!(sb)
    full = _ad_phase(:full); sf = Sequence(full); initialize!(sf)

    x = get_decision_vector(sf)
    set_decision_vector!(sb, copy(x)); get_functions(sb)
    set_decision_vector!(sf, copy(x)); get_functions(sf)

    gb = get_objective_gradient(sb)
    gf = get_objective_gradient(sf)
    @test length(gb) == length(gf) == length(x)
    @test maximum(abs, gb .- gf) < 1e-8

    g_fd = similar(x)
    h = 1e-6
    for j in eachindex(x)
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(sb, xp); get_functions(sb); op = get_objective(sb)
        set_decision_vector!(sb, xm); get_functions(sb); om = get_objective(sb)
        g_fd[j] = (op - om) / (2h)
    end
    set_decision_vector!(sb, x); get_functions(sb)

    @test maximum(abs, gb .- g_fd) < 1e-6
end

# ── A Bolza objective, both halves differentiated ────────────────────────────

struct _AdHuState{T}   <: AbstractState;   x::T end
struct _AdHuControl{T} <: AbstractControl; u::T end

_adhu_dyn!(dy, y::_AdHuState, u::_AdHuControl, p, t, model) = (dy[1] = u.u)
_adhu_initial(c) = [state(c).x]
_adhu_terminal(c) = 2.5 * state(c).x^2
_adhu_running(c)  = 0.5 * control(c).u^2

"""
The Hull Bolza problem with no partial declared on either half of its objective. `order` says
which half is declared first.
"""
function _adhu_phase(; n_steps = 10, order = :mayer_first)
    phase = CollocationPhase(
        name          = :hull_bare,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _adhu_dyn!,
        state         = _AdHuState,     control = _AdHuControl,
        tspan         = (0.0, 1.0))
    Vary(state,   phase; guess = [1.5 0.25],
                         lower_bound = [-10.0], upper_bound = [10.0])
    Vary(control, phase; guess = reshape([-1.0 -1.0], 1, 2),
                         lower_bound = [-10.0], upper_bound = [10.0])
    Constraint(_adhu_initial, phase; equals = [1.5], at = Initial())
    if order === :mayer_first
        Objective(_adhu_terminal, phase; sense = Min())
        Objective(_adhu_running,  phase; sense = Min(), at = Path())
    else
        Objective(_adhu_running,  phase; sense = Min(), at = Path())
        Objective(_adhu_terminal, phase; sense = Min())
    end
    return phase
end

@testset "AD fallback — a Bolza objective's Mayer half is differentiated" begin
    # The Mayer half used to contribute nothing to the gradient when it had no declared partial:
    # the Bolza gradient differentiated the running cost alone. The value was right, so the solve
    # reported success with the terminal cost unoptimised, x(1) = 1.5 where the optimum is 0.25.
    # Each half now takes its declared partial or is differentiated, and the two are summed.
    phase = _adhu_phase()
    seq = Sequence(phase); initialize!(seq)
    x = get_decision_vector(seq); get_functions(seq)

    # The value includes both halves. The guess interpolates the state from 1.5 down to 0.25, so
    # the terminal term is 2.5(0.25)² and the running term ½∫u² over a unit interval at u = -1.
    @test get_objective(seq) ≈ 2.5 * phase._Y[1, end]^2 + 0.5 atol = 1e-6

    g = get_objective_gradient(seq)
    g_fd = similar(x)
    h = 1e-6
    for j in eachindex(x)
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(seq, xp); get_functions(seq); op = get_objective(seq)
        set_decision_vector!(seq, xm); get_functions(seq); om = get_objective(seq)
        g_fd[j] = (op - om) / (2h)
    end
    set_decision_vector!(seq, x); get_functions(seq)
    @test maximum(abs, g .- g_fd) < 1e-6

    solved = _adhu_phase()
    r = solve!(Sequence(solved))
    @test r.info === :Solve_Succeeded
    @test get_final_state(solved).x ≈ 0.25 atol = 1e-5
    @test r.objective ≈ 0.9375 atol = 1e-6
end

@testset "AD fallback — the order the two halves are declared in does not matter" begin
    # Declaring the running cost first, then the terminal cost, used to replace the Bolza
    # objective with a Mayer one. The running cost was dropped and the solve reported success at
    # J = 0 with x(1) = 0. The Objective docstring says order does not matter; this holds it to it.
    for order in (:mayer_first, :lagrange_first)
        phase = _adhu_phase(order = order)
        r = solve!(Sequence(phase))
        @test r.info === :Solve_Succeeded
        @test get_final_state(phase).x ≈ 0.25 atol = 1e-5
        @test r.objective ≈ 0.9375 atol = 1e-6
    end
end
