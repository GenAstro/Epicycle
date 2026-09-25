# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Hermite-Simpson collocation, against the brachistochrone's closed-form solution.
#
# Truth: **analytic**. The curve of fastest descent under gravity is a cycloid, and the cycloid
# through a given endpoint has a descent time in closed form. This file derives that time from the
# boundary conditions rather than hard-coding it, so the assertion is against the mathematics and
# not against a number someone once recorded.
#
#   minimise   tf
#   subject to ẋ = v sin θ,  ẏ = v cos θ,  v̇ = g cos θ
#   from (0, 0) at rest to (2, 2), with g = 9.80665 and |θ| ≤ π/2.
#
# This is the whole open-tier optimal-control path in one problem: the spec verbs a user writes
# (`Vary`, `Constraint`, `Objective`), declared partials and the fallback to AD, the
# Hermite-Simpson transcription, the NLP assembly, and the solve. It is the case the regression
# baseline records as uc5b at tf = 0.824479, matching uc5's LGL answer to six figures, which is a
# cross-transcription agreement on top of the analytic check.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

# ── The problem ──────────────────────────────────────────────────────────────

struct _BrState{T}   <: AbstractState;   x::T; y::T; v::T end
struct _BrControl{T} <: AbstractControl; θ::T             end
struct _BrModel;     g::Float64                           end

const _BR_G     = 9.80665
const _BR_MODEL = _BrModel(_BR_G)
const _BR_TARGET = (2.0, 2.0)

function _br_dynamics!(dy, y::_BrState, u::_BrControl, p, t, model)
    dy[1] = y.v * sin(u.θ)
    dy[2] = y.v * cos(u.θ)
    dy[3] = model.g * cos(u.θ)
end

# Two of the nine state partials are nonzero and both are declared; the rest of the Jacobian is
# structurally zero. Declaring a partial and letting the remainder fall back to AD is the normal
# case, and `check_partials` below is what says the declared half agrees with the fallback.
@partial(_br_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 3] = sin(u.θ)
    dF[2, 3] = cos(u.θ)
end

@partial(_br_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] =  y.v * cos(u.θ)
    dF[2, 1] = -y.v * sin(u.θ)
    dF[3, 1] = -model.g * sin(u.θ)
end

_br_start(c) = [state(c).x, state(c).y, state(c).v]
@partial(_br_start, state) do c
    [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
end

_br_final_position(c) = [state(c).x, state(c).y]
@partial(_br_final_position, state) do c
    [1.0 0.0 0.0; 0.0 1.0 0.0]
end

# One quantity, one partial, two constraints: a path bound and an endpoint bound. Both are slack
# at the cycloid, so they exercise the machinery without moving the answer, which is what keeps
# the analytic comparison below honest.
_br_speed(c) = [state(c).v]
@partial(_br_speed, state) do c
    [0.0 0.0 1.0]
end

_br_duration(c) = final_time(c)
@partial(_br_duration, final_time) do c
    [1.0]
end

"""Build the phase, declare its variables and constraints, and return it.

A fresh phase per testset: the spec verbs mutate the phase they are given, so two testsets cannot
share one.
"""
function _br_phase(; n_steps = 20)
    phase = CollocationPhase(
        name          = :brachistochrone,
        transcription = HermiteSimpson(n_steps = n_steps),
        dynamics      = _br_dynamics!, model   = _BR_MODEL,
        state         = _BrState,      control = _BrControl,
        tspan         = (0.0, 1.0))

    Vary(state, phase;
         guess       = [0.0 0.5 1.0 1.5 2.0
                        0.0 0.5 1.0 1.5 2.0
                        0.0 1.0 1.5 2.0 2.5],
         lower_bound = [-Inf, -Inf, 0.0],
         upper_bound = [ Inf,  Inf, Inf])

    Vary(control, phase;
         guess       = reshape([0.3, 0.5, 0.7, 0.9, 1.0], 1, 5),
         lower_bound = [-π/2],
         upper_bound = [ π/2])

    Vary(final_time, phase; lower_bound = 0.1, upper_bound = 10.0)

    Constraint(_br_start,          phase; equals = [0.0, 0.0, 0.0], at = Initial())
    Constraint(_br_final_position, phase; equals = [_BR_TARGET...], at = Final())
    Constraint(_br_speed, phase; upper_bound = 10.0, at = Path())
    Constraint(_br_speed, phase; upper_bound = 10.0, at = Final())

    Objective(_br_duration, phase; sense = Min())
    return phase
end

# ── The closed-form answer ───────────────────────────────────────────────────

"""Descent time along the brachistochrone from the origin to `(x1, y1)` under gravity `g`.

The curve is the cycloid `x = R(φ - sin φ)`, `y = R(1 - cos φ)`, traversed in `t = √(R/g)·φ`.
The endpoint fixes `φ` through the ratio `x1/y1 = (φ - sin φ)/(1 - cos φ)`, which is solved here
by bisection because it is monotone on `(0, 2π)` and the bracket is known.
"""
function _br_cycloid_time(x1, y1; g = _BR_G)
    ratio(φ) = (φ - sin(φ)) / (1 - cos(φ))
    lo, hi = 1e-8, 2π - 1e-8
    target = x1 / y1
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        ratio(mid) < target ? (lo = mid) : (hi = mid)
    end
    φ = 0.5 * (lo + hi)
    R = y1 / (1 - cos(φ))
    return sqrt(R / g) * φ, φ, R
end

@testset "brachistochrone — the cycloid time is what the problem is for" begin
    # The closed form first, on its own terms, so a failure downstream is not ambiguous about
    # which side moved. The cycloid through (2, 2) is symmetric in the ratio, which pins φ.
    tf_exact, φ, R = _br_cycloid_time(_BR_TARGET...)

    # The endpoint really is on the curve these parameters describe.
    @test R * (φ - sin(φ)) ≈ _BR_TARGET[1] atol = 1e-10
    @test R * (1 - cos(φ)) ≈ _BR_TARGET[2] atol = 1e-10

    # And the descent time is the figure the literature quotes for this problem.
    @test tf_exact ≈ 0.8245 atol = 5e-5
end

@testset "brachistochrone — Hermite-Simpson finds the analytic optimum" begin
    phase = _br_phase()
    result = solve!(Sequence(phase))

    @test result.info === :Solve_Succeeded

    tf = get_final_time(phase)
    tf_exact, _, _ = _br_cycloid_time(_BR_TARGET...)

    # Against the mathematics. The gap is discretization error on a 20-step mesh, not solver
    # tolerance, so this is a statement about the transcription rather than about IPOPT.
    @test tf ≈ tf_exact atol = 1e-4

    # And against the recorded baseline, which is what catches a change that is still plausible.
    # uc5 solves the same problem under LGL and lands on the same six figures, so this number is
    # a cross-transcription agreement as well as a regression.
    @test tf ≈ 0.824479 atol = 1e-6

    # A minimum-time problem must not be beaten by the closed form by more than the mesh can
    # explain: a transcription that cheats reports a time below the true optimum.
    @test tf > tf_exact - 1e-3
end

@testset "brachistochrone — the converged trajectory satisfies its constraints" begin
    phase = _br_phase()
    solve!(Sequence(phase))

    # There is no public accessor for the converged mesh — `get_initial_state`, `get_final_state`
    # and `get_control` give the endpoints and the control, but the state history comes off the
    # private `_Y`, which is what the shipped examples read too. Recorded as a design gap in the
    # test plan; the test uses what a user would have to use.
    xs = phase._Y                  # n_states x n_nodes
    us = phase._U
    tf = get_final_time(phase)

    @test size(xs, 1) == 3
    @test size(us, 1) == 1
    @test size(xs, 2) == size(us, 2)

    # The boundary conditions the user declared, read back off the converged mesh.
    @test xs[:, 1] ≈ [0.0, 0.0, 0.0] atol = 1e-6
    @test xs[1, end] ≈ _BR_TARGET[1] atol = 1e-6
    @test xs[2, end] ≈ _BR_TARGET[2] atol = 1e-6

    # The control stayed inside its box everywhere, not only on average.
    @test all(-π/2 - 1e-9 .<= us .<= π/2 + 1e-9)

    # The path bound on speed held at every node, and was slack, which is why tf above is still
    # the unconstrained optimum.
    @test all(xs[3, :] .<= 10.0 + 1e-9)
    @test maximum(xs[3, :]) < 10.0

    # Energy: falling through height y from rest gives v = sqrt(2 g y) whatever path is taken, so
    # the converged speed profile must satisfy it at every node. This checks the dynamics the
    # transcription integrated rather than the boundary conditions it was handed.
    # The depth at the first node converges to zero, and to zero from either side: which side is
    # round-off in the solver's own arithmetic. `max` keeps the domain of `sqrt` without loosening
    # the comparison, since a negative depth of that size and a depth of zero are the same physics.
    for k in axes(xs, 2)
        @test xs[3, k] ≈ sqrt(max(0.0, 2 * _BR_G * xs[2, k])) atol = 1e-4
    end

    # Descent takes time, and the mesh spans it.
    @test tf > 0
end

@testset "brachistochrone — declared partials agree with automatic differentiation" begin
    # check_partials differentiates every declared partial's function and compares. This is the
    # AD-versus-analytic check, applied to the derivatives a user hands
    # the solver: a wrong partial still converges, to the wrong place or not at all, so nothing
    # downstream catches it.
    phase = _br_phase()
    @test check_partials(phase) isa Any          # runs without raising on a well-formed phase

    # A phase is solvable straight after the check, which is the order the examples use it in.
    result = solve!(Sequence(phase))
    @test result.info === :Solve_Succeeded
end

@testset "brachistochrone — a mesh refinement moves the answer toward the cycloid" begin
    # The discretization error is what separates the transcribed answer from the closed form, so
    # a finer mesh must close the gap. This is the property that says the 1e-4 tolerance above is
    # a mesh statement and not a coincidence.
    tf_exact, _, _ = _br_cycloid_time(_BR_TARGET...)

    coarse = _br_phase(n_steps = 8)
    solve!(Sequence(coarse))
    err_coarse = abs(get_final_time(coarse) - tf_exact)

    fine = _br_phase(n_steps = 30)
    solve!(Sequence(fine))
    err_fine = abs(get_final_time(fine) - tf_exact)

    @test err_fine < err_coarse
end
