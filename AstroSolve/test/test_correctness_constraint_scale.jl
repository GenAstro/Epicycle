# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Constraint scaling on a collocation phase.
#
# `Constraint(...; scale = s)` divides the constraint by `s` before the solver sees it, and its
# bounds with it, the convention the shooting phases already use. On a collocation phase the
# keyword was accepted and dropped, so a script that scaled its constraints got none of it and was
# not told.
#
# Truth: **two formulations of one problem**. Scaled and unscaled, the rows the solver sees differ
# by exactly the scale, their bounds and Jacobian rows likewise, and the answer is the same.

using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: get_constraint_bounds, get_decision_vector, get_functions, get_jacobian,
    initialize!, set_decision_vector!
using Test

struct _CsState{T}   <: AbstractState;   x::T end
struct _CsControl{T} <: AbstractControl; u::T end

_cs_dyn!(dy, y::_CsState, u::_CsControl, p, t, model) = (dy[1] = u.u)
_cs_start(c) = [state(c).x]
@partial(_cs_start, state) do c
    [1.0]
end
_cs_rate(c) = [control(c).u]
@partial(_cs_rate, control) do c
    [1.0]
end

"""The Hull problem, with its start and its control bound scaled by `s_start` and `s_rate`."""
function _cs_phase(; s_start = nothing, s_rate = nothing)
    phase = CollocationPhase(name = :scaled, transcription = HermiteSimpson(n_steps = 8),
                             dynamics = _cs_dyn!, state = _CsState, control = _CsControl,
                             tspan = (0.0, 1.0))
    Vary(state,   phase; guess = [1.5 0.25], lower_bound = [-10.0], upper_bound = [10.0])
    Vary(control, phase; guess = [-1.25 -1.25], lower_bound = [-10.0], upper_bound = [10.0])
    Constraint(_cs_start, phase; equals = [1.5], at = Initial(), scale = s_start)
    Constraint(_cs_rate, phase; lower_bound = -2.0, upper_bound = 2.0, at = Path(), scale = s_rate)
    Objective(c -> 2.5 * state(c).x^2, phase; sense = Min())
    Objective(c -> 0.5 * control(c).u^2, phase; sense = Min(), at = Path())
    return phase
end

function _cs_rows(phase)
    seq = Sequence(phase); initialize!(seq)
    set_decision_vector!(seq, get_decision_vector(seq))
    F = copy(get_functions(seq)); lb, ub = get_constraint_bounds(seq)
    return F, copy(lb), copy(ub), Matrix(get_jacobian(seq))
end

@testset "Constraint scale — the rows, bounds and Jacobian are divided by the scale" begin
    F0, lb0, ub0, J0 = _cs_rows(_cs_phase())
    F1, lb1, ub1, J1 = _cs_rows(_cs_phase(s_start = 10.0, s_rate = 4.0))
    @test length(F0) == length(F1)

    # Every row either matches or is divided by one of the two scales, and both scales occur.
    ratio = [F0[i] == 0 ? (F1[i] == 0 ? 1.0 : NaN) : F1[i] / F0[i] for i in eachindex(F0)]
    scaled_rows = findall(r -> !(r ≈ 1.0), ratio)
    @test !isempty(scaled_rows)
    @test all(r -> r ≈ 1.0 || r ≈ 0.1 || r ≈ 0.25, ratio)
    @test any(r -> r ≈ 0.1, ratio) && any(r -> r ≈ 0.25, ratio)

    # The bounds and the Jacobian rows scale with their rows.
    for i in scaled_rows
        s = 1 / ratio[i]
        @test lb1[i] ≈ lb0[i] / s
        @test ub1[i] ≈ ub0[i] / s
        @test J1[i, :] ≈ J0[i, :] ./ s
    end
end

@testset "Constraint scale — the scaled problem has the unscaled answer" begin
    a = _cs_phase()
    b = _cs_phase(s_start = 10.0, s_rate = 4.0)
    ra = solve!(Sequence(a); method = Optimize(print_level = 0))
    rb = solve!(Sequence(b); method = Optimize(print_level = 0))
    @test ra.info === :Solve_Succeeded
    @test rb.info === :Solve_Succeeded
    @test rb.objective ≈ ra.objective rtol = 1e-8
    # To IPOPT's tolerance, which a scale changes: dividing a row by 10 loosens its convergence test
    # tenfold in physical units. The two land 2e-6 apart; the optimum is 0.25.
    @test get_final_state(b).x ≈ get_final_state(a).x atol = 1e-5
    @test get_final_state(b).x ≈ 0.25 atol = 1e-5
end

@testset "Constraint scale — a wrong length or a non-positive entry is refused" begin
    phase = _cs_phase()
    @test_throws ArgumentError Constraint(_cs_start, phase; equals = [1.5], at = Initial(),
                                          scale = [1.0, 2.0])
    @test_throws ArgumentError Constraint(_cs_start, phase; equals = [1.5], at = Initial(),
                                          scale = 0.0)
    @test_throws ArgumentError Constraint(_cs_rate, phase; lower_bound = -2.0, upper_bound = 2.0,
                                          at = Path(), scale = -1.0)
end
