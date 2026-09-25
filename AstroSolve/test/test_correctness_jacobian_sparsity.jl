# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The declared constraint Jacobian sparsity pattern contains every entry the Jacobian can hold.
#
# Truth: **self-consistency**. There is no external reference for a sparsity pattern. What makes one
# correct is that it never omits an entry the dense Jacobian produces. An omitted entry does not
# raise: the solver never asks for that derivative, nothing reports it, and the run converges to the
# wrong answer or fails for reasons that look like a badly posed problem.
#
# Evaluating the dense Jacobian at one point does not establish what the pattern must hold. A
# structurally nonzero entry can be numerically zero at any given point — the partial of
# `u.u - sigma^2` with respect to the thrust is exactly zero wherever the thrust is, which is where a
# minimum-fuel guess starts. So this file unions the nonzero positions over many points inside the
# variable bounds and requires the pattern to contain the union. Reading the Jacobian at one point is
# exactly the mistake `build_sparsity_pattern` makes, which is why that function is unused.
#
# No solve happens here. Containment is a property of the assembly, so the check needs only
# `evaluate!` and `get_jacobian`, which keeps it fast enough to run over several problem shapes.

using Random
using LinearAlgebra
using EpicycleBase
using AstroSolve
using Test

# ── A problem carrying every row type the assembly produces ──────────────────
#
# A double integrator: position and speed driven by acceleration. Small, linear, and exactly
# differentiable, so nothing here is about whether the partials are right. What it is for is shape —
# Hermite-Simpson defect rows, a boundary constraint at each end, a path constraint, and a varied
# final time so the time columns are live.

struct _SpState{T}   <: AbstractState;   x::T; v::T end
struct _SpControl{T} <: AbstractControl; a::T       end

function _sp_dynamics!(dy, y::_SpState, u::_SpControl, p, t, model)
    dy[1] = y.v
    dy[2] = u.a
end

@partial(_sp_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 2] = 1.0
end

@partial(_sp_dynamics!, control) do dF, y, u, p, t, model
    dF[2, 1] = 1.0
end

_sp_both(c)  = [state(c).x, state(c).v]
@partial(_sp_both, state) do c
    Matrix{Float64}(I, 2, 2)
end

_sp_speed(c) = [state(c).v]
@partial(_sp_speed, state) do c
    [0.0 1.0]
end

# An objective is scalar valued, where a constraint is a vector, so the fixed-time shape needs its
# own cost rather than reusing the path bound `_sp_speed`.
_sp_cost(c) = state(c).v^2
@partial(_sp_cost, state) do c
    [0.0  2 * state(c).v]
end

_sp_duration(c) = final_time(c)
@partial(_sp_duration, final_time) do c
    [1.0]
end

"""Build the phase fresh. The spec verbs mutate the phase they are given, so it cannot be shared."""
function _sp_phase(; n_steps = 5, with_path = true, vary_time = true)
    phase = CollocationPhase(name = :double_integrator,
                             transcription = HermiteSimpson(n_steps = n_steps),
                             dynamics = _sp_dynamics!,
                             state = _SpState, control = _SpControl,
                             tspan = (0.0, 1.0))

    Vary(state, phase;
         guess = [0.0 0.25 0.5 0.75 1.0
                  0.0 1.00 1.0 1.00 0.0],
         lower_bound = [-10.0, -10.0], upper_bound = [10.0, 10.0])

    Vary(control, phase;
         guess = reshape([1.0, 0.5, 0.0, -0.5, -1.0], 1, 5),
         lower_bound = [-5.0], upper_bound = [5.0])

    # A fixed tspan and a varied final time are different variable layouts, and the pattern has to
    # index both. Leaving the final time alone is the common case in a library example.
    vary_time && Vary(final_time, phase; lower_bound = 0.5, upper_bound = 5.0)

    Constraint(_sp_both, phase; equals = [0.0, 0.0], at = Initial())
    Constraint(_sp_both, phase; equals = [1.0, 0.0], at = Final())
    with_path && Constraint(_sp_speed, phase; upper_bound = 3.0, at = Path())

    # Minimising the duration needs the duration to be free. With the time fixed there is nothing to
    # trade, so the cost becomes the integrated acceleration instead.
    vary_time ? Objective(_sp_duration, phase; sense = Min()) :
                Objective(_sp_cost, phase; sense = Min(), at = Path())
    return phase
end

# ── Measuring what the Jacobian actually holds ────────────────────────────────

"""Positions nonzero in the dense Jacobian at any of several points inside the variable bounds.

Sampling several points is the point of this helper. One evaluation reports the entries that happen
to be nonzero there, which is a subset of the structure and not the structure.
"""
function _nonzero_union(seq; n_samples = 24, seed = 20260922)
    rng    = MersenneTwister(seed)
    x0     = AstroSolve.get_decision_vector(seq)
    lx, ux = AstroSolve.get_variable_bounds(seq)
    found  = Set{Tuple{Int, Int}}()

    for k in 0:n_samples
        # Stay inside the declared bounds. An out-of-bounds point is not wrong for a structural
        # question, but it can put a function where it was never meant to be evaluated, and a NaN
        # there reads as a structural zero.
        x = if k == 0
            copy(x0)
        else
            clamp.(x0 .+ max.(abs.(x0), 1.0) .* (2 .* rand(rng, length(x0)) .- 1), lx, ux)
        end
        AstroSolve.evaluate!(seq, x)
        J = AstroSolve.get_jacobian(seq)
        for c in axes(J, 2), r in axes(J, 1)
            J[r, c] != 0.0 && push!(found, (r, c))
        end
    end
    AstroSolve.evaluate!(seq, x0)          # leave the sequence where it started
    J = AstroSolve.get_jacobian(seq)
    return found, size(J)
end

"""Check one problem shape, and report what the declaration cost."""
function _check_containment(name, seq)
    actual, (ng, nx) = _nonzero_union(seq)
    pattern  = AstroSolve.jacobian_pattern(seq, ng, nx)
    declared = Set(zip(pattern.rows, pattern.cols))
    omitted  = setdiff(actual, declared)

    @info name ng nx dense = ng * nx actual = length(actual) declared = length(declared) omitted = length(omitted) fraction = round(length(declared) / (ng * nx), digits = 4)

    # Every declared position has to be a position. An index outside the problem is how a mistake in
    # the variable-block layout shows up, and it is otherwise invisible: the solver is handed a
    # pattern it cannot use and the only trace is a nonzero count that does not add up.
    @test all(r -> 1 <= r <= ng, pattern.rows)
    @test all(c -> 1 <= c <= nx, pattern.cols)

    # Nothing the Jacobian can produce may be missing. This is the assertion the solver's
    # correctness rests on; the density assertion below is about cost only.
    @test isempty(omitted)

    # And the declaration has to be smaller than dense, or there was no reason to make one.
    @test length(declared) < ng * nx
    return nothing
end

@testset "constraint Jacobian sparsity — the pattern contains every entry the Jacobian holds" begin
    _check_containment("double integrator, 5 steps",        Sequence(_sp_phase(n_steps = 5)))
    _check_containment("double integrator, 20 steps",       Sequence(_sp_phase(n_steps = 20)))
    _check_containment("double integrator, no path bound",  Sequence(_sp_phase(n_steps = 8, with_path = false)))
    _check_containment("double integrator, fixed final time", Sequence(_sp_phase(n_steps = 10, vary_time = false)))
end
