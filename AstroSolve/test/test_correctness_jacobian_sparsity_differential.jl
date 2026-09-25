# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The declared sparsity pattern reaches the same answer as declaring the Jacobian dense.
#
# Truth: **self-consistency, end to end**. `test_correctness_jacobian_sparsity.jl` checks that the
# declared pattern holds every position the Jacobian produces. That is necessary and not sufficient:
# it says nothing about whether the values are delivered to the positions they belong to. A pattern
# whose rows and columns are each individually plausible but paired wrongly passes containment and
# gives the solver a transposed or shifted Jacobian.
#
# So this file solves each problem twice, once with the declared pattern and once with everything
# declared, and requires the two answers to agree. Identical arithmetic reaches the solver either
# way; only the amount of it IPOPT is told about differs, so a disagreement is a mistake in the
# pattern and not a tolerance.
#
# The problems are the library's own, reached through the phase builders the other correctness files
# already define, so the shapes are the ones that ship rather than ones invented here. Between them
# they cover a varied final time and a fixed one, path constraints and none, a Mayer cost, varied
# static parameters, bang-bang control, and a phase carrying no declared partials at all, which
# reaches the same assembly through automatic differentiation.
#
# Not covered here, and worth knowing: several collocation phases joined by linkages. The suite's
# multi-phase cases go through the shooting or optimal-control managers rather than
# `_solve_phases!`, so they exercise a path this pattern does not touch, and `test_correctness_mixed_links.jl`
# additionally shares its `_ML_` constant prefix with the moon landing file, so the two cannot be
# loaded together. Linkage rows are declared dense against both phases they join, which is the
# conservative choice, but it is untested.
#
# Including those files runs their testsets, which is deliberate: it re-checks each suite under the
# declared pattern in the same process.

using LinearAlgebra
using EpicycleBase
using AstroSolve
using Printf
using Test

const _SPD_DIR = @__DIR__

# Each builder comes from the correctness file that owns that problem, and each is taken only if the
# suite has not already loaded it. Run on its own, this file pulls in what it needs; run from
# `runtests.jl`, which loads all six earlier, it re-runs none of their testsets.
for (builder, file) in ((:_br_phase,  "test_correctness_brachistochrone.jl"),    # varied tf, path + ends
                        (:_oa_phase,  "test_correctness_obstacle_avoidance.jl"), # fixed tf, path constraints
                        (:_hu_phase,  "test_correctness_bolza.jl"),              # Mayer plus Lagrange
                        (:_ml_phase,  "test_correctness_moon_landing.jl"),       # bang-bang control
                        (:_pid_phase, "test_correctness_parameter_id.jl"),       # varied static parameters
                        (:_ad_phase,  "test_correctness_ad_fallback.jl"))        # no declared partials
    isdefined(@__MODULE__, builder) || include(joinpath(_SPD_DIR, file))
end

"""Solve `build()` twice, with the declared pattern and with everything declared, and compare.

`build` returns whatever `Sequence` takes, so a single phase or a tuple of linked ones both work.
"""
function _differential(name, build; max_iter = 500)
    results = Dict{Bool, Any}()
    omitted = -1
    for dense in (false, true)
        dense_jacobian!(dense)
        seq = Sequence(build()...)
        dense || (omitted = check_sparsity(seq; verbose = false))
        r = solve!(seq; method = Optimize(print_level = 0, max_iter = max_iter))
        results[dense] = (obj = r.objective, x = copy(r.variables), info = r.info)
    end
    dense_jacobian!(false)

    sparse_r, dense_r = results[false], results[true]
    d_obj = abs(sparse_r.obj - dense_r.obj)
    d_x   = maximum(abs, sparse_r.x .- dense_r.x)
    both  = sparse_r.info === dense_r.info

    @printf("  %-26s  %-18s omitted %-3d  dobj %.2e  max dx %.2e\n",
            name, string(sparse_r.info), omitted, d_obj, d_x)

    # The pattern must never omit a position, whatever the solve then does with it.
    @test omitted == 0

    # Agreement is only meaningful when both runs converged. An unconverged run stops at whatever
    # iterate it reached, and two different factorisations reach different ones legitimately, so
    # asserting on those would be testing the solver's path rather than the pattern.
    if sparse_r.info === :Solve_Succeeded && dense_r.info === :Solve_Succeeded
        @test d_obj <= 1e-6 * max(1.0, abs(dense_r.obj))
        @test d_x   <= 1e-4 * max(1.0, maximum(abs, dense_r.x))
    else
        @test both      # at least the two paths agree about what happened
    end
    return nothing
end

@testset "constraint Jacobian sparsity — the declared pattern gives the dense answer" begin
    println()
    _differential("brachistochrone",       () -> (_br_phase(n_steps = 20),))
    _differential("obstacle avoidance",    () -> (_oa_phase(n_steps = 35),))
    _differential("Hull problem",          () -> (_hu_phase(n_steps = 20),))
    _differential("moon landing",          () -> (_ml_phase(n_steps = 30),))
    _differential("parameter id",          () -> (_pid_phase(n_steps = 20),))
    _differential("AD fallback",           () -> (_ad_phase(:bare; n_steps = 24),))
    println()
end
