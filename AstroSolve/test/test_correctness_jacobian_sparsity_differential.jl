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
# pattern and not a tolerance. That holds only where the optimum is unique. Obstacle avoidance has
# a family of equally optimal paths, and is checked at fixed points instead (the last testset).
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
using Random
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
    _differential("Hull problem",          () -> (_hu_phase(n_steps = 20),))
    _differential("moon landing",          () -> (_ml_phase(n_steps = 30),))
    _differential("parameter id",          () -> (_pid_phase(n_steps = 20),))
    _differential("AD fallback",           () -> (_ad_phase(:bare; n_steps = 24),))
    println()
end

# Obstacle avoidance is not in the solve-twice comparison above, because its optimum is not
# unique. The cost is V² tf whatever the heading does (test_correctness_obstacle_avoidance.jl),
# so every feasible path is optimal, and two solves that differ only in how much of the Jacobian
# IPOPT is told about can legitimately stop on different paths. They did: the two solutions
# differed by 4.4e-2 with identical objectives, on Windows under Julia 1.13, having agreed on
# Linux under 1.12. Where a solve stops is a property of the platform's floating point, so the
# comparison was testing that rather than the pattern.
#
# So for this problem the check is made where the answer is unique: at a fixed point, before any
# solve. What the solver receives through the declared pattern has to be the dense Jacobian:
#
#   - the pattern is in the order the value vector is filled in, column-major, because SNOW
#     pairs value k with pattern position k;
#   - each value `get_jacobian_values!` writes equals the dense entry at its position. The two are
#     separate fills of the same chunks, which is why comparing them means something;
#   - the dense Jacobian is zero everywhere the pattern leaves out, so nothing the solver needs
#     is withheld.
#
# At several points inside the variable bounds, for the reason test_correctness_jacobian_sparsity.jl
# gives: an entry can be zero at one point and not at another.
@testset "constraint Jacobian sparsity — obstacle avoidance delivers the dense Jacobian" begin
    seq    = Sequence(_oa_phase(n_steps = 35))
    @test check_sparsity(seq; verbose = false) == 0

    x0     = AstroSolve.get_decision_vector(seq)
    lx, ux = AstroSolve.get_variable_bounds(seq)
    ng     = length(first(AstroSolve.get_constraint_bounds(seq)))
    nx     = length(x0)

    pattern = AstroSolve.jacobian_pattern(seq, ng, nx; dense = false)
    index   = AstroSolve.JacobianIndex(pattern, ng, nx)
    @test issorted(collect(zip(pattern.cols, pattern.rows)))

    # The column of each slot, from the compressed columns.
    slot_col = similar(index.rowval)
    for c in 1:nx, k in index.colptr[c]:(index.colptr[c + 1] - 1)
        slot_col[k] = c
    end
    declared = falses(ng, nx)
    for k in eachindex(index.rowval)
        declared[index.rowval[k], slot_col[k]] = true
    end

    rng  = Random.MersenneTwister(20260929)
    vals = zeros(length(index.rowval))
    for s in 0:8
        x = s == 0 ? copy(x0) :
            clamp.(x0 .+ max.(abs.(x0), 1.0) .* (2 .* rand(rng, nx) .- 1), lx, ux)
        AstroSolve.evaluate!(seq, x)
        J = AstroSolve.get_jacobian(seq)
        AstroSolve.get_jacobian_values!(vals, seq, index)

        @test all(k -> vals[k] == J[index.rowval[k], slot_col[k]], eachindex(vals))
        @test all(iszero, J[.!declared])
    end
end
