# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What declaring the constraint Jacobian sparsity is worth, per problem.
#
# The two ratio columns are dense over sparse, so a larger number is better: 8x means the declared
# pattern runs eight times faster, or allocates an eighth as much, as declaring everything.
#
# Each problem is solved with the declared pattern and with everything declared, and the pair is run
# twice so the first run pays for compilation and the second is the one reported. The answers are the
# same either way — `test_correctness_jacobian_sparsity_differential.jl` is what establishes that —
# so the only thing changing here is how much linear algebra IPOPT does per iteration.
#
# Run with the shared environment active:
#
#   using TestEnv; TestEnv.activate("AstroSolve")
#   include("AstroSolve/test/Benchmark_JacobianSparsity.jl")

using LinearAlgebra
using EpicycleBase
using AstroSolve
using Printf

const _BJS_DIR = @__DIR__

include(joinpath(_BJS_DIR, "test_correctness_brachistochrone.jl"))
include(joinpath(_BJS_DIR, "test_correctness_obstacle_avoidance.jl"))
include(joinpath(_BJS_DIR, "test_correctness_bolza.jl"))
include(joinpath(_BJS_DIR, "test_correctness_moon_landing.jl"))
include(joinpath(_BJS_DIR, "test_correctness_parameter_id.jl"))

"""Time one solve of `build()` under the current declaration, and report what was declared."""
function _timed_solve(build; max_iter = 500)
    seq = Sequence(build()...)
    AstroSolve.evaluate!(seq, AstroSolve.get_decision_vector(seq))
    ng, nx = size(AstroSolve.get_jacobian(seq))
    nnz = length(AstroSolve.jacobian_pattern(seq, ng, nx).rows)
    t = @timed solve!(seq; method = Optimize(print_level = 0, max_iter = max_iter))
    return (seconds = t.time, bytes = t.bytes, nnz = nnz, ng = ng, nx = nx,
            obj = t.value.objective, info = t.value.info)
end

function _compare(name, build; max_iter = 500)
    best = Dict{Bool, Any}()
    # Twice, keeping the faster. The first pass compiles; the second is the measurement.
    for _ in 1:2, dense in (false, true)
        dense_jacobian!(dense)
        run = _timed_solve(build; max_iter = max_iter)
        prev = get(best, dense, nothing)
        (prev === nothing || run.seconds < prev.seconds) && (best[dense] = run)
    end
    dense_jacobian!(false)

    sp, dn = best[false], best[true]
    @printf("  %-20s %4d x %-5d %6.2f%%  %6.2f / %6.2f s %5.1fx  %8.1f / %9.1f MB %6.1fx
",
            name, sp.ng, sp.nx, 100 * sp.nnz / dn.nnz,
            sp.seconds, dn.seconds, dn.seconds / max(sp.seconds, 1e-9),
            sp.bytes / 2^20, dn.bytes / 2^20, dn.bytes / max(sp.bytes, 1))
    return nothing
end

println()
@printf("  %-20s %-13s %-7s %-16s %-6s %-20s
",
        "problem", "Jacobian", "declared", "time  sparse / dense", "", "allocated  sparse / dense")
println("  ", "-"^112)
_compare("brachistochrone 20",  () -> (_br_phase(n_steps = 20),))
_compare("brachistochrone 60",  () -> (_br_phase(n_steps = 60),))
_compare("obstacle avoidance",  () -> (_oa_phase(n_steps = 35),))
_compare("Hull problem",        () -> (_hu_phase(n_steps = 20),))
_compare("moon landing",        () -> (_ml_phase(n_steps = 30),))
_compare("parameter id",        () -> (_pid_phase(n_steps = 20),))
println()
