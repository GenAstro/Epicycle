# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Several shooting phases assembled into one NLP.
#
# `ShootingManager` is what turns a list of phases plus the constraints that tie them together
# into a single decision vector, constraint vector and Jacobian. Nothing about that is physics:
# it is bookkeeping, and every claim here is an exact arithmetic identity between what the
# manager reports and what the phases it was built from report separately.
#
# Truth: **analytic**, in the sense that every figure follows from the phase counts, plus an
# **AD-versus-analytic** check on the cross-phase constraint's Jacobian, which is the one part
# that is a derivative rather than a concatenation.
#
# The bookkeeping is where multi-phase problems go wrong, and the failure is quiet. A column range
# placed one entry off, or a variable shared between two phases counted twice, produces a Jacobian
# that is the wrong shape in a way the solver absorbs: it takes a worse step, not an error. So the
# assertions are about offsets, lengths and overlaps rather than about an answer.
#
# This file depends on the fixture in test_correctness_sims_flanagan.jl, which runtests.jl
# includes first.

using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: build_sparsity_pattern, get_constraint_bounds, get_decision_vector,
    get_functions, get_jacobian, get_jacobian_values!, get_objective, get_objective_gradient,
    get_variable_bounds, n_constraints, nlp_length, set_decision_vector!, variable_list
using Test

const _SM = AstroSolve

"""Two independent Sims-Flanagan phases and the manager that assembles them."""
function _sm_setup(; seq_constraints = _SM.SequenceConstraint[])
    p1, nh1, _ = _sf_phase(n_segments = 20, throttle = 0.5)
    p2, nh2, _ = _sf_phase(n_segments = 12, throttle = 0.3)
    get_functions(p1); get_functions(p2)          # warm the propagated arcs
    sm = _SM.ShootingManager(Any[p1, p2], seq_constraints)
    return sm, p1, p2, nh1, nh2
end

@testset "shooting manager — the assembled problem is the sum of its phases" begin
    sm, p1, p2, nh1, nh2 = _sm_setup()

    # Variables concatenate. The two phases declare disjoint variables, so nothing is deduplicated
    # and the global count is the plain sum.
    @test sm.n_vars == nlp_length(p1) + nlp_length(p2)
    @test length(sm.ordered_vars) == length(variable_list(p1)) + length(variable_list(p2))

    # Constraints concatenate, in phase order, with no sequence constraints to add.
    @test sm.n_funs == n_constraints(p1) + n_constraints(p2)
    @test sm.phase_fun_offsets == [1, n_constraints(p1) + 1]

    # Every variable gets a column range, the ranges tile the decision vector exactly, and none of
    # them overlaps. This is the bookkeeping a misplaced offset breaks.
    ranges = [sm.global_var_ranges[objectid(v)] for v in sm.ordered_vars]
    @test sum(length, ranges) == sm.n_vars
    @test sort(reduce(vcat, collect.(ranges))) == collect(1:sm.n_vars)
end

@testset "shooting manager — the decision vector round-trips through the phases" begin
    sm, p1, p2, _, _ = _sm_setup()

    x = get_decision_vector(sm)
    @test length(x) == sm.n_vars

    # Reading the manager's vector and writing it back leaves every phase where it was, which is
    # what lets a solver own the vector while the phases own the physics.
    f_before = get_functions(sm)
    set_decision_vector!(sm, copy(x))
    @test get_decision_vector(sm) ≈ x atol = 0.0
    @test get_functions(sm) ≈ f_before atol = 0.0

    # And the manager's constraint vector really is the two phases' vectors end to end.
    F  = get_functions(sm)
    F1 = get_functions(p1)
    F2 = get_functions(p2)
    @test length(F) == length(F1) + length(F2)
    @test F[1:length(F1)] ≈ F1 atol = 0.0
    @test F[length(F1)+1:end] ≈ F2 atol = 0.0

    # Changing one phase's slice must move only that phase's rows. A column range that reached
    # into the wrong phase would show up here and nowhere else.
    x2 = copy(x)
    x2[1] += 0.05
    set_decision_vector!(sm, x2)
    F_moved = get_functions(sm)
    @test F_moved[length(F1)+1:end] ≈ F2 atol = 1e-12
    @test !(F_moved[1:length(F1)] ≈ F1)
end

@testset "shooting manager — bounds come from the phases that own the variables" begin
    sm, p1, p2, _, _ = _sm_setup()

    lx, ux = get_variable_bounds(sm)
    @test length(lx) == length(ux) == sm.n_vars
    @test all(lx .<= ux)

    # Each phase's own bounds appear in its own columns.
    l1, u1 = get_variable_bounds(p1)
    l2, u2 = get_variable_bounds(p2)
    @test lx[1:nlp_length(p1)] ≈ l1 atol = 0.0
    @test ux[1:nlp_length(p1)] ≈ u1 atol = 0.0
    @test lx[nlp_length(p1)+1:end] ≈ l2 atol = 0.0
    @test ux[nlp_length(p1)+1:end] ≈ u2 atol = 0.0

    clo, chi = get_constraint_bounds(sm)
    @test length(clo) == length(chi) == sm.n_funs
    @test all(clo .<= chi)
end

@testset "shooting manager — a cross-phase constraint adds rows and couples columns" begin
    # Mass continuity across the join: the second phase must start with the mass the first one
    # ends with. That is the simplest real sequence constraint and it touches variables in both
    # phases, which is what makes it a test of the coupling rather than of a phase.
    mass_continuity(c1, c2) = [c2.m0 - c1.mf]

    p1, _, _ = _sf_phase(n_segments = 20, throttle = 0.5)
    p2, _, _ = _sf_phase(n_segments = 12, throttle = 0.3)
    get_functions(p1); get_functions(p2)

    sc = _SM.SequenceConstraint(p1, p2, mass_continuity, [0.0], [0.0], "mass", Dict{UInt64,Function}())
    sm = _SM.ShootingManager(Any[p1, p2], [sc])

    # One more row than the phases alone, placed after both of them.
    @test sm.n_funs == n_constraints(p1) + n_constraints(p2) + 1
    @test sm.seq_con_offsets == [sm.n_funs]

    # The row evaluates to what the function says it should.
    F = get_functions(sm)
    @test length(F) == sm.n_funs
    @test F[end] ≈ p2._m0 - p1._mf atol = 1e-9

    # Its bounds are the equality that was declared.
    clo, chi = get_constraint_bounds(sm)
    @test clo[end] == 0.0 && chi[end] == 0.0

    # The variable count did not change: a sequence constraint adds rows, not columns.
    @test sm.n_vars == nlp_length(p1) + nlp_length(p2)
end

@testset "shooting manager — the assembled Jacobian is the derivative it claims to be" begin
    # The manager places each phase's Jacobian block and each sequence constraint's row into one
    # matrix. A block placed at the wrong offset is not an error: the solver takes a worse step
    # and still converges somewhere, so only differencing the assembled matrix catches it.
    mass_continuity(c1, c2) = [c2.m0 - c1.mf]
    p1, _, _ = _sf_phase(n_segments = 12, throttle = 0.5)
    p2, _, _ = _sf_phase(n_segments = 12, throttle = 0.3)
    get_functions(p1); get_functions(p2)
    sc = _SM.SequenceConstraint(p1, p2, mass_continuity, [0.0], [0.0], "mass",
                                Dict{UInt64,Function}())
    sm = _SM.ShootingManager(Any[p1, p2], [sc])

    x = get_decision_vector(sm)
    J = get_jacobian(sm)
    @test size(J) == (sm.n_funs, sm.n_vars)
    @test all(isfinite, J)

    # Central differences of the same constraint vector, over every column.
    J_fd = zeros(sm.n_funs, sm.n_vars)
    h = 1e-6
    for j in 1:sm.n_vars
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(sm, xp); Fp = copy(get_functions(sm))
        set_decision_vector!(sm, xm); Fm = copy(get_functions(sm))
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(sm, x)

    scale = max(maximum(abs, J_fd), 1.0)
    @test maximum(abs, J .- J_fd) / scale < 1e-5

    # The sequence-constraint row in particular, which is the part no single-phase test can reach.
    @test maximum(abs, J[end, :] .- J_fd[end, :]) < 1e-5

    # Mass continuity reaches phase 1 only, because phase 2's initial mass is pinned rather than
    # varied — a sequence constraint differentiates only what is actually a variable.
    @test any(abs.(J[end, 1:nlp_length(p1)]) .> 0)
    @test all(J[end, nlp_length(p1)+1:end] .== 0.0)

    # A constraint on both final masses does reach both, which is the two-sided coupling.
    both_masses(c1, c2) = [c2.mf - c1.mf]
    q1, _, _ = _sf_phase(n_segments = 12, throttle = 0.5)
    q2, _, _ = _sf_phase(n_segments = 12, throttle = 0.3)
    get_functions(q1); get_functions(q2)
    sc2 = _SM.SequenceConstraint(q1, q2, both_masses, [0.0], [0.0], "dm",
                                 Dict{UInt64,Function}())
    sm2 = _SM.ShootingManager(Any[q1, q2], [sc2])
    J2 = get_jacobian(sm2)
    @test any(abs.(J2[end, 1:nlp_length(q1)]) .> 0)
    @test any(abs.(J2[end, nlp_length(q1)+1:end]) .> 0)
end

@testset "shooting manager — the sparsity pattern claims every nonzero" begin
    # A sparsity pattern is one-sided in the direction that matters: an entry the pattern calls
    # zero and the derivative calls nonzero never reaches the solver at all, and the solver finds
    # its way anyway, so the answer is the same and nothing fails. That is the failure this
    # checks for, and it is why the comparison is not symmetric.
    mass_continuity(c1, c2) = [c2.m0 - c1.mf]
    p1, _, _ = _sf_phase(n_segments = 12, throttle = 0.5)
    p2, _, _ = _sf_phase(n_segments = 12, throttle = 0.3)
    get_functions(p1); get_functions(p2)
    sc = _SM.SequenceConstraint(p1, p2, mass_continuity, [0.0], [0.0], "mass",
                                Dict{UInt64,Function}())
    sm = _SM.ShootingManager(Any[p1, p2], [sc])

    rows, cols, vals = build_sparsity_pattern(sm)
    @test length(rows) == length(cols) == length(vals)
    @test all(1 .<= rows .<= sm.n_funs)
    @test all(1 .<= cols .<= sm.n_vars)
    @test all(vals .!= 0.0)

    # Refilling the same pattern reproduces the values it was built from.
    refilled = zeros(Float64, length(rows))
    get_jacobian_values!(refilled, sm, rows, cols)
    @test refilled ≈ vals atol = 0.0

    # Every nonzero the dense Jacobian has is claimed by the pattern.
    J = get_jacobian(sm)
    claimed = Set(zip(rows, cols))
    unclaimed = count((i, j) for j in 1:sm.n_vars, i in 1:sm.n_funs
                      if abs(J[i, j]) > 0.0 && (i, j) ∉ claimed)
    @test unclaimed == 0
end

@testset "shooting manager — the objective is the phases' objective" begin
    sm, p1, p2, _, _ = _sm_setup()

    obj = get_objective(sm)
    @test isfinite(obj)

    g = get_objective_gradient(sm)
    @test length(g) == sm.n_vars
    @test all(isfinite, g)

    # The gradient is a derivative, so it must agree with a difference of the objective.
    x = get_decision_vector(sm)
    h = 1e-6
    for j in (1, sm.n_vars ÷ 2, sm.n_vars)
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(sm, xp); op = get_objective(sm)
        set_decision_vector!(sm, xm); om = get_objective(sm)
        set_decision_vector!(sm, x)
        @test g[j] ≈ (op - om) / (2h) atol = 1e-4
    end
end

@testset "shooting manager — it refuses to assemble nothing" begin
    # The exception type and its message.
    @test_throws ArgumentError _SM.ShootingManager(Any[])
    msg = try; _SM.ShootingManager(Any[]); catch e; sprint(showerror, e); end
    @test occursin("at least one phase", msg)

    sm, _, _, _, _ = _sm_setup()
    @test_throws ArgumentError set_decision_vector!(sm, zeros(sm.n_vars + 1))
    msg2 = try; set_decision_vector!(sm, zeros(sm.n_vars + 1))
           catch e; sprint(showerror, e); end
    @test occursin("one entry per NLP variable", msg2)
    @test occursin(string(sm.n_vars), msg2)
end
