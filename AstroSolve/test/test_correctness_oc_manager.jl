# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# One NLP from phases of different kinds.
#
# `OCManager` is what lets a trajectory be transcribed one way on one leg and another way on the
# next — a collocated arc handed to a Sims-Flanagan arc, joined by continuity. That is the
# capability `ShootingManager` does not have, and the piece that makes it work is
# `_unified_boundary_context`: a collocation phase carries mass inside its state and a shooting
# phase carries it beside, so each maps its own endpoint cache into one shape and a link between
# them is then a plain equality over every component rather than a special case per pair.
#
# Truth: **analytic**, in that every count follows from the phases, plus a differencing of the
# assembled Jacobian. Nothing here is physics. The failure mode it guards against is quiet: a
# column range placed one entry off produces a Jacobian the solver absorbs as a worse step rather
# than rejecting as an error.
#
# This file depends on the fixtures in test_correctness_brachistochrone.jl and
# test_correctness_sims_flanagan.jl, which runtests.jl includes first.

using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: add_continuity!, build_sparsity_pattern, get_constraint_bounds,
    get_decision_vector, get_functions, get_jacobian, get_jacobian_values!, get_objective,
    get_objective_gradient, get_variable_bounds, n_constraints, nlp_length, set_decision_vector!
using Test

const _OC = AstroSolve

"""A collocation phase and a Sims-Flanagan phase, each ready to be read at its endpoints."""
function _oc_phases()
    coll = _br_phase(n_steps = 6)
    solve!(Sequence(coll))                 # populates the node cache the context reads
    shoot, _, _ = _sf_phase(n_segments = 12, throttle = 0.4)
    get_functions(shoot)
    return coll, shoot
end

@testset "OC manager — a unified context reports both phase kinds in one shape" begin
    coll, shoot = _oc_phases()

    c_coll  = _OC._unified_boundary_context(coll)
    c_shoot = _OC._unified_boundary_context(shoot)

    # A collocation phase reports its state as it stores it.
    @test length(c_coll.y0) == 3
    @test c_coll.y0 ≈ coll._Y[:, 1] atol = 0.0
    @test c_coll.yf ≈ coll._Y[:, end] atol = 0.0
    @test c_coll.t0 ≈ coll._t0 atol = 0.0
    @test c_coll.tf ≈ coll._tf atol = 0.0

    # A Sims-Flanagan phase keeps mass outside the state, so the unified context appends it and
    # the two kinds become comparable component by component. This is the whole reason the
    # unified context exists.
    @test length(c_shoot.y0) == 7
    @test c_shoot.y0[7] ≈ shoot._m0 atol = 1e-12
    @test c_shoot.yf[7] ≈ shoot._mf atol = 1e-12
    @test c_shoot.t0 ≈ shoot._t0 atol = 0.0
    @test c_shoot.tf ≈ shoot._tf atol = 0.0

    # Both report epochs that run forward.
    @test c_coll.tf > c_coll.t0
    @test c_shoot.tf > c_shoot.t0
end

@testset "OC manager — a mixed sequence assembles into one problem" begin
    coll, shoot = _oc_phases()
    om = _OC.OCManager(Any[coll, shoot])

    @test om.n_vars == nlp_length(coll) + nlp_length(shoot)
    @test om.n_funs == n_constraints(coll) + n_constraints(shoot)
    @test om.phase_fun_offsets == [1, n_constraints(coll) + 1]

    # Column ranges tile the decision vector exactly and do not overlap.
    ranges = [om.global_var_ranges[objectid(v)] for v in om.ordered_vars]
    @test sum(length, ranges) == om.n_vars
    @test sort(reduce(vcat, collect.(ranges))) == collect(1:om.n_vars)

    # The decision vector round-trips, and the constraint vector is the two phases end to end.
    x = get_decision_vector(om)
    @test length(x) == om.n_vars
    F = get_functions(om)
    @test length(F) == om.n_funs
    set_decision_vector!(om, copy(x))
    @test get_decision_vector(om) ≈ x atol = 0.0
    @test get_functions(om) ≈ F atol = 0.0

    lx, ux = get_variable_bounds(om)
    @test length(lx) == length(ux) == om.n_vars
    @test all(lx .<= ux)
    clo, chi = get_constraint_bounds(om)
    @test length(clo) == length(chi) == om.n_funs
    @test all(clo .<= chi)
end

@testset "OC manager — a continuity link adds rows and reads both contexts" begin
    coll, shoot = _oc_phases()

    # The two phases carry different state dimensions, which the generic continuity helper would
    # refuse. Writing the link directly is what a mixed problem does when the legs do not have the
    # same state: here the epochs must meet and the collocated arc's last speed must match the
    # shooting arc's first, which is a two-row condition rather than a full state equality.
    joined(c1, c2) = [c1.tf - c2.t0, c1.yf[3] - c2.y0[4]]

    lk = _OC.ContinuityLink(coll, shoot, joined, [0.0, 0.0], [0.0, 0.0], "join")
    om = _OC.OCManager(Any[coll, shoot], [lk])

    @test om.n_funs == n_constraints(coll) + n_constraints(shoot) + 2
    @test om.link_fun_offsets == [om.n_funs - 1]

    F = get_functions(om)
    c1 = _OC._unified_boundary_context(coll)
    c2 = _OC._unified_boundary_context(shoot)
    @test F[end-1] ≈ c1.tf - c2.t0 atol = 1e-9
    @test F[end]   ≈ c1.yf[3] - c2.y0[4] atol = 1e-9

    # Its bounds are the equalities declared, and a link adds rows rather than columns.
    clo, chi = get_constraint_bounds(om)
    @test clo[end-1:end] == [0.0, 0.0]
    @test chi[end-1:end] == [0.0, 0.0]
    @test om.n_vars == nlp_length(coll) + nlp_length(shoot)

    # And it shows itself by the phases it joins, which is what a reader needs from it.
    @test occursin("brachistochrone", sprint(show, lk))
    @test occursin("earth_mars", sprint(show, lk))
end

@testset "OC manager — the assembled Jacobian is the derivative it claims to be" begin
    coll, shoot = _oc_phases()
    joined(c1, c2) = [c1.tf - c2.t0, c1.yf[3] - c2.y0[4]]
    lk = _OC.ContinuityLink(coll, shoot, joined, [0.0, 0.0], [0.0, 0.0], "join")
    om = _OC.OCManager(Any[coll, shoot], [lk])

    x = get_decision_vector(om)
    J = get_jacobian(om)
    @test size(J) == (om.n_funs, om.n_vars)
    @test all(isfinite, J)

    # Per-column steps, because this decision vector mixes throttles of order one with epochs of
    # order 1e7 seconds — the same span that defeats a single absolute step in the MGAnDSMs file.
    J_fd = zeros(om.n_funs, om.n_vars)
    for j in eachindex(x)
        h  = max(1e-6, 1e-7 * abs(x[j]))
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(om, xp); Fp = copy(get_functions(om))
        set_decision_vector!(om, xm); Fm = copy(get_functions(om))
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(om, x)

    worst = 0.0
    for j in axes(J, 2), i in axes(J, 1)
        denom = max(abs(J[i, j]), abs(J_fd[i, j]), 1e-8)
        worst = max(worst, abs(J[i, j] - J_fd[i, j]) / denom)
    end
    @test worst < 1e-3

    # The link rows in particular, which are what no single-phase test reaches.
    @test maximum(abs, J[end-1:end, :] .- J_fd[end-1:end, :]) < 1e-3
end

@testset "OC manager — the sparsity pattern claims every nonzero" begin
    coll, shoot = _oc_phases()
    joined(c1, c2) = [c1.tf - c2.t0]
    lk = _OC.ContinuityLink(coll, shoot, joined, [0.0], [0.0], "epochs")
    om = _OC.OCManager(Any[coll, shoot], [lk])

    rows, cols, vals = build_sparsity_pattern(om)
    @test length(rows) == length(cols) == length(vals)
    @test all(1 .<= rows .<= om.n_funs)
    @test all(1 .<= cols .<= om.n_vars)

    refilled = zeros(Float64, length(rows))
    get_jacobian_values!(refilled, om, rows, cols)
    @test refilled ≈ vals atol = 0.0

    # An entry the pattern calls zero and the derivative calls nonzero never reaches the solver,
    # and the solver finds its way anyway, so nothing fails. That is what this checks for.
    J = get_jacobian(om)
    claimed = Set(zip(rows, cols))
    unclaimed = count((i, j) for j in 1:om.n_vars, i in 1:om.n_funs
                      if abs(J[i, j]) > 0.0 && (i, j) ∉ claimed)
    @test unclaimed == 0
end

@testset "OC manager — the objective and its gradient come off the phases" begin
    coll, shoot = _oc_phases()
    om = _OC.OCManager(Any[coll, shoot])

    obj = get_objective(om)
    @test isfinite(obj)

    g = get_objective_gradient(om)
    @test length(g) == om.n_vars
    @test all(isfinite, g)

    # The gradient is a derivative, checked where it is not structurally zero.
    x = get_decision_vector(om)
    nz = findall(!=(0.0), g)
    @test !isempty(nz)
    for j in nz[1:min(3, length(nz))]
        h  = max(1e-6, 1e-7 * abs(x[j]))
        xp = copy(x); xp[j] += h
        xm = copy(x); xm[j] -= h
        set_decision_vector!(om, xp); op = get_objective(om)
        set_decision_vector!(om, xm); omv = get_objective(om)
        set_decision_vector!(om, x)
        @test g[j] ≈ (op - omv) / (2h) rtol = 1e-3
    end
end

@testset "OC manager — it refuses what it cannot assemble" begin
    # The exception type and its message.
    @test_throws ArgumentError _OC.OCManager(Any[])
    msg = try; _OC.OCManager(Any[]); catch e; sprint(showerror, e); end
    @test occursin("at least one phase", msg)

    coll, shoot = _oc_phases()
    om = _OC.OCManager(Any[coll, shoot])
    @test_throws ArgumentError set_decision_vector!(om, zeros(om.n_vars + 3))
    msg2 = try; set_decision_vector!(om, zeros(om.n_vars + 3))
           catch e; sprint(showerror, e); end
    @test occursin("one entry per NLP variable", msg2)
    @test occursin(string(om.n_vars), msg2)

    # A continuity link between phases of different state dimension is refused with both
    # dimensions named, rather than producing a mis-sized constraint.
    #
    # This used to raise a FieldError listing SimsFlanaganPhase's internal fields, because
    # `_state_dim` read `p.n_states` and no shooting phase declares that field — so the default
    # continuity path had never worked for a mixed sequence, which is the one thing OCManager is
    # for. It now reads the dimension off the unified boundary context, which is the same number
    # the link itself compares.
    seq = Sequence(coll)
    @test_throws ArgumentError add_continuity!(seq, coll, shoot)
    msg3 = try; add_continuity!(seq, coll, shoot); catch e; sprint(showerror, e); end
    @test occursin("same number of states", msg3)
    @test occursin("3", msg3) && occursin("7", msg3)

    # And the dimension it reports is the one the unified context actually produces, for both
    # kinds of phase, which is what keeps the check and the link in agreement.
    @test _OC._state_dim(coll)  == length(_OC._unified_boundary_context(coll).y0)
    @test _OC._state_dim(shoot) == length(_OC._unified_boundary_context(shoot).y0)
    @test _OC._state_dim(shoot) == 7

    # And mismatched bound lengths on a hand-written link.
    @test_throws ArgumentError add_continuity!((c1, c2) -> [0.0], seq, coll, shoot;
                                               lower_bounds = [0.0, 0.0], upper_bounds = [0.0])
end
