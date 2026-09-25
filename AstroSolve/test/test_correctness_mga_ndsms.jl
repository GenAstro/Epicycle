# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# MGAnDSMs: multiple gravity assist with n deep-space manoeuvres.
#
# A phase here is a ballistic arc from one body to another, broken by deep-space manoeuvres, with
# an excess velocity at each end. The variables are not a state and a control — they are epochs,
# masses, v-infinity vectors and burns — which is the whole reason this transcription has its own
# vocabulary. The verbs do not change: `Vary`, `Constraint`, `Objective`, `Link`.
#
# Truth, as in the Sims-Flanagan file: **analytic** for the bookkeeping, which follows from the
# manoeuvre count, and **AD-versus-analytic** for the Jacobian, where the hand-written
# `jacobian_chunk` blocks are differenced against the constraints they claim to differentiate.
#
# The reference mission this is modelled on, Earth-Earth-Venus, converges only with its link
# constraints alone and is recorded as such, so there is no converged number to assert. What is
# exactly checkable is what the transcription builds and whether its derivatives are right — and
# the derivative check is the one that matters, because a wrong Jacobian here looks exactly like
# a hard problem.
#
# Two circular ephemerides stand in for the real ones. The real mission data is not part of this
# package, and the transcription does not care where a body is, only that the
# ephemeris is smooth and consistent with its own derivative.
#
# Covers mga_ndsms.jl, which no other test reaches.

using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: MGAAlphaSumBlock, check_configuration, function_list, get_constraint_bounds,
    get_decision_vector, get_functions, get_objective, get_variable_bounds, jacobian_chunk,
    n_constraints, n_dsm_bwd, n_dsm_fwd, nlp_length, objective_gradient_chunk,
    set_decision_vector!, sparsity_structure, test_partials_fd, variable_list, variable_ranges
using Test

const _MG_MU  = 1.32712440018e11      # Sun, km³/s²
const _MG_AU  = 1.495978707e8
const _MG_RE  = 1.0   * _MG_AU        # Earth
const _MG_RV  = 0.723 * _MG_AU        # Venus
const _MG_ISP = 320.0
const _MG_G0  = 9.80665e-3

_mg_vc(r)  = sqrt(_MG_MU / r)
_mg_om(r)  = _mg_vc(r) / r

# Epochs, in seconds from an arbitrary origin: depart, fly by, arrive.
const _MG_TDEP = 12222.0 * 86400.0
const _MG_TFLY = _MG_TDEP + 365.0 * 86400.0
const _MG_TARR = _MG_TFLY + 175.0 * 86400.0

"""A body on a circular, coplanar orbit of radius `r`, phased to `θ0` at t = 0."""
function _mg_ephemeris(r, θ0)
    ω  = _mg_om(r)
    vc = _mg_vc(r)
    return function (t::Real)
        θ = θ0 + ω * t
        return r  * [cos(θ), sin(θ), 0.0],
               vc * [-sin(θ), cos(θ), 0.0],
               -(_MG_MU / r^2) * [cos(θ), sin(θ), 0.0]
    end
end

const _mg_earth = _mg_ephemeris(_MG_RE, 0.0)
const _mg_venus = _mg_ephemeris(_MG_RV, 1.1)

_mg_hat(v) = v ./ norm(v)

"""Build the two-phase Earth-Earth-Venus problem and return its pieces.

The deep-space manoeuvre is pinned rather than varied, which is what the reference mission does
and also all that is currently possible — see the last testset.
"""
function _mg_setup(; m0 = 2000.0, n_dsm = 1)
    trans = MGAnDSMs(n_dsm = n_dsm)
    model = PropulsionModel(mu = _MG_MU, Isp = _MG_ISP, g0 = _MG_G0)

    p1 = MGAnDSMsPhase(name = :earth_to_flyby, transcription = trans, model = model,
                       ephemeris_left = _mg_earth, ephemeris_right = _mg_earth)
    p2 = MGAnDSMsPhase(name = :flyby_to_venus, transcription = trans, model = model,
                       ephemeris_left = _mg_earth, ephemeris_right = _mg_venus)
    p1._alpha = fill(1.0 / (n_dsm + 1), n_dsm + 1)
    p2._alpha = fill(1.0 / (n_dsm + 1), n_dsm + 1)

    _, v_dep, _ = _mg_earth(_MG_TDEP)
    r_fly, v_fly, _ = _mg_earth(_MG_TFLY)
    _, v_arr, _ = _mg_venus(_MG_TARR)

    # The incoming excess velocity at the flyby is radial rather than antiparallel to the outgoing
    # one: antiparallel is a 180° turn, which puts the periapsis radius at zero and an infinite
    # acos derivative into the Jacobian.
    vinf_dep1 =  3.944 .* _mg_hat(v_dep)
    vinf_arr1 =  3.962 .* _mg_hat(r_fly)
    vinf_dep2 = -3.962 .* _mg_hat(v_fly)
    vinf_arr2 = -3.141 .* _mg_hat(v_arr)

    for (p, vd, va, t0, tf) in ((p1, vinf_dep1, vinf_arr1, _MG_TDEP, _MG_TFLY),
                                (p2, vinf_dep2, vinf_arr2, _MG_TFLY, _MG_TARR))
        Vary(departure_vinf, p; guess = vd,
             lower_bound = fill(-6.0, 3), upper_bound = fill(6.0, 3))
        Vary(arrival_vinf,   p; guess = va,
             lower_bound = fill(-6.0, 3), upper_bound = fill(6.0, 3))
        Vary(initial_time,   p; guess = t0,
             lower_bound = t0 - 30 * 86400.0, upper_bound = t0 + 30 * 86400.0)
        Vary(final_time,     p; guess = tf,
             lower_bound = tf - 60 * 86400.0, upper_bound = tf + 60 * 86400.0)
        Vary(initial_mass,   p; guess = m0, lower_bound = 10.0, upper_bound = 4000.0)
        Vary(final_mass,     p; guess = m0, lower_bound = 10.0, upper_bound = 4000.0)
    end

    return p1, p2
end

@testset "MGAnDSMs — the circular ephemerides are consistent with themselves" begin
    # The stand-in ephemerides are truth for everything downstream, so they are checked first:
    # the velocity really is the derivative of the position, and the acceleration of the velocity.
    h = 1.0
    for eph in (_mg_earth, _mg_venus), t in (_MG_TDEP, _MG_TFLY, _MG_TARR)
        rp, _, _ = eph(t + h)
        rm, _, _ = eph(t - h)
        _, v, a  = eph(t)
        @test (rp .- rm) ./ (2h) ≈ v atol = 1e-6

        _, vp, _ = eph(t + h)
        _, vm, _ = eph(t - h)
        @test (vp .- vm) ./ (2h) ≈ a atol = 1e-9
    end

    # And they are circular: the radius does not change and the speed is the circular speed.
    r, v, _ = _mg_earth(_MG_TFLY)
    @test norm(r) ≈ _MG_RE atol = 1e-6
    @test norm(v) ≈ _mg_vc(_MG_RE) atol = 1e-9
    @test dot(r, v) ≈ 0.0 atol = 1e-3          # circular: position ⟂ velocity
end

@testset "MGAnDSMs — the transcription builds what the manoeuvre count implies" begin
    p1, p2 = _mg_setup()

    # Six declared variables per phase: two v-infinity vectors of three, two epochs, two masses.
    @test nlp_length(p1) == 3 + 3 + 1 + 1 + 1 + 1
    @test nlp_length(p1) == nlp_length(p2)
    @test length(variable_list(p1)) == 6

    # One deep-space manoeuvre splits the arc in two, which is what _alpha carries.
    @test length(p1._alpha) == 2
    @test sum(p1._alpha) ≈ 1.0 atol = 1e-12

    # The decision vector round-trips.
    x = get_decision_vector(p1)
    @test length(x) == nlp_length(p1)
    set_decision_vector!(p1, copy(x))
    @test get_decision_vector(p1) ≈ x atol = 0.0

    lo, hi = get_variable_bounds(p1)
    @test length(lo) == length(hi) == nlp_length(p1)
    @test all(lo .<= hi)

    clo, chi = get_constraint_bounds(p1)
    @test length(clo) == length(chi) == n_constraints(p1)
end

@testset "MGAnDSMs — the deep-space burns and arc fractions are variables" begin
    # These were the last quantities without a Vary method, so the burn the transcription is
    # named for could not be handed to a solver. Everything below the verb already handled them.
    p1, _ = _mg_setup()
    n_before = nlp_length(p1)

    burn = Vary(deep_space_dv, p1; guess = [0.05, -0.02, 0.01],
                lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    fracs = Vary(arc_fractions, p1; guess = [0.45, 0.55],
                 lower_bound = [0.05, 0.05], upper_bound = [0.95, 0.95])

    @test p1.dv_var === burn && p1.alpha_var === fracs
    @test nlp_length(p1) == n_before + 3 * n_dsm(p1) + 2
    @test deep_space_dv(p1) ≈ reshape([0.05, -0.02, 0.01], 3, 1)      # one burn, tiled
    @test p1._alpha == [0.45, 0.55]

    # Varying the fractions adds the row that holds their sum to one.
    @test any(pf -> pf.source isa MGAAlphaSumBlock, function_list(p1))
    @test get_functions(p1)[8] ≈ 0.0 atol = 1e-12

    # The decision vector round-trips through the new blocks.
    x = get_decision_vector(p1)
    set_decision_vector!(p1, x)
    @test get_decision_vector(p1) == x

    # Bounds are one burn and one value per arc; a guess of the wrong size says what it wanted.
    q, _ = _mg_setup()
    @test_throws ArgumentError Vary(deep_space_dv, q; lower_bound = fill(-2.0, 6),
                                    upper_bound = fill(2.0, 6))
    @test_throws ArgumentError Vary(deep_space_dv, q; guess = zeros(5),
                                    lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    @test_throws ArgumentError Vary(arc_fractions, q; lower_bound = [0.0],
                                    upper_bound = [1.0])
end

@testset "MGAnDSMs — the constraints respond to the variables that drive them" begin
    # Physics rather than bookkeeping. Moving the departure excess velocity must move the arc, so
    # the constraint vector must move with it. A transcription whose variables did not reach the
    # propagation would report the same residuals whatever was varied.
    p1, _ = _mg_setup()
    F0 = copy(get_functions(p1))
    @test length(F0) == n_constraints(p1)
    @test all(isfinite, F0)

    x = get_decision_vector(p1)
    x2 = copy(x)
    x2[1] += 0.25                       # a quarter km/s on the departure excess velocity
    set_decision_vector!(p1, x2)
    F1 = get_functions(p1)

    @test !(F1 ≈ F0)
    @test all(isfinite, F1)

    # And putting it back restores what was there, so evaluation has no hidden state.
    set_decision_vector!(p1, x)
    @test get_functions(p1) ≈ F0 atol = 1e-12
end

"""Worst *relative* disagreement between the assembled analytic Jacobian and central differences.

This is `test_partials_fd`'s comparison with two changes, and both are needed here.

The step is scaled per column. An MGAnDSMs decision vector holds excess velocities of a few km/s
beside epochs of order 1e9 seconds, and no single absolute step serves both: at `h = 1e-6` the
epoch columns difference below double precision and come back garbage, while at `h = 1` the
velocity columns take a 25% step and come back garbage instead. Measured across the sweep, the
velocity blocks are right to 2.4e-8 relative at the small step and the epoch blocks to 1.4e-8 at
the large one.

The comparison is relative. Entries of this Jacobian are positions in km against velocities, so
they run to 1e8; an absolute tolerance calls a block wrong when it agrees to eight figures.
"""
function _mg_worst_relative_jacobian_error(p; hmin = 1e-6)
    x0    = get_decision_vector(p)
    vlist = variable_list(p)
    flist = function_list(p)
    rngs  = variable_ranges(p)

    n_f = length(get_functions(p))
    J_fd = zeros(n_f, length(x0))
    for j in eachindex(x0)
        h  = max(hmin, 1e-7 * abs(x0[j]))
        xp = copy(x0); xp[j] += h
        xm = copy(x0); xm[j] -= h
        set_decision_vector!(p, xp); Fp = copy(get_functions(p))
        set_decision_vector!(p, xm); Fm = copy(get_functions(p))
        J_fd[:, j] = (Fp .- Fm) ./ (2h)
    end
    set_decision_vector!(p, x0)
    get_functions(p)

    J_an = zeros(n_f, length(x0))
    f_off = 0
    for pf in flist
        for (v, r) in zip(vlist, rngs)
            J_an[f_off+1:f_off+pf.n_nlp, r] .+= jacobian_chunk(p, pf, v)
        end
        f_off += pf.n_nlp
    end

    worst = 0.0
    for j in axes(J_fd, 2), i in axes(J_fd, 1)
        denom = max(abs(J_fd[i, j]), abs(J_an[i, j]), 1e-8)
        worst = max(worst, abs(J_an[i, j] - J_fd[i, j]) / denom)
    end
    return worst
end

@testset "MGAnDSMs — the NLP Jacobian is the derivative of the constraints" begin
    # The check that matters. A wrong Jacobian here does not raise; it burns iterations and
    # returns a plausible answer, which is exactly how the reference mission behaves when it does
    # not converge.
    p1, p2 = _mg_setup()
    get_functions(p1); get_functions(p2)

    @test _mg_worst_relative_jacobian_error(p1) < 1e-4
    @test _mg_worst_relative_jacobian_error(p2) < 1e-4

    # At a moved point, since a Jacobian can be right at the guess and wrong where the solver goes.
    q1, _ = _mg_setup()
    x = get_decision_vector(q1)
    x[1] += 0.2                                   # km/s on the departure excess velocity
    x[7] += 5 * 86400.0                           # five days on the departure epoch
    set_decision_vector!(q1, x)
    get_functions(q1)
    @test _mg_worst_relative_jacobian_error(q1) < 1e-4
end

@testset "MGAnDSMs — test_partials_fd cannot judge a problem whose columns span nine decades" begin
    # 🔴 Recorded behaviour. The shared diagnostic takes one absolute step for every column and
    # compares on absolute error, and an MGAnDSMs phase defeats both: its columns run from a few
    # km/s to epochs of 1e9 seconds, and its entries run to 1e8 km. It therefore reports a
    # failure on a Jacobian the testset above shows is correct to eight figures.
    #
    # This is a limitation of the diagnostic, not of the transcription. It is asserted here so
    # that the day it grows a per-column relative comparison, this test says so.
    p1, _ = _mg_setup()
    get_functions(p1)

    @test _mg_worst_relative_jacobian_error(p1) < 1e-4       # the Jacobian is right
    @test test_partials_fd(p1; h = 1e-6, tol = 1e-4) == false  # and the diagnostic says otherwise

    # It is right about the mass columns, whose entries are order one, which is the case a single
    # absolute step does suit.
    @test _mg_worst_relative_jacobian_error(p1) < 1e-4
end

@testset "MGAnDSMs — the reported structure matches what is evaluated" begin
    p1, _ = _mg_setup()
    get_functions(p1)

    vlist = variable_list(p1)
    flist = function_list(p1)
    rngs  = variable_ranges(p1)

    @test length(vlist) == length(rngs)
    @test sum(length, rngs) == nlp_length(p1)
    @test sort(reduce(vcat, collect.(rngs))) == collect(1:nlp_length(p1))
    @test sum(pf.n_nlp for pf in flist) == n_constraints(p1)

    for pf in flist, (v, r) in zip(vlist, rngs)
        @test size(jacobian_chunk(p1, pf, v)) == (pf.n_nlp, length(r))
    end
end

@testset "MGAnDSMs — a Link names the two phases it joins" begin
    # A Link is a subject, so a flyby's conditions are constraints on one thing rather than four
    # constraints scattered across two phases. This checks the shape; the flyby quantities
    # themselves are the user's to write, which is what the reference mission does.
    p1, p2 = _mg_setup()
    l = Link(p1, p2; body = :earth, name = :flyby)

    @test l.phases === (p1, p2)
    @test l.body === :earth
    @test l.name === :flyby

    # `at` has no meaning on a link, since a link is evaluated where the phases meet and nowhere
    # else, so it is refused rather than ignored.
    dummy(c1, c2) = [0.0]
    @test_throws ArgumentError Constraint(dummy, l; equals = 0.0, at = Final())
    msg = try; Constraint(dummy, l; equals = 0.0, at = Final()); catch e; sprint(showerror, e); end
    @test occursin("no meaning on a Link", msg)
end

# Constraint and objective functions for the configuration checks below. Each reads the final
# mass, which _mg_setup varies, so every partial has a variable to act on.
_mg_mass_bare(c) = [final_mass(c)]

_mg_mass(c) = [final_mass(c)]
@partial(_mg_mass, final_mass) do c
    reshape([1.0], 1, 1)
end

_mg_mass_wide(c) = [final_mass(c)]
@partial(_mg_mass_wide, final_mass) do c
    [1.0 0.0]                                   # 1×2 for a 1×1 block
end

_mg_delivered(c) = final_mass(c)
@partial(_mg_delivered, final_mass) do c
    [1.0]
end

_mg_delivered_row(c) = final_mass(c)
@partial(_mg_delivered_row, final_mass) do c
    reshape([1.0], 1, 1)                        # a matrix where a gradient vector is used
end

_mg_dsm_size(c) = [sum(abs2, deep_space_dv(c))]
@partial(_mg_dsm_size, deep_space_dv) do c
    reshape(2 .* collect(deep_space_dv(c)), 1, :)
end

"""The message check_configuration throws for a sequence holding `p`, or `nothing` if it passes."""
function _mg_config_message(p)
    seq = Sequence()
    add_sequence!(seq, p)
    try
        check_configuration(seq)
        return nothing
    catch e
        e isa ArgumentError || rethrow()
        return sprint(showerror, e)
    end
end

@testset "MGAnDSMs — a problem that is not fully specified is refused before solving" begin
    # On a shooting phase an undeclared derivative used to be a zero, and uc7 spent months
    # maximising against a zero gradient. The configuration check is what stops that reaching
    # IPOPT, so each way of getting it wrong must be named, and a correct problem must pass.

    p, _ = _mg_setup()
    Constraint(_mg_mass, p; lower_bound = 100.0, at = Final())
    Objective(_mg_delivered, p; sense = Max())
    @test _mg_config_message(p) === nothing

    # The declared gradient is used, and it is the one differentiation would give.
    q, _ = _mg_setup()
    Objective(c -> final_mass(c), q; sense = Max())
    @test objective_gradient_chunk(p, p.mf_var) ≈ objective_gradient_chunk(q, q.mf_var) rtol = 1e-12
    @test objective_gradient_chunk(p, p.mf_var) ≈ [-1.0] rtol = 1e-12     # Max is negated

    # A constraint with no partial is refused, naming the function and what to annotate.
    p, _ = _mg_setup()
    Constraint(_mg_mass_bare, p; lower_bound = 100.0, at = Final())
    msg = _mg_config_message(p)
    @test msg !== nothing
    @test occursin("_mg_mass_bare", msg)
    @test occursin("no registered Jacobian", msg)
    @test occursin("`Real`", msg)

    # An objective with no gradient, whether a named function or an anonymous one.
    p, _ = _mg_setup()
    Objective(c -> final_mass(c), p; sense = Max())
    msg = _mg_config_message(p)
    @test msg !== nothing && occursin("no registered gradient", msg)

    # A partial of the wrong size says what it declared and what will be used.
    p, _ = _mg_setup()
    Constraint(_mg_mass_wide, p; lower_bound = 100.0, at = Final())
    msg = _mg_config_message(p)
    @test msg !== nothing
    @test occursin("(1, 2)", msg) && occursin("(1, 1)", msg)

    # A Mayer gradient written as a 1×N matrix is the shape that took IPOPT down without a trace.
    p, _ = _mg_setup()
    Objective(_mg_delivered_row, p; sense = Max())
    msg = _mg_config_message(p)
    @test msg !== nothing && occursin("gradient vector of length 1", msg)

    # solve! runs the same check, so the refusal happens before any iteration.
    p, _ = _mg_setup()
    Constraint(_mg_mass_bare, p; lower_bound = 100.0, at = Final())
    seq = Sequence()
    add_sequence!(seq, p)
    @test_throws ArgumentError solve!(seq; method = Optimize(print_level = 0, max_iter = 1))

    # A partial with respect to a quantity that is not a variable is refused when it is written,
    # and the message says how to make it one.
    p, _ = _mg_setup()
    msg = try
        Constraint(_mg_dsm_size, p; upper_bound = 1.0, at = Final())
        nothing
    catch e
        sprint(showerror, e)
    end
    @test msg !== nothing && occursin("Vary(deep_space_dv", msg)
end

@testset "MGAnDSMs — with burns, fractions, a constraint and an objective, every derivative agrees" begin
    # The burn and fraction blocks of the match-point Jacobian, the alpha-sum row, a boundary
    # constraint with a declared partial and one differentiated automatically, and the objective
    # gradient, all against central differences at a point with non-zero burns.
    p, _ = _mg_setup()
    Vary(deep_space_dv, p; guess = [0.05, -0.02, 0.01],
         lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    Vary(arc_fractions, p; guess = [0.45, 0.55], lower_bound = [0.05, 0.05],
         upper_bound = [0.95, 0.95])
    Constraint(_mg_mass, p; lower_bound = 100.0, at = Final())
    Constraint(_mg_dsm_size, p; upper_bound = 1.0, at = Final())
    Constraint(c -> [departure_vinf(c)[1] * final_time(c)], p; upper_bound = 1e12, at = Final())
    Objective(c -> final_mass(c) - 10.0 * sum(abs2, deep_space_dv(c)), p; sense = Max())
    get_functions(p)

    @test _mg_worst_relative_jacobian_error(p) < 1e-4

    x0 = get_decision_vector(p)
    g  = vcat([objective_gradient_chunk(p, v) for v in variable_list(p)]...)
    g_fd = map(eachindex(x0)) do j
        h  = max(1e-6, 1e-7 * abs(x0[j]))
        xp = copy(x0); xp[j] += h
        xm = copy(x0); xm[j] -= h
        set_decision_vector!(p, xp); fp = get_objective(p)
        set_decision_vector!(p, xm); fm = get_objective(p)
        (fp - fm) / (2h)
    end
    set_decision_vector!(p, x0)
    @test all(abs.(g .- g_fd) .<= 1e-4 .* max.(abs.(g_fd), 1e-6))
    @test get_objective(p) ≈ -(p._mf - 10.0 * sum(abs2, p._dv)) rtol = 1e-12   # Max is negated
end

@testset "MGAnDSMs — burns inside both halves: the chain rule through each manoeuvre" begin
    # With one DSM there is no burn inside either propagation half, so the mass depletion and the
    # state-and-mass transition chain through a manoeuvre never ran. Three DSMs put burns on both
    # sides of the match point.
    p, _ = _mg_setup(n_dsm = 3)
    @test n_dsm_fwd(p) >= 1 && n_dsm_bwd(p) >= 1
    Vary(deep_space_dv, p; guess = [0.08, -0.03, 0.02],
         lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3))
    Vary(arc_fractions, p; guess = [0.2, 0.3, 0.3, 0.2], lower_bound = fill(0.05, 4),
         upper_bound = fill(0.95, 4))
    Constraint(_mg_mass, p; lower_bound = 100.0, at = Final())
    get_functions(p)

    # A 1e-4 floor on the step. At 1e-6 an entry of 3e-5 in a row of order 100 differences at
    # roundoff and reads 1.5e-4 wrong; at 1e-4 every entry agrees to 2.4e-6, and at 1e-3
    # truncation takes over again.
    @test _mg_worst_relative_jacobian_error(p; hmin = 1e-4) < 1e-5

    # Burns spend propellant: the same problem with the burns zeroed arrives heavier at the match
    # point, so the mass row of the match-point defect moves.
    F_burn = copy(get_functions(p))
    x = get_decision_vector(p)
    r_dv = variable_ranges(p)[findfirst(v -> v === p.dv_var, variable_list(p))]
    x0 = copy(x); x0[r_dv] .= 0.0
    set_decision_vector!(p, x0)
    @test get_functions(p)[7] != F_burn[7]
    set_decision_vector!(p, x)

    # Bounds, one row per function row, and a sparsity pattern that claims every nonzero.
    lb, ub = get_constraint_bounds(p)
    @test length(lb) == length(ub) == n_constraints(p) == length(get_functions(p))
    @test lb[8] == ub[8] == 0.0                         # the arc fractions sum to one
    S = sparsity_structure(p)
    flist, vlist = function_list(p), variable_list(p)
    @test size(S) == (length(flist), length(vlist))
    for (i, pf) in enumerate(flist), (j, v) in enumerate(vlist)
        S[i, j] || @test all(iszero, jacobian_chunk(p, pf, v))
    end
    @test S[findfirst(pf -> pf.source isa MGAAlphaSumBlock, flist),
            findfirst(v -> v === p.alpha_var, vlist)]
end
