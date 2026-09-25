# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Sims-Flanagan transcription: assembly, physics and derivatives, without a solve.
#
# The low-thrust Earth-to-Mars benchmark this is built from does not converge — it reaches
# IPOPT's iteration cap at a final mass of 1153.5848 kg, which the regression suite records as a
# baseline to match rather than as an optimum. Running it here would cost minutes and would assert
# a number that is not the answer to anything.
#
# So this file tests the parts that are exactly checkable and fast: what the transcription builds,
# what the physics does at a stated point, and whether the Jacobian the solver is handed is the
# derivative of the constraints it is handed. That last one is the check that matters most,
# because a wrong NLP Jacobian does not fail — it burns iterations and returns something
# plausible, which is exactly what the unconverged benchmark looks like.
#
# Truth is a mix, and each testset says which it is using:
#
#   · **analytic** for the Hohmann reference the problem is posed against, and for the segment
#     bookkeeping, both of which are closed-form;
#   · **AD-versus-analytic**, in the §2 sense, for the Jacobian: the hand-written
#     `jacobian_chunk` blocks against central differences of the same functions.
#
# Covers sims_flanagan.jl and the shooting half of diagnostics.jl, neither of which any other
# open-tier test reaches.

using SNOW
using LinearAlgebra
using EpicycleBase
using AstroSolve
using AstroSolve: MayerObjective, dump_partials_fd, function_list, get_constraint_bounds,
    get_decision_vector, get_functions, get_variable_bounds, jacobian_chunk, matchpoint_defect,
    n_constraints, n_fwd, nlp_length, report_functions, report_variables, set_decision_vector!,
    set_matchpoint_scale!, set_objective!, test_partials_fd, variable_list, variable_ranges
using Test

# ── The mission, in km, kg and seconds ───────────────────────────────────────

const _SF_MU   = 1.32712440018e11     # Sun, km³/s²
const _SF_AU   = 1.495978707e8        # km
const _SF_RE   = 1.0   * _SF_AU
const _SF_RM   = 1.524 * _SF_AU

const _SF_VE   = sqrt(_SF_MU / _SF_RE)
const _SF_VM   = sqrt(_SF_MU / _SF_RM)
const _SF_OME  = _SF_VE / _SF_RE
const _SF_OMM  = _SF_VM / _SF_RM

const _SF_M0   = 1500.0               # kg
const _SF_ISP  = 3000.0               # s
const _SF_TMAX = 1.0e-3               # kN, i.e. 1 N
const _SF_G0   = 9.80665e-3           # km/s²

# Hohmann transfer: half the period of the transfer ellipse.
const _SF_A    = (_SF_RE + _SF_RM) / 2
const _SF_TOF  = π * sqrt(_SF_A^3 / _SF_MU)

# Mars is placed so that it arrives where the transfer ellipse does.
const _SF_THM0 = π - _SF_OMM * _SF_TOF

function _sf_earth(t::Real)
    θ = _SF_OME * t
    return _SF_RE * [cos(θ), sin(θ), 0.0],
           _SF_VE * [-sin(θ), cos(θ), 0.0],
           -(_SF_MU / _SF_RE^2) * [cos(θ), sin(θ), 0.0]
end

function _sf_mars(t::Real)
    θ = _SF_THM0 + _SF_OMM * t
    return _SF_RM * [cos(θ), sin(θ), 0.0],
           _SF_VM * [-sin(θ), cos(θ), 0.0],
           -(_SF_MU / _SF_RM^2) * [cos(θ), sin(θ), 0.0]
end

const _SF_NSEG = 20        # fewer than the benchmark's 60: this file never solves
const _SF_EPS  = 0.001     # the regularization weight the benchmark uses

_sf_thrust_ball(c) = [dot(control(c), control(c))]

"""Build the phase, declare its variables, and seed it with the prograde guess.

Returns `(phase, n_half, vars)`. Each call builds fresh objects, since the spec verbs mutate the
phase they are given and the decision vector is written back onto it.
"""
function _sf_phase(; n_segments = _SF_NSEG, throttle = 0.5)
    trans = SimsFlanagan(n_segments = n_segments)
    nh    = n_fwd(trans)

    phase = SimsFlanaganPhase(
        name            = :earth_mars,
        transcription   = trans,
        model           = PropulsionModel(mu = _SF_MU, Isp = _SF_ISP,
                                          Tmax = _SF_TMAX, g0 = _SF_G0),
        ephemeris_left  = _sf_earth,
        ephemeris_right = _sf_mars)

    uf = Vary(forward_control,  phase;
              lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3), name = "u_fwd")
    ub = Vary(backward_control, phase;
              lower_bound = fill(-2.0, 3), upper_bound = fill(2.0, 3), name = "u_bwd")
    mf = Vary(final_mass, phase;
              lower_bound = 100.0, upper_bound = _SF_M0, scale = _SF_M0, name = "mf")

    set_matchpoint_scale!(phase, [_SF_RE, _SF_RE, _SF_RE, _SF_VE, _SF_VE, _SF_VE, _SF_M0])

    # A ball, not the box `Vary` gives: a box admits √3 at its corners.
    Constraint(_sf_thrust_ball, phase; lower_bound = 0.0, upper_bound = 1.0,
               at = Path(), name = :thrust_ball)

    obj = MayerObjective(phase; sense = :Max) do ctx
        ctx.mf - _SF_EPS * (sum(dot(phase._u_fwd[:, k], phase._u_fwd[:, k]) for k in 1:nh) +
                            sum(dot(phase._u_bwd[:, k], phase._u_bwd[:, k]) for k in 1:nh))
    end
    set_objective!(phase, obj)
    AstroSolve.add_objective_jacobian!(() -> [1.0], obj, mf)
    AstroSolve.add_objective_jacobian!(() -> -2.0 * _SF_EPS .* vec(phase._u_fwd), obj, uf)
    AstroSolve.add_objective_jacobian!(() -> -2.0 * _SF_EPS .* vec(phase._u_bwd), obj, ub)

    # Pinned rather than varied: the epochs and the wet mass are given, and v∞ is zero at both
    # ends, so none of them is an NLP variable.
    phase._t0 = 0.0
    phase._tf = _SF_TOF
    phase._m0 = _SF_M0

    # The standard seed for an energy-raising transfer: thrust prograde on both halves.
    _, v0, _ = _sf_earth(0.0)
    _, vf, _ = _sf_mars(_SF_TOF)
    for k in 1:nh
        phase._u_fwd[:, k] = throttle * v0 / norm(v0)
        phase._u_bwd[:, k] = throttle * vf / norm(vf)
    end
    phase._mf = _SF_M0 * 0.85

    return phase, nh, (uf, ub, mf)
end

@testset "Sims-Flanagan — the Hohmann reference the problem is posed against" begin
    # Analytic, and worth checking before anything is asserted relative to it: the impulsive
    # two-burn transfer is the upper bound a low-thrust solution is measured against, and it is
    # Tsiolkovsky applied to two closed-form burns.
    dv1 = _SF_VE * (sqrt(2 * _SF_RM / (_SF_RE + _SF_RM)) - 1)
    dv2 = _SF_VM * (1 - sqrt(2 * _SF_RE / (_SF_RE + _SF_RM)))

    # The standard figures for an Earth-to-Mars Hohmann transfer, for a Mars radius of 1.524 AU.
    @test dv1 ≈ 2.9461 atol = 1e-3
    @test dv2 ≈ 2.6500 atol = 1e-3
    @test dv1 + dv2 ≈ 5.5961 atol = 2e-3
    @test _SF_TOF / 86400 ≈ 258.9 atol = 0.2

    # Tsiolkovsky on the two burns. The shipped example's comment calls this bound 1249 kg, but
    # its own formula on its own constants gives 1240.2 — a stale comment rather than a different
    # calculation, and recorded in the test plan.
    mf = _SF_M0 * exp(-(dv1 + dv2) / (_SF_ISP * _SF_G0))
    @test mf ≈ 1240.175 atol = 0.05

    # A low-thrust transfer cannot beat the impulsive bound, which is what makes this the ceiling
    # the recorded benchmark of 1153.58 kg sits under.
    @test 1153.5848 < mf
end

@testset "Sims-Flanagan — the transcription builds what the segments imply" begin
    # Analytic bookkeeping. Every count here follows from the segment count, so a change in the
    # assembly shows up as an arithmetic disagreement rather than as a solver that behaves oddly.
    phase, nh, _ = _sf_phase()

    @test nh == _SF_NSEG ÷ 2

    # Three throttle components per segment on each half, plus the final mass.
    @test nlp_length(phase) == 3 * nh + 3 * nh + 1

    # Seven match-point defects — three position, three velocity, one mass — and one ball
    # constraint per segment.
    @test n_constraints(phase) == 7 + 2 * nh

    # The decision vector round-trips through the phase unchanged.
    x0 = get_decision_vector(phase)
    @test length(x0) == nlp_length(phase)
    set_decision_vector!(phase, copy(x0))
    @test get_decision_vector(phase) ≈ x0 atol = 0.0

    # Bounds are the ones that were declared, in the same order.
    lo, hi = get_variable_bounds(phase)
    @test length(lo) == length(hi) == nlp_length(phase)
    @test all(lo[1:6nh] .== -2.0)
    @test all(hi[1:6nh] .==  2.0)
    @test lo[end] == 100.0 / _SF_M0      # scaled by the `scale = _SF_M0` on the variable
    @test hi[end] == 1.0

    # The constraint bounds: match-point defects are equalities, the ball is an inequality.
    clo, chi = get_constraint_bounds(phase)
    @test length(clo) == length(chi) == n_constraints(phase)
    @test all(clo[1:7] .== 0.0) && all(chi[1:7] .== 0.0)
    @test all(chi[8:end] .== 1.0)
end

@testset "Sims-Flanagan — the match point responds to the throttle" begin
    # Physics rather than bookkeeping. The seed does not close the transfer, so the defect is
    # large; thrusting harder must change it. A transcription whose controls did not reach the
    # propagation would report the same defect whatever the throttle, and would then "converge"
    # by driving the mass instead.
    coast, _, _ = _sf_phase(throttle = 0.0)
    get_functions(coast)
    d_coast = matchpoint_defect(coast)

    push_, _, _ = _sf_phase(throttle = 1.0)
    get_functions(push_)
    d_push = matchpoint_defect(push_)

    @test length(d_coast) == 7
    @test length(d_push) == 7

    # The two differ, and by a margin far outside round-off.
    @test norm(d_coast[1:3] .- d_push[1:3]) > 1.0          # km
    @test norm(d_coast[4:6] .- d_push[4:6]) > 1e-4         # km/s

    # Thrusting spends propellant, so the mass defect moves too.
    @test d_coast[7] != d_push[7]

    # And the seed really is far from closed, which is why the benchmark needs a solver.
    @test norm(d_coast[1:3]) > 1e5
end

@testset "Sims-Flanagan — the NLP Jacobian is the derivative of the constraints" begin
    # AD-versus-analytic, against central differences. This
    # is the check that matters most for a shooting transcription: a wrong Jacobian does not
    # raise, it burns iterations and returns a plausible answer, which is indistinguishable from
    # a hard problem. `test_partials_fd` assembles every `jacobian_chunk` block and compares.
    phase, _, _ = _sf_phase()
    get_functions(phase)          # warm up the match-point evaluation

    @test test_partials_fd(phase; h = 1e-6, tol = 1e-4)

    # The same at a different point in the decision space, since a Jacobian can be right at the
    # seed and wrong where the solver actually goes.
    moved, _, _ = _sf_phase(throttle = 0.9)
    get_functions(moved)
    @test test_partials_fd(moved; h = 1e-6, tol = 1e-4)
end

@testset "Sims-Flanagan — the reported structure matches what is evaluated" begin
    # `variable_list`, `function_list` and `variable_ranges` are what the manager uses to place
    # each Jacobian block, and `test_partials_fd` above trusts them. Checking them directly means
    # a mismatch is reported as a structural disagreement rather than as a derivative error.
    phase, nh, _ = _sf_phase()

    # `jacobian_chunk` reads the propagated arcs, which only exist once the constraints have been
    # evaluated once. Without this warm-up it fails from inside the Kepler propagator with
    # "the initial position must have non-zero magnitude", which is a §9.9 boundary-guard gap
    # rather than a wrong answer; recorded in the test plan.
    get_functions(phase)

    vlist = variable_list(phase)
    flist = function_list(phase)
    rngs  = variable_ranges(phase)

    @test length(vlist) == length(rngs) == 3            # u_fwd, u_bwd, mf
    @test sum(length, rngs) == nlp_length(phase)
    @test isempty(intersect(rngs[1], rngs[2]))          # blocks do not overlap
    @test sum(pf.n_nlp for pf in flist) == n_constraints(phase)

    # Every declared block has the shape the assembly will place it into.
    for pf in flist, (v, r) in zip(vlist, rngs)
        chunk = jacobian_chunk(phase, pf, v)
        @test size(chunk) == (pf.n_nlp, length(r))
    end

    # And the functions evaluate to the length they advertise.
    F = get_functions(phase)
    @test length(F) == n_constraints(phase)
    @test all(isfinite, F)
end

"""Run `f` and return what it printed."""
function _sf_capture(f)
    mktemp() do path, io
        redirect_stdout(f, io)
        close(io)
        read(path, String)
    end
end

@testset "Sims-Flanagan — the diagnostic tables report what the phase holds" begin
    # report_variables, report_functions and dump_partials_fd are what a user reaches for when a
    # solve misbehaves, so what they print has to agree with the phase and leave it untouched.
    phase, _, _ = _sf_phase()
    get_functions(phase)
    x0 = copy(get_decision_vector(phase))
    n_x, n_f = nlp_length(phase), n_constraints(phase)

    # One row per scalar variable, and a flag on any at a bound.
    vars = _sf_capture(() -> report_variables(phase))
    rows = filter(l -> startswith(l, "  ") && !occursin(r"\bname\b", l), split(vars, '\n'))
    @test length(rows) == n_x
    for v in variable_list(phase)
        @test occursin(v.name, vars)
    end

    # One row per constraint, with the violation flagged where the seed leaves the match point
    # open, which is every match-point row at the seed.
    funcs = _sf_capture(() -> report_functions(phase))
    frows = filter(l -> startswith(l, "  ") && !occursin(r"\bname\b", l), split(funcs, '\n'))
    @test length(frows) == n_f
    @test occursin("← viol", funcs)

    # The full dump: every Jacobian entry, then the objective gradient, each with its error. The
    # errors it prints must agree with test_partials_fd, which passed above at tol = 1e-4.
    dump = _sf_capture(() -> dump_partials_fd(phase))
    jac_part, grad_part = split(dump, "Objective gradient dump")
    number = r"[-+]?\d+(\.\d+)?([eE][-+]?\d+)?"
    entry_rows(block) = filter(l -> occursin(number, l) && !occursin("analytic", l) &&
                                    !occursin("h=", l), split(block, '\n'))
    jrows = entry_rows(jac_part)
    grows = entry_rows(grad_part)
    @test length(jrows) == n_f * n_x
    @test length(grows) == n_x
    abs_err(l) = parse(Float64, split(strip(l))[end - 1])
    @test maximum(abs_err, grows) < 1e-4

    # Differencing moves the decision vector; the dump must put it back.
    @test get_decision_vector(phase) == x0
end
