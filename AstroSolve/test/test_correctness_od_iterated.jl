# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Iterated EKF + RTS smoother, and the model-based EKF entry point.
#
# `run_iterated_rts!` is Gauss-Newton on the MAP problem: filter forward, smooth backward,
# back-propagate the smoothed estimate at the first observation to the reference epoch, and
# repeat from there. Its reason to exist is a poor initial linearization, which is the case a
# plain EKF handles worst, so that is the case tested here.
#
# Three properties are exact and are what the assertions rest on. One sweep must reproduce
# `run_ekf!` followed by `run_rts` on the same data, which the docstring claims outright. The
# iterated estimate at the reference epoch and a batch least-squares fit of the same data are
# both MAP estimates of the same initial state under the same prior, so they must agree — a
# cross-method check, and independent of the filter
# tests because the batch path shares no code with the smoother. And a vector covariance and the
# matrix with that diagonal state the same prior, so they must give the same answer.
#
# This file depends on the fixture in test_correctness_od_batch.jl, which runtests.jl includes
# first, the same way test_correctness_rts.jl depends on test_correctness_ekf.jl.

using AstroSolve
using AstroSolve: run_ekf!, run_iterated_rts!, run_rts
using AstroModels
using AstroProp
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase
using LinearAlgebra
using Test

const _ORTS_EKF = AstroSolve.ExtendedKalmanFilter

# Half an hour of range and Doppler. Shorter than the batch fixture's full orbit because each
# outer sweep is a complete filter pass with propagation, and the checks here are either
# exact identities or comparisons between two estimators on the same data, neither of which
# needs the observability that the absolute recovery test in the batch file needs.
const _ORTS_OFFSETS = collect(60.0:60.0:1800.0)

# A prior wide enough not to drive the answer, in the shapes the estimators accept.
const _ORTS_COV_VEC = [1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]
const _ORTS_COV_MAT = Diagonal(_ORTS_COV_VEC)

"""Build a spacecraft at `guess`, the variable that estimates its state, and the problem.

Every call makes fresh objects. `run_iterated_rts!` and `run_ekf!` both write the estimate back
onto the spacecraft through `assign!`, so two runs cannot share one.
"""
function _orts_setup(guess; covariance = _ORTS_COV_MAT)
    gs = _od_station()
    sc = _od_sc(guess)
    y0 = Vary(state, sc; guess = guess, covariance = covariance)
    problem = ODProblem(spacecraft = sc, propagator = _od_prop(),
                        measurements = _od_specs(gs, sc),
                        solve_for = [y0])
    return sc, y0, problem
end

# A starting error large enough that the first sweep's linearization is visibly wrong, which is
# the condition the iteration exists for.
const _ORTS_ERR = [3.0, -2.0, 1.5, 2e-3, -3e-3, 1e-3]

# The two results this file checks belong to different epochs, and comparing either against the
# wrong one is off by a whole orbit's worth of motion — 7.66 km/s times the 60 s to the first
# observation is 460 km, which looks like divergence and is not.
#
#   `rts.X_hat`      the smoothed state at the FIRST OBSERVATION epoch, `rts.t[1]`
#   `to_posvel(sc)`  the back-propagated estimate at the REFERENCE epoch, t = 0
#
# The second is what `run_iterated_rts!` exists to produce: the t = 0 linearization point the
# next sweep starts from, which is also the initial state a user wants out of an OD fit.
"""Propagate the truth to `dt` seconds after the reference epoch."""
function _orts_truth_at(dt)
    sc = _od_sc(_OD_TRUTH)
    dt == 0 && return to_posvel(sc)
    propagate!(_od_prop(), sc, StopAt(sc, PropDurationSeconds(), dt))
    return to_posvel(sc)
end

@testset "iterated RTS — recovers the truth and reports its own convergence" begin
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR
    sc, y0, problem = _orts_setup(guess)

    res = run_iterated_rts!([y0], recs; model = problem, max_iters = 6, tol = 1e-8)

    @test res isa AstroSolve.ExtendedKalmanFilter.IteratedRTSResult
    @test 1 <= res.iters_run <= 6
    @test length(res.delta_history) == res.iters_run

    # The first sweep has nothing to compare against, so it seeds the history with NaN rather
    # than a relative change it cannot compute.
    @test isnan(res.delta_history[1])
    @test all(isfinite, res.delta_history[2:end])

    # Later sweeps move the reference-epoch estimate less than earlier ones. This is what
    # Gauss-Newton on a well-posed problem does, and it is the property that fails first if the
    # back-propagation gain C0 has a sign or transpose error.
    if res.iters_run >= 3
        @test res.delta_history[end] <= res.delta_history[2]
    end

    # The smoothed estimate is inside three times the uncertainty it reports, compared against
    # the truth at the epoch it belongs to. Same reasoning as the batch fixture: the geometry
    # sets what is achievable, so the covariance is the honest yardstick.
    X = res.rts.X_hat
    @test all(abs.(X .- _orts_truth_at(res.rts.t[1])) .<= 3 .* res.rts.sigma)

    # The reference-epoch estimate is the product of the iteration, and it is written onto the
    # spacecraft, which is where a user reads it. Assert it against the t = 0 truth directly:
    # its own covariance lives at the first observation epoch, so a sigma bound would be the
    # wrong yardstick here, and a large improvement factor is the claim that matters.
    X0 = to_posvel(sc)
    @test X0 ≈ _ORTS_EKF._flatten_svs([y0]) atol = 1e-12
    @test norm(X0[1:3] .- _OD_TRUTH[1:3]) < norm(_ORTS_ERR[1:3]) / 10
    @test norm(X0[4:6] .- _OD_TRUTH[4:6]) < norm(_ORTS_ERR[4:6]) / 10
end

@testset "iterated RTS — one sweep is the filter followed by the smoother" begin
    # The docstring says max_iters = 1 gives the same result as run_ekf! plus run_rts. Two
    # setups from the same guess, one driven each way, must agree to round-off.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR

    _, y_a, prob_a = _orts_setup(guess)
    one = run_iterated_rts!([y_a], recs; model = prob_a, max_iters = 1)

    _, y_b, prob_b = _orts_setup(guess)
    filt = run_ekf!([y_b], recs; model = prob_b)
    sm   = run_rts(filt)

    @test one.iters_run == 1
    @test one.converged == false          # one sweep has nothing to converge against
    @test one.ekf.X_hat ≈ filt.X_hat atol = 1e-12
    @test one.ekf.P_hat ≈ filt.P_hat atol = 1e-12
    @test one.rts.X_hat ≈ sm.X_hat   atol = 1e-12
    @test one.rts.P_hat ≈ sm.P_hat   atol = 1e-12
    @test length(one.rts.y_smooth) == length(sm.y_smooth)
end

@testset "iterated RTS — the smoother properties hold through the model path" begin
    # The same two exact properties test_correctness_rts.jl establishes for the low-level
    # smoother, re-checked on a result produced through ODProblem, where the covariance comes
    # off Vary and the dynamics off the propagator rather than from hand-written closures.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    sc, y0, problem = _orts_setup(_OD_TRUTH .+ _ORTS_ERR)
    res = run_iterated_rts!([y0], recs; model = problem, max_iters = 3)

    sm, filt = res.rts, res.ekf

    # No future data at the last epoch, so the smoother returns the filter's answer.
    @test sm.y_smooth[end] ≈ filt.records[end].y_post atol = 1e-12
    @test sm.P_smooth[end] ≈ filt.records[end].P_post atol = 1e-12

    # Smoothing cannot lose information anywhere.
    for k in eachindex(sm.P_smooth)
        @test tr(sm.P_smooth[k]) <= tr(filt.records[k].P_post) + 1e-9
        @test isposdef(Symmetric(sm.P_smooth[k]))
    end

    # And it gains some somewhere, or the smoother did nothing.
    @test any(tr(sm.P_smooth[k]) < tr(filt.records[k].P_post) - 1e-12
              for k in 1:length(sm.P_smooth)-1)
end

@testset "iterated RTS — agrees with the batch fit of the same data" begin
    # Cross-method. Both are MAP estimates of the state at the reference epoch from all of the
    # data under the same prior, so they answer the same question by different arithmetic: one
    # filters and smooths, the other accumulates a normal equation. The batch path shares no
    # code with the smoother, which is what makes this independent of the checks above.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR

    sc_i, y_i, prob_i = _orts_setup(guess)
    run_iterated_rts!([y_i], recs; model = prob_i, max_iters = 6, tol = 1e-10)
    X0_iterated = to_posvel(sc_i)          # reference epoch, which is where batch reports too

    _, y_b, prob_b = _orts_setup(guess)
    bt = _OD_BLS.solve_batch_ls!([y_b], recs; model = prob_b, n_iters = 10)

    # Agreement is asserted relative to what the data supports, not as an absolute distance:
    # over a weakly observable arc both estimators sit in the same broad valley, and the claim
    # worth making is that they sit in the same place in it.
    @test all(abs.(X0_iterated .- bt.X_hat) .<= 0.5 .* bt.sigma)
end

@testset "iterated RTS — a vector covariance is the matrix with that diagonal" begin
    # _build_initial_P and init_ekf both accept a Real variance, a vector of variances, or a
    # covariance matrix. A vector and the Diagonal built from it state the same prior, so the
    # two must produce the same estimate; anything else means one branch mis-indexes.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR

    _, y_v, prob_v = _orts_setup(guess; covariance = _ORTS_COV_VEC)
    vec_res = run_iterated_rts!([y_v], recs; model = prob_v, max_iters = 3)

    _, y_m, prob_m = _orts_setup(guess; covariance = _ORTS_COV_MAT)
    mat_res = run_iterated_rts!([y_m], recs; model = prob_m, max_iters = 3)

    @test vec_res.rts.X_hat ≈ mat_res.rts.X_hat atol = 1e-12
    @test vec_res.rts.P_hat ≈ mat_res.rts.P_hat atol = 1e-12
    @test vec_res.iters_run == mat_res.iters_run

    # A covariance of a shape neither branch handles is refused, and the message enumerates the
    # three that are.
    _, y_bad, prob_bad = _orts_setup(guess; covariance = "wide")
    err = try
        run_iterated_rts!([y_bad], recs; model = prob_bad, max_iters = 1)
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("covariance", msg)
    @test occursin("String", msg)          # names what was actually given
end

@testset "iterated RTS — input validation" begin
    # A caller's bad argument is a DomainError naming the value, not an
    # assertion or a silently clamped default.
    _, recs = _od_truth_records(_ORTS_OFFSETS[1:3])
    _, y0, problem = _orts_setup(_OD_TRUTH .+ _ORTS_ERR)

    @test_throws DomainError run_iterated_rts!([y0], recs; model = problem, max_iters = 0)
    @test_throws DomainError run_iterated_rts!([y0], recs; model = problem, step_size = 0.0)
    @test_throws DomainError run_iterated_rts!([y0], recs; model = problem, step_size = 1.5)
    @test_throws DomainError run_iterated_rts!([y0], recs; model = problem, step_size = -0.5)

    # The messages name the constraint and the value, which is what makes them actionable.
    m1 = try; run_iterated_rts!([y0], recs; model = problem, max_iters = 0)
         catch e; sprint(showerror, e); end
    @test occursin("max_iters", m1) && occursin("at least 1", m1)

    m2 = try; run_iterated_rts!([y0], recs; model = problem, step_size = 1.5)
         catch e; sprint(showerror, e); end
    @test occursin("step_size", m2) && occursin("(0, 1]", m2)

    # step_size = 1.0 is the boundary and is allowed — pure Gauss-Newton.
    @test run_iterated_rts!([y0], recs; model = problem,
                            max_iters = 1, step_size = 1.0).iters_run == 1
end

@testset "iterated RTS — damping changes the path, not the answer" begin
    # step_size < 1 takes a shorter step each sweep, so it needs more of them to reach the same
    # place. Both must land together, or the damping is doing something other than scaling.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR

    sc_f, y_f, prob_f = _orts_setup(guess)
    full = run_iterated_rts!([y_f], recs; model = prob_f, max_iters = 12, tol = 1e-10)

    sc_d, y_d, prob_d = _orts_setup(guess)
    damped = run_iterated_rts!([y_d], recs; model = prob_d,
                               max_iters = 12, tol = 1e-10, step_size = 0.5)

    @test damped.iters_run >= full.iters_run
    @test all(abs.(to_posvel(sc_d) .- to_posvel(sc_f)) .<= 0.5 .* full.rts.sigma)
end

@testset "EKF — the ODProblem entry point matches the closure form" begin
    # run_ekf!(svs, records; model) exists so a caller writes no physics or geometry. It builds
    # the closures and forwards, so driving the low-level form with the closures that
    # build_od_closures produces must give the identical filter.
    _, recs = _od_truth_records(_ORTS_OFFSETS)
    guess = _OD_TRUTH .+ _ORTS_ERR

    _, y_m, prob_m = _orts_setup(guess)
    by_model = run_ekf!([y_m], recs; model = prob_m)

    _, y_c, prob_c = _orts_setup(guess)
    dyn, meas, obs_times, obs_data, R_per_obs = _OD_BLS.build_od_closures(recs, prob_c)
    by_closure = run_ekf!([y_c], obs_times, obs_data, dyn, meas;
                          model = prob_c, R = R_per_obs[1],
                          R_per_obs = R_per_obs, t0 = 0.0)

    @test by_model.X_hat ≈ by_closure.X_hat atol = 1e-12
    @test by_model.P_hat ≈ by_closure.P_hat atol = 1e-12
    @test length(by_model.records) == length(by_closure.records)
    @test by_model.records[end].t ≈ by_closure.records[end].t atol = 1e-9

    # Per-observation R is what makes a mixed range and Doppler stream work: the two kinds carry
    # sigmas six orders of magnitude apart, and a single R would weight them alike.
    @test length(unique(R[1, 1] for R in R_per_obs)) == 2
end

@testset "iterated RTS — verbose logging does not change the answer" begin
    # The verbose branches are @info lines in run_iterated_rts! and run_ekf!. Asserting the log
    # text would pin a format; asserting the result is unchanged is the property that matters.
    _, recs = _od_truth_records(_ORTS_OFFSETS[1:10])
    guess = _OD_TRUTH .+ _ORTS_ERR

    _, y_q, prob_q = _orts_setup(guess)
    quiet = run_iterated_rts!([y_q], recs; model = prob_q, max_iters = 2)

    _, y_v, prob_v = _orts_setup(guess)
    loud = @test_logs (:info,) match_mode = :any run_iterated_rts!(
        [y_v], recs; model = prob_v, max_iters = 2, verbose = true)

    @test loud.rts.X_hat ≈ quiet.rts.X_hat atol = 1e-12
    @test loud.iters_run == quiet.iters_run
end
