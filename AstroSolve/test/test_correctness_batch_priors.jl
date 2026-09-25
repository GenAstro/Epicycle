# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Variable bookkeeping in both estimators: how a prior is stated, what happens to a variable that
# is not being solved for, and what each refuses.
#
# The OD tests drive `solve_batch_ls!` through `ODProblem` with one six-element state and a matrix
# prior, which is one shape out of several the estimator accepts. This file drives the low-level
# form directly with the analytic spring-mass system, where the variables are scalars and the truth
# is a matrix exponential rather than a propagation. That combination reaches the paths the OD
# fixture cannot: a one-component variable, a prior given as a variance or as a vector, a variable
# held fixed while another is estimated, and a prior left unstated altogether.
#
# The claims are exact identities rather than tolerances. A variance `σ²`, the vector `[σ², σ²]`
# and `Diagonal([σ², σ²])` state the same prior, so the three must give the same estimate to
# round-off. A prior wide enough to carry no information must give the same estimate as no prior at
# all. A fixed variable must come back untouched.
#
# The last testsets do the same for `init_ekf`, whose prior and process-noise branches are the
# same shape one file over, and for the two small guards that sit outside either estimator.
#
# This file depends on the fixture in test_correctness_ekf.jl, which runtests.jl includes first.

using AstroSolve
using EpicycleBase
using LinearAlgebra
using Test

const _BP_BLS = AstroSolve.BatchLeastSquares

# The truth this file fits, and the data it fits. Both come off the shared spring-mass fixture.
const _BP_TRUE = [1.15, -0.42]
const _BP_GUESS = (1.15 * 1.05, -0.42 * 0.95)

"""Observations of the truth, and a fresh oscillator to fit them with."""
function _bp_case(; cov_pos, cov_vel, role_vel = SolveFor(), vel_guess = _BP_GUESS[2])
    osc0, _ = _setup2()
    obs = _obs2(_BP_TRUE, _T2, osc0)

    osc = _Osc2(_BP_GUESS[1], vel_guess, osc0.k, osc0.m, osc0.h)
    svs = [Vary(_pos2, osc; guess = _BP_GUESS[1], covariance = cov_pos, name = "x"),
           Vary(_vel2, osc; guess = vel_guess, covariance = cov_vel,
                role = role_vel, name = "v")]
    return osc, svs, obs
end

@testset "batch — scalar variables are estimated one component at a time" begin
    # Every other estimation test declares either a six-element state or goes through the
    # spring-mass wrapper. Driving solve_batch_ls! directly with two one-component variables is
    # what exercises the single-component branches in the reference-state assembly and in the
    # write-back, which index `rng[1]` rather than a range.
    osc, svs, obs = _bp_case(cov_pos = 1e6, cov_vel = 1e6)
    res = _BP_BLS.solve_batch_ls!(svs, _T2, obs, _vec_bundles()...;
                                  model = osc, R = _R2, n_iters = 8, t0 = 0.0)

    @test length(res.X_hat) == 2
    @test res.X_hat ≈ _BP_TRUE atol = 1e-5

    # The estimate is written back onto the subject, one scalar per variable.
    @test osc.position ≈ res.X_hat[1] atol = 1e-12
    @test osc.velocity ≈ res.X_hat[2] atol = 1e-12

    # And the covariance is a covariance over the two scalars.
    @test size(res.P_hat) == (2, 2)
    @test res.sigma ≈ sqrt.(diag(res.P_hat)) atol = 1e-12
    @test isposdef(Symmetric(res.P_hat))
end

@testset "batch — a variance, a vector and a matrix state the same prior" begin
    # The estimator accepts three shapes for the a priori covariance. For a one-component variable
    # they are the same number written three ways, so the estimates must agree to round-off.
    # Anything else means one branch mis-indexes or inverts differently.
    σ² = 1e-2

    o_r, svs_r, obs = _bp_case(cov_pos = σ², cov_vel = σ²)
    r_scalar = _BP_BLS.solve_batch_ls!(svs_r, _T2, obs, _vec_bundles()...;
                                       model = o_r, R = _R2, n_iters = 8, t0 = 0.0)

    o_v, svs_v, _ = _bp_case(cov_pos = [σ²], cov_vel = [σ²])
    r_vector = _BP_BLS.solve_batch_ls!(svs_v, _T2, obs, _vec_bundles()...;
                                       model = o_v, R = _R2, n_iters = 8, t0 = 0.0)

    o_m, svs_m, _ = _bp_case(cov_pos = fill(σ², 1, 1), cov_vel = fill(σ², 1, 1))
    r_matrix = _BP_BLS.solve_batch_ls!(svs_m, _T2, obs, _vec_bundles()...;
                                       model = o_m, R = _R2, n_iters = 8, t0 = 0.0)

    @test r_scalar.X_hat ≈ r_vector.X_hat atol = 1e-12
    @test r_scalar.X_hat ≈ r_matrix.X_hat atol = 1e-12
    @test r_scalar.P_hat ≈ r_vector.P_hat atol = 1e-12
    @test r_scalar.P_hat ≈ r_matrix.P_hat atol = 1e-12

    # A tight prior really does pull toward the guess, which is what makes the agreement above a
    # statement about the three shapes rather than about a prior that does nothing.
    o_w, svs_w, _ = _bp_case(cov_pos = 1e12, cov_vel = 1e12)
    r_weak = _BP_BLS.solve_batch_ls!(svs_w, _T2, obs, _vec_bundles()...;
                                     model = o_w, R = _R2, n_iters = 8, t0 = 0.0)
    @test norm(r_weak.X_hat .- _BP_TRUE) < norm(r_scalar.X_hat .- _BP_TRUE)
end

@testset "batch — no prior is maximum likelihood, not a default prior" begin
    # A variable declared without a covariance contributes nothing to the normal equation, which
    # is a different thing from a wide prior only in that it is exact. With noise-free data both
    # land on the truth, and a prior so wide it carries no information must match.
    osc, svs, obs = _bp_case(cov_pos = nothing, cov_vel = nothing)
    ml = _BP_BLS.solve_batch_ls!(svs, _T2, obs, _vec_bundles()...;
                                 model = osc, R = _R2, n_iters = 8, t0 = 0.0)

    @test ml.X_hat ≈ _BP_TRUE atol = 1e-6

    o_w, svs_w, _ = _bp_case(cov_pos = 1e14, cov_vel = 1e14)
    wide = _BP_BLS.solve_batch_ls!(svs_w, _T2, obs, _vec_bundles()...;
                                   model = o_w, R = _R2, n_iters = 8, t0 = 0.0)
    @test ml.X_hat ≈ wide.X_hat atol = 1e-6
end

@testset "batch — a fixed variable is held and left out of the estimate" begin
    # Fixed says the variable is part of the state the dynamics need but is not being solved for.
    # Pin velocity at the truth so position can still reach it; pinned at the guess instead, the
    # fit lands 5e-4 away, which is the estimator being right about a model it was handed wrong.
    osc, svs, obs = _bp_case(cov_pos = 1e6, cov_vel = nothing,
                             role_vel = Fixed(), vel_guess = _BP_TRUE[2])
    v_before = osc.velocity

    res = _BP_BLS.solve_batch_ls!(svs, _T2, obs, _vec_bundles()...;
                                  model = osc, R = _R2, n_iters = 8, t0 = 0.0)

    # X_hat is the full state, fixed slots included, while P_hat and sigma cover only what was
    # estimated. solve_for_idx is the bridge between the two, and reading X_hat without it is the
    # mistake this asserts against.
    @test length(res.X_hat) == 2
    @test res.solve_for_idx == [1]
    @test size(res.P_hat) == (1, 1)
    @test length(res.sigma) == 1
    @test size(res.corr) == (1, 1)

    # The fixed variable did not move, in the subject or in the returned state.
    @test osc.velocity == v_before
    @test res.X_hat[2] == v_before

    # And the estimated one found the truth, which it can because the measurement sees it and the
    # velocity it was pinned at is right.
    @test osc.position ≈ _BP_TRUE[1] atol = 1e-4
    @test res.X_hat[res.solve_for_idx[1]] ≈ osc.position atol = 1e-12
end

@testset "batch — a prior of the wrong shape is refused" begin
    # Name the variable, the shape required, and what was given.
    # A prior silently reshaped or ignored would bias every estimate that used it.
    msg_of(f) = try; f(); ""; catch e; sprint(showerror, e); end

    # A vector whose length does not match the variable's component count.
    o, svs, obs = _bp_case(cov_pos = [1e-2, 1e-2], cov_vel = 1e-2)
    e1 = try
        _BP_BLS.solve_batch_ls!(svs, _T2, obs, _vec_bundles()...;
                                model = o, R = _R2, n_iters = 1, t0 = 0.0)
    catch e
        e
    end
    @test e1 isa ArgumentError
    m1 = sprint(showerror, e1)
    @test occursin("\"x\"", m1) || occursin("x", m1)
    @test occursin("length 1", m1) && occursin("got 2", m1)

    # A matrix whose size does not match.
    o2, svs2, _ = _bp_case(cov_pos = Matrix{Float64}(I, 2, 2), cov_vel = 1e-2)
    m2 = msg_of(() -> _BP_BLS.solve_batch_ls!(svs2, _T2, obs, _vec_bundles()...;
                                              model = o2, R = _R2, n_iters = 1, t0 = 0.0))
    @test occursin("(1, 1)", m2) && occursin("(2, 2)", m2)

    # Something that is none of the three, with the three enumerated.
    o3, svs3, _ = _bp_case(cov_pos = "wide", cov_vel = 1e-2)
    m3 = msg_of(() -> _BP_BLS.solve_batch_ls!(svs3, _T2, obs, _vec_bundles()...;
                                              model = o3, R = _R2, n_iters = 1, t0 = 0.0))
    @test occursin("Real variance", m3)
    @test occursin("vector of variances", m3)
    @test occursin("covariance matrix", m3)
    @test occursin("String", m3)
end

@testset "batch — verbose logging does not change the answer" begin
    # The verbose branch is an @info per iteration. Pinning its text would pin a format; the
    # property that matters is that asking for it changes nothing about the fit.
    o_q, svs_q, obs = _bp_case(cov_pos = 1e-2, cov_vel = 1e-2)
    quiet = _BP_BLS.solve_batch_ls!(svs_q, _T2, obs, _vec_bundles()...;
                                    model = o_q, R = _R2, n_iters = 4, t0 = 0.0)

    o_v, svs_v, _ = _bp_case(cov_pos = 1e-2, cov_vel = 1e-2)
    loud = @test_logs (:info,) match_mode = :any _BP_BLS.solve_batch_ls!(
        svs_v, _T2, obs, _vec_bundles()...;
        model = o_v, R = _R2, n_iters = 4, t0 = 0.0, verbose = true)

    @test loud.X_hat ≈ quiet.X_hat atol = 1e-12
    @test loud.P_hat ≈ quiet.P_hat atol = 1e-12
end

@testset "batch — dynamics without a registered Jacobian fall back to AD" begin
    # A DynamicsFunction may register an analytic state Jacobian or leave the estimator to
    # differentiate it. Both paths must fit the same data to the same answer, which is the
    # AD-versus-analytic check done through the estimator rather than
    # on the Jacobian in isolation.
    osc, svs, obs = _bp_case(cov_pos = 1e-2, cov_vel = 1e-2)
    analytic = _BP_BLS.solve_batch_ls!(svs, _T2, obs, _vec_bundles()...;
                                       model = osc, R = _R2, n_iters = 8, t0 = 0.0)

    # The same bundles with no dynamics Jacobian registered. The measurement keeps its own, so
    # only the dynamics side changes.
    dyn_ad  = AstroSolve.DynamicsFunction(_dyn_v; name = :spring_ad)
    meas_ad = AstroSolve.MeasurementFunction(_meas_v; name = :range_rate)
    AstroSolve.add_jacobian!(_mjac_v, meas_ad, AstroSolve.State())

    o2, svs2, _ = _bp_case(cov_pos = 1e-2, cov_vel = 1e-2)
    by_ad = _BP_BLS.solve_batch_ls!(svs2, _T2, obs, dyn_ad, meas_ad;
                                    model = o2, R = _R2, n_iters = 8, t0 = 0.0)

    @test by_ad.X_hat ≈ analytic.X_hat atol = 1e-8
    @test by_ad.P_hat ≈ analytic.P_hat rtol = 1e-6
end

@testset "EKF — init_ekf refuses a prior or a process noise it cannot use" begin
    # The same three covariance shapes and the same enumerated message as the batch estimator,
    # one file over. These are the branches a caller reaches by passing the wrong thing to Vary,
    # which Vary itself does not check because it does not know which estimator will read it.
    msg_of(f) = try; f(); ""; catch e; sprint(showerror, e); end
    dyn, meas = _vec_bundles()

    # A covariance that is none of the three shapes.
    osc = _Osc2(1.2, -0.4, 3.0, 1.5, 5.4)
    bad_cov = [Vary(_pos2, osc; guess = 1.2, covariance = "wide", name = "x")]
    @test_throws ArgumentError _EKF.init_ekf(bad_cov, dyn, meas;
                                             model = osc, R = _R2, t0 = 0.0)
    m = msg_of(() -> _EKF.init_ekf(bad_cov, dyn, meas; model = osc, R = _R2, t0 = 0.0))
    @test occursin("Real variance", m) && occursin("covariance matrix", m) && occursin("String", m)

    # Process noise that is not a ProcessNoiseModel. Vary insists it comes with a covariance but
    # does not check its type, so the filter is where this is caught.
    osc2 = _Osc2(1.2, -0.4, 3.0, 1.5, 5.4)
    bad_pn = [Vary(_pos2, osc2; guess = 1.2, covariance = 1e-2, process_noise = 0.5, name = "x")]
    @test_throws ArgumentError _EKF.init_ekf(bad_pn, dyn, meas;
                                             model = osc2, R = _R2, t0 = 0.0)
    m = msg_of(() -> _EKF.init_ekf(bad_pn, dyn, meas; model = osc2, R = _R2, t0 = 0.0))
    @test occursin("process_noise", m) || occursin("ProcessNoiseModel", m)

    # A real ProcessNoiseModel is accepted and reaches the time update, where it grows the
    # covariance between observations rather than being carried and ignored.
    # One entry per component, and each of these variables is a single scalar.
    osc3, svs3 = _setup2(guess = (1.2, -0.4), cov = (1e-2, 1e-2),
                         pn = AstroSolve.DiagonalSNC([1e-6]))
    ekf = _EKF.init_ekf(svs3, dyn, meas; model = osc3, R = _R2, t0 = 0.0)
    P_before = _EKF.current_covariance(ekf)
    _EKF.time_update!(ekf, 1.0)
    @test tr(_EKF.current_covariance(ekf)) > tr(P_before)
end

@testset "EKF — a time update to the current epoch is a no-op" begin
    # Two observations at the same epoch is ordinary: a station reports range and Doppler
    # together. The second update must not re-propagate, and the state transition it reports for
    # that step must be the identity, or the smoother would compound a step that did not happen.
    osc, svs = _setup2(guess = (1.2, -0.4), cov = (1e-2, 1e-2))
    dyn, meas = _vec_bundles()
    ekf = _EKF.init_ekf(svs, dyn, meas; model = osc, R = _R2, t0 = 0.0)

    _EKF.time_update!(ekf, 1.0)
    y_after = _EKF.current_state(ekf)
    P_after = _EKF.current_covariance(ekf)

    _EKF.time_update!(ekf, 1.0)                      # same epoch again
    @test _EKF.current_time(ekf) ≈ 1.0 atol = 1e-14
    @test _EKF.current_state(ekf) == y_after
    @test _EKF.current_covariance(ekf) ≈ P_after atol = 1e-14
    @test ekf.Phi_last ≈ Matrix{Float64}(I, size(ekf.Phi_last, 1), size(ekf.Phi_last, 2))

    # An epoch in the past is treated the same way rather than propagating backwards.
    _EKF.time_update!(ekf, 0.5)
    @test _EKF.current_time(ekf) ≈ 1.0 atol = 1e-14
end

@testset "EKF — verbose logging does not change the answer" begin
    # The verbose branch reports prefit and postfit RMS per observation, which is what an analyst
    # watches a filter through. Pinning the text would pin a format; the property is that the
    # numbers it reports on are the same either way.
    osc0, _ = _setup2()
    obs = _obs2(_BP_TRUE, _T2, osc0)
    dyn, meas = _vec_bundles()

    o_q, svs_q = _setup2(guess = (1.2, -0.4), cov = (1e-2, 1e-2))
    quiet = _EKF.run_ekf!(svs_q, _T2, obs, dyn, meas; model = o_q, R = _R2, t0 = 0.0)

    o_v, svs_v = _setup2(guess = (1.2, -0.4), cov = (1e-2, 1e-2))
    loud = @test_logs (:info,) match_mode = :any _EKF.run_ekf!(
        svs_v, _T2, obs, dyn, meas; model = o_v, R = _R2, t0 = 0.0, verbose = true)

    @test loud.X_hat ≈ quiet.X_hat atol = 1e-12
    @test loud.P_hat ≈ quiet.P_hat atol = 1e-12

    # mean_sq, which the log line calls, is defined for an empty vector so a filter with no
    # residuals to report logs a zero rather than a NaN.
    @test _EKF.mean_sq(Float64[]) == 0.0
    @test _EKF.mean_sq([3.0, 4.0]) ≈ 12.5
end

@testset "batch — a TDM kind with no measurement type is refused" begin
    # build_od_closures maps a record's kind to the measurement spec type that predicts it. The
    # map is a small closed set, so an unmapped kind names the ones that work rather than sending
    # the reader to the docstring.
    @test_throws ArgumentError _BP_BLS._spec_type_for_tdm(:ANGLE_1)
    m = try; _BP_BLS._spec_type_for_tdm(:ANGLE_1); catch e; sprint(showerror, e); end
    @test occursin("ANGLE_1", m) && occursin("RANGE", m) && occursin("DOPPLER", m)

    # The two that are mapped resolve to the specs that predict them.
    @test _BP_BLS._spec_type_for_tdm(:RANGE)   === TwoWayRange
    @test _BP_BLS._spec_type_for_tdm(:DOPPLER) === TwoWayDoppler
end
