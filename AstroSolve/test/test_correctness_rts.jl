# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Rauch-Tung-Striebel smoother, and the stateful filter API.
#
# A smoother uses every observation to estimate the state at every epoch, where the filter at
# epoch k has only seen data up to k. Two consequences are exact and are what is asserted here:
# the smoothed estimate at the last epoch is the filtered one, because there is no future data
# to add, and the smoothed covariance is never larger than the filtered covariance anywhere,
# because information cannot be negative. On noise-free data the smoothed trajectory must also
# sit on the truth at every epoch, not only at the end.
#
# The smoother's estimate at the first epoch is an estimate of the initial state, which makes it
# comparable to a batch fit of the same data — a second cross-method check, independent of the
# one in the filter tests.
#
# This file depends on the fixture in test_correctness_ekf.jl, which runtests.jl includes first.

using AstroSolve
using EpicycleBase
using LinearAlgebra
using Test

@testset "RTS — smoothing never loses information" begin
    osc0, _ = _setup2()
    X_true  = [1.15, -0.42]
    obs     = _obs2(X_true, _T2, osc0)

    near      = (X_true[1] * 1.02, X_true[2] * 0.98)
    osc, svs  = _setup2(guess = near, cov = (1e-2, 1e-2))
    dyn, meas = _vec_bundles()
    filt = _EKF.run_ekf!(svs, _T2, obs, dyn, meas; model = osc, R = _R2, t0 = 0.0)

    sm = _EKF.run_rts(filt)

    @test length(sm.t) == length(_T2)
    @test length(sm.y_smooth) == length(_T2)
    @test length(sm.P_smooth) == length(_T2)
    @test sm.t == [r.t for r in filt.records]

    # At the last epoch there is no future data, so the smoother returns the filter's answer.
    @test sm.y_smooth[end] ≈ filt.records[end].y_post atol = 1e-14
    @test sm.P_smooth[end] ≈ filt.records[end].P_post atol = 1e-14

    # Everywhere else the smoothed covariance is no larger than the filtered one. This is the
    # property that makes a smoother worth running, and a sign error in the gain breaks it.
    for k in eachindex(sm.P_smooth)
        @test tr(sm.P_smooth[k]) <= tr(filt.records[k].P_post) + 1e-9
        @test isposdef(Symmetric(sm.P_smooth[k]))
        @test sm.P_smooth[k] ≈ sm.P_smooth[k]' atol = 1e-10
    end

    # And strictly smaller somewhere, or the smoother did nothing.
    @test any(tr(sm.P_smooth[k]) < tr(filt.records[k].P_post) - 1e-12
              for k in 1:length(sm.P_smooth)-1)

    # X_hat is the smoothed state at the first epoch, with sigma read off its covariance.
    @test sm.X_hat ≈ sm.y_smooth[1][sm.solve_for_idx] atol = 1e-12
    @test sm.P_hat ≈ sm.P_smooth[1] atol = 1e-14
    @test sm.sigma ≈ sqrt.(diag(sm.P_hat)) atol = 1e-14
    @test all(sm.sigma .> 0)
end

@testset "RTS — the smoothed trajectory sits on the truth" begin
    # Noise-free data, so every smoothed epoch should match the propagated truth. The filter
    # cannot do this at early epochs because it has not yet seen the data that constrains them,
    # which is the difference the smoother exists to make.
    osc0, _ = _setup2()
    X_true  = [1.15, -0.42]
    obs     = _obs2(X_true, _T2, osc0)

    near      = (X_true[1] * 1.05, X_true[2] * 0.95)
    osc, svs  = _setup2(guess = near, cov = (1e-1, 1e-1))
    dyn, meas = _vec_bundles()
    filt = _EKF.run_ekf!(svs, _T2, obs, dyn, meas; model = osc, R = _R2, t0 = 0.0)
    sm   = _EKF.run_rts(filt)

    truth_at = [_truth2(X_true, t, osc0) for t in _T2]

    # The smoothed estimate is on the truth everywhere.
    for k in eachindex(_T2)
        @test sm.y_smooth[k] ≈ truth_at[k] atol = 1e-4
    end

    # And it beats the filter at the first epoch, which is where the filter has least data.
    filt_err = norm(filt.records[1].y_post .- truth_at[1])
    sm_err   = norm(sm.y_smooth[1] .- truth_at[1])
    @test sm_err < filt_err
end

@testset "RTS — smoothed initial state agrees with the batch fit" begin
    # Cross-method, and independent of the filter-versus-batch check: the smoother's estimate at
    # the first epoch and a batch least-squares fit of the same data are both estimates of the
    # initial state from all of the data, so they must agree.
    osc0, _ = _setup2()
    X_true  = [1.15, -0.42]
    obs     = _obs2(X_true, _T2, osc0)
    near    = (X_true[1] * 1.02, X_true[2] * 0.98)

    osc_f, svs_f = _setup2(guess = near, cov = (1e6, 1e6))
    dynv, measv  = _vec_bundles()
    sm = _EKF.run_rts(_EKF.run_ekf!(svs_f, _T2, obs, dynv, measv;
                                    model = osc_f, R = _R2, t0 = 0.0))

    osc_b, svs_b = _setup2(guess = near, cov = (1e6, 1e6))
    dynn, measn  = _nt_bundles()
    batch = _SMB.solve_spring_mass_batch!(svs_b, _T2, obs, dynn, measn;
                                          model = osc_b, n_iters = 8)

    # The smoother's first epoch is t = 0.5, not 0, so propagate the batch answer there.
    batch_at_first = _truth2(batch.X_hat, _T2[1], osc0)
    @test sm.y_smooth[1] ≈ batch_at_first atol = 1e-3
end

@testset "EKF — the stateful API steps by hand" begin
    # init_ekf / time_update! / measurement_update! is the same filter driven one step at a
    # time, which is what a user needs when observations arrive live rather than as a list.
    # Driving it by hand must reproduce what run_ekf! does in one call.
    osc0, _ = _setup2()
    X_true  = [1.1, -0.4]
    obs     = _obs2(X_true, _T2, osc0)
    near    = (X_true[1] * 1.02, X_true[2] * 0.98)

    osc_a, svs_a = _setup2(guess = near, cov = (1e-2, 1e-2))
    dyn_a, meas_a = _vec_bundles()
    batched = _EKF.run_ekf!(svs_a, _T2, obs, dyn_a, meas_a; model = osc_a, R = _R2, t0 = 0.0)

    osc_b, svs_b = _setup2(guess = near, cov = (1e-2, 1e-2))
    dyn_b, meas_b = _vec_bundles()
    ekf = _EKF.init_ekf(svs_b, dyn_b, meas_b; model = osc_b, R = _R2, t0 = 0.0)

    @test _EKF.current_time(ekf) ≈ 0.0 atol = 1e-14
    @test length(_EKF.current_state(ekf)) == 2
    @test size(_EKF.current_covariance(ekf)) == (2, 2)

    for k in eachindex(_T2)
        _EKF.time_update!(ekf, _T2[k])
        @test _EKF.current_time(ekf) ≈ _T2[k] atol = 1e-12
        _EKF.measurement_update!(ekf, obs[k])
    end

    # Same inputs, same arithmetic, so the two paths must agree to round-off.
    @test _EKF.current_state(ekf) ≈ batched.X_hat atol = 1e-10
    @test _EKF.current_covariance(ekf) ≈ batched.P_hat atol = 1e-10

    # Committing writes the estimate onto the subjects, which the batched driver does for you.
    # The filter already holds its variables, so the call takes no second argument.
    @test osc_b.position != _EKF.current_state(ekf)[1]      # not yet written
    _EKF.commit_to_svs!(ekf)
    @test osc_b.position ≈ _EKF.current_state(ekf)[1] atol = 1e-12
    @test osc_b.velocity ≈ _EKF.current_state(ekf)[2] atol = 1e-12
end

@testset "RTS — input validation" begin
    # An empty record stream is an ArgumentError that says what
    # to do about it rather than only that something was empty.
    empty_result = _EKF.EKFResult(Float64[], zeros(0, 0), Float64[],
                                  _EKF.EKFRecord[], Int[])
    @test_throws ArgumentError _EKF.run_rts(empty_result)
    msg = try; _EKF.run_rts(empty_result); catch e; sprint(showerror, e); end
    @test occursin("run_ekf!", msg)          # names the call that produces records
end
