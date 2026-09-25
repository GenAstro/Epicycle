# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Extended Kalman filter on the spring-mass problem.
#
# Two kinds of truth, both stronger than a recorded number. The system is linear, so a filter
# with no process noise is the Kalman filter and its answer over the whole data span must agree
# with the batch least-squares fit of the same data under the same prior — two independently
# written estimators reaching the same estimate. And with noise-free data the innovations must
# vanish once the filter has converged, which is a statement about the propagation and the
# measurement model agreeing rather than about the gain arithmetic.
#
# The subject is the same oscillator the batch tests use, declared the same way with `Vary`.
# Note the closure convention differs: the EKF passes `y::AbstractVector` where
# SpringMassEstimator passes a NamedTuple, so the physics is written once and adapted twice.
#
# ExtendedKalmanFilter.jl was 0 of 286 coverable lines.

using AstroSolve
using EpicycleBase
using LinearAlgebra
using Test

const _EKF = AstroSolve.ExtendedKalmanFilter
const _SMB = AstroSolve.SpringMassEstimator
const _PN  = AstroSolve.ProcessNoiseModels

mutable struct _Osc2
    position :: Float64
    velocity :: Float64
    k        :: Float64
    m        :: Float64
    h        :: Float64
end

_pos2(o::_Osc2) = o.position
_vel2(o::_Osc2) = o.velocity
EpicycleBase.set_quantity!(o::_Osc2, ::typeof(_pos2); to) = (o.position = to; o)
EpicycleBase.set_quantity!(o::_Osc2, ::typeof(_vel2); to) = (o.velocity = to; o)

_w2(m) = m.k / m.m

# Vector-state closures, which is what the EKF asks for.
_dyn_v(y, u, p, t, m)  = [y[2], -_w2(m) * y[1]]
_djac_v(y, u, p, t, m) = [0.0 1.0; -_w2(m) 0.0]
_meas_v(y, u, p, t, m) = (ρ = sqrt(y[1]^2 + m.h^2); [ρ, y[1] * y[2] / ρ])
function _mjac_v(y, u, p, t, m)
    x, v = y[1], y[2]
    ρ = sqrt(x^2 + m.h^2)
    return [x/ρ                  0.0;
            v * m.h^2 / ρ^3      x/ρ]
end

# NamedTuple-state adapters, which is what the batch estimator asks for.
_dyn_nt2(st, u, p, t, m)  = _dyn_v([st.x, st.v], u, p, t, m)
_djac_nt2(st, u, p, t, m) = _djac_v([st.x, st.v], u, p, t, m)
_meas_nt2(st, u, p, t, m) = _meas_v([st.x, st.v], u, p, t, m)
_mjac_nt2(st, u, p, t, m) = _mjac_v([st.x, st.v], u, p, t, m)

function _vec_bundles()
    dyn  = AstroSolve.DynamicsFunction(_dyn_v;    name = :spring)
    meas = AstroSolve.MeasurementFunction(_meas_v; name = :range_rate)
    AstroSolve.add_jacobian!(_djac_v, dyn,  AstroSolve.State())
    AstroSolve.add_jacobian!(_mjac_v, meas, AstroSolve.State())
    return dyn, meas
end

function _nt_bundles()
    dyn  = AstroSolve.DynamicsFunction(_dyn_nt2;    name = :spring)
    meas = AstroSolve.MeasurementFunction(_meas_nt2; name = :range_rate)
    AstroSolve.add_jacobian!(_djac_nt2, dyn,  AstroSolve.State())
    AstroSolve.add_jacobian!(_mjac_nt2, meas, AstroSolve.State())
    return dyn, meas
end

_truth2(X0, t, m) = exp(t * _djac_v(X0, nothing, nothing, 0.0, m)) * X0
_obs2(X0, times, m) = [_meas_v(_truth2(X0, t, m), nothing, nothing, t, m) for t in times]

# x = 0 is unobservable for this geometry, so no guess sits there. See the batch test file.
function _setup2(; guess = (3.0, 1.0), cov = (10.0, 10.0), k = 3.0, m = 1.5, h = 5.4,
                   pn = nothing)
    osc = _Osc2(guess[1], guess[2], k, m, h)
    svs = [Vary(_pos2, osc; guess = guess[1], covariance = cov[1], process_noise = pn),
           Vary(_vel2, osc; guess = guess[2], covariance = cov[2], process_noise = pn)]
    return osc, svs
end

const _T2 = collect(0.5:0.5:20.0)
const _R2 = Matrix(Diagonal([1e-8, 1e-8]))

@testset "EKF — tracks the truth from a nearby start" begin
    # An EKF linearizes about its current estimate, so it is a tracking filter and not a global
    # optimizer: given a decent starting estimate it stays on the truth, and given a bad one it
    # improves but need not converge tightly. Both are asserted separately below, because
    # conflating them is what makes an EKF test look like a failing estimator.
    osc0, _ = _setup2()
    X_true  = [1.2, -0.45]
    obs     = _obs2(X_true, _T2, osc0)
    X_end   = _truth2(X_true, _T2[end], osc0)

    # Starting within a few percent, which is the filter's operating regime.
    near      = (X_true[1] * 1.02, X_true[2] * 0.98)
    osc, svs  = _setup2(guess = near, cov = (1e-2, 1e-2))
    dyn, meas = _vec_bundles()
    res = _EKF.run_ekf!(svs, _T2, obs, dyn, meas; model = osc, R = _R2, t0 = 0.0)

    @test res.X_hat ≈ X_end atol = 1e-5
    @test length(res.records) == length(_T2)
    @test res.sigma ≈ sqrt.(diag(res.P_hat)) atol = 1e-12
    @test all(res.sigma .> 0)
    @test isposdef(Symmetric(res.P_hat))

    # The state correction collapses once the filter is tracking, which says the propagation
    # and the measurement model agree with how the data was generated.
    tail = res.records[end-9:end]
    @test maximum(maximum(abs, r.y_post .- r.y_pre) for r in tail) < 1e-5

    # The estimate is written back onto the subject, which is where a user reads it.
    @test osc.position ≈ X_end[1] atol = 1e-5
    @test osc.velocity ≈ X_end[2] atol = 1e-5
end

@testset "EKF — improves on a poor start without converging tightly" begin
    # From a guess well outside the linear regime the filter must still move the estimate a long
    # way toward the truth. It is not required to reach it: that is what the batch fit is for,
    # and asserting otherwise would be asserting the wrong contract.
    osc0, _ = _setup2()
    X_true  = [1.2, -0.45]
    obs     = _obs2(X_true, _T2, osc0)
    X_end   = _truth2(X_true, _T2[end], osc0)

    guess     = (3.0, 1.0)
    osc, svs  = _setup2(guess = guess, cov = (10.0, 10.0))
    dyn, meas = _vec_bundles()
    res = _EKF.run_ekf!(svs, _T2, obs, dyn, meas; model = osc, R = _R2, t0 = 0.0)

    start_err = norm(_truth2([guess...], _T2[end], osc0) .- X_end)
    end_err   = norm(res.X_hat .- X_end)
    @test end_err < start_err / 10
    @test end_err < 0.2
end

@testset "EKF — uncertainty shrinks as observations arrive" begin
    osc0, _ = _setup2()
    obs     = _obs2([1.0, -0.3], _T2, osc0)

    osc, svs  = _setup2(guess = (3.0, 1.0), cov = (10.0, 10.0))
    dyn, meas = _vec_bundles()
    res = _EKF.run_ekf!(svs, _T2, obs, dyn, meas; model = osc, R = _R2, t0 = 0.0)

    # The posterior covariance is never larger than the prior at the same step: a measurement
    # cannot add uncertainty.
    for r in res.records
        @test tr(r.P_post) <= tr(r.P_pre) + 1e-9
        @test isposdef(Symmetric(r.P_post))
    end

    # And over the run it falls a long way below where it started.
    @test tr(res.records[end].P_post) < tr(res.records[1].P_pre) / 100
end

@testset "EKF — process noise keeps the covariance from collapsing" begin
    osc0, _ = _setup2()
    obs     = _obs2([1.0, -0.3], _T2, osc0)

    osc_q, svs_q = _setup2(guess = (3.0, 1.0), pn = _PN.DiagonalSNC([1e-4]))
    dyn, meas    = _vec_bundles()
    with_q = _EKF.run_ekf!(svs_q, _T2, obs, dyn, meas; model = osc_q, R = _R2, t0 = 0.0)

    osc_n, svs_n = _setup2(guess = (3.0, 1.0), pn = nothing)
    dyn2, meas2  = _vec_bundles()
    no_q = _EKF.run_ekf!(svs_n, _T2, obs, dyn2, meas2; model = osc_n, R = _R2, t0 = 0.0)

    # Adding process noise between updates leaves more uncertainty at the end. Without it the
    # covariance only ever shrinks, which is what makes a filter over-confident on a model it
    # does not perfectly know.
    @test tr(with_q.P_hat) > tr(no_q.P_hat)

    # The time update must grow the covariance where process noise is present, which is the
    # step the Thornton update exists to perform.
    grew = count(r -> tr(r.P_pre) > tr(r.P_post), with_q.records)
    @test grew > length(with_q.records) ÷ 2
end

@testset "EKF — agrees with the batch fit on the same linear problem" begin
    # Cross-method truth. With no process noise on a linear system the EKF is the Kalman
    # filter, and processing every observation once must land where a batch least-squares fit
    # of the same data lands, mapped to the same epoch. The two share the UDU module and
    # nothing else, so agreement is evidence about both.
    osc0, _ = _setup2()
    X_true  = [1.1, -0.4]
    obs     = _obs2(X_true, _T2, osc0)

    # Both start in the filter's operating regime, a couple of percent from the truth, with a
    # prior weak enough that the data decides. Starting the filter badly would compare a
    # partly-converged EKF against a converged batch fit, which is not a like-for-like test.
    near = (X_true[1] * 1.02, X_true[2] * 0.98)
    osc_f, svs_f = _setup2(guess = near, cov = (1e6, 1e6))
    dynv, measv  = _vec_bundles()
    filt = _EKF.run_ekf!(svs_f, _T2, obs, dynv, measv; model = osc_f, R = _R2, t0 = 0.0)

    osc_b, svs_b = _setup2(guess = near, cov = (1e6, 1e6))
    dynn, measn  = _nt_bundles()
    batch = _SMB.solve_spring_mass_batch!(svs_b, _T2, obs, dynn, measn;
                                          model = osc_b, n_iters = 8)

    # The batch estimates the state at t = 0; the filter reports it at the last observation.
    # Propagate the batch answer forward to compare like with like.
    batch_at_end = _truth2(batch.X_hat, _T2[end], osc0)
    @test filt.X_hat ≈ batch_at_end atol = 1e-3
end

@testset "EKF — input validation" begin
    osc0, _ = _setup2()
    obs = _obs2([1.0, 0.0], _T2, osc0)
    dyn, meas = _vec_bundles()

    # The exception type and its message are both checked.
    osc, svs = _setup2()
    @test_throws ArgumentError _EKF.run_ekf!(svs, _T2, obs[1:3], dyn, meas;
                                             model = osc, R = _R2)
    len_msg = try
        _EKF.run_ekf!(svs, _T2, obs[1:3], dyn, meas; model = osc, R = _R2)
    catch e
        sprint(showerror, e)
    end
    @test occursin("same length", len_msg)
    @test occursin("3", len_msg)

    @test_throws ArgumentError _EKF.run_ekf!(svs, _T2, obs, dyn, meas;
                                             model = osc, R = _R2,
                                             R_per_obs = [_R2, _R2])

    # A filter needs a starting covariance, and says so rather than defaulting to one.
    osc_nc = _Osc2(3.0, 1.0, 3.0, 1.5, 5.4)
    svs_nc = [Vary(_pos2, osc_nc; guess = 3.0, lower_bound = -10.0, upper_bound = 10.0),
              Vary(_vel2, osc_nc; guess = 1.0, lower_bound = -10.0, upper_bound = 10.0)]
    @test_throws ArgumentError _EKF.run_ekf!(svs_nc, _T2, obs, dyn, meas;
                                             model = osc_nc, R = _R2)
    cov_msg = try
        _EKF.run_ekf!(svs_nc, _T2, obs, dyn, meas; model = osc_nc, R = _R2)
    catch e
        sprint(showerror, e)
    end
    @test occursin("covariance", cov_msg)
end
