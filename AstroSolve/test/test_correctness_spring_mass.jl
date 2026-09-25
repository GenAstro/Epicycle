# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Batch least squares on the spring-mass problem — analytic truth.
#
# The system is linear and time-invariant, so noise-free observations generated from a known
# initial state must be recovered to numerical precision, not approximately. Every assertion
# here is therefore exact and no tolerance is chosen by judgement. The problem is Tapley,
# Schutz & Born, *Statistical Orbit Determination*, §4.8.2.
#
# The variables are declared with `Vary`, which is the verb a user writes. That is deliberate:
# building a `SolverVariable` by hand needs a `calc` satisfying get_calc/set_calc!, which is
# not the documented path and not what a script looks like. Testing through `Vary` means this
# file breaks if the user-facing surface breaks.
#
# SpringMassEstimator.jl was 0 of 51 coverable lines.

using AstroSolve
using EpicycleBase
using LinearAlgebra
using Test

const _SM = AstroSolve.SpringMassEstimator

# ── The subject ──────────────────────────────────────────────────────────────
# A mass on a spring: ẍ = −(k/m)·x, observed in range and range-rate from a station offset
# perpendicular to the motion by h, so the measurement is smooth everywhere.

mutable struct _Oscillator
    position :: Float64
    velocity :: Float64
    k        :: Float64
    m        :: Float64
    h        :: Float64
end

_position(o::_Oscillator) = o.position
_velocity(o::_Oscillator) = o.velocity

# Settable on an oscillator specifically. Declaring these settable in the abstract would claim
# it of every subject in the system.
EpicycleBase.set_quantity!(o::_Oscillator, ::typeof(_position); to) = (o.position = to; o)
EpicycleBase.set_quantity!(o::_Oscillator, ::typeof(_velocity); to) = (o.velocity = to; o)

# ── The physics ──────────────────────────────────────────────────────────────
_ω²(m) = m.k / m.m
_spring_dyn(y, u, p, t, m)     = [y[2], -_ω²(m) * y[1]]
_spring_jac(y, u, p, t, m)     = [0.0 1.0; -_ω²(m) 0.0]

# ρ = √(x² + h²), ρ̇ = x·ẋ / ρ
_range_rate(y, u, p, t, m) = (ρ = sqrt(y[1]^2 + m.h^2); [ρ, y[1] * y[2] / ρ])
function _range_rate_jac(y, u, p, t, m)
    x, v = y[1], y[2]
    ρ = sqrt(x^2 + m.h^2)
    return [x/ρ                       0.0;
            v * m.h^2 / ρ^3           x/ρ]
end

# The estimator calls closures as f(state, u, p, t, model) where `state` is a NamedTuple with
# the fields this problem solves for. Adapt to the vector form the physics is written in.
_as_vec(st) = [st.x, st.v]
_dyn_nt(st, u, p, t, m)  = _spring_dyn(_as_vec(st), u, p, t, m)
_jac_nt(st, u, p, t, m)  = _spring_jac(_as_vec(st), u, p, t, m)
_meas_nt(st, u, p, t, m) = _range_rate(_as_vec(st), u, p, t, m)
_mjac_nt(st, u, p, t, m) = _range_rate_jac(_as_vec(st), u, p, t, m)

function _bundles()
    dyn  = AstroSolve.DynamicsFunction(_dyn_nt;    name = :spring)
    meas = AstroSolve.MeasurementFunction(_meas_nt; name = :range_rate)
    AstroSolve.add_jacobian!(_jac_nt,  dyn,  AstroSolve.State())
    AstroSolve.add_jacobian!(_mjac_nt, meas, AstroSolve.State())
    return dyn, meas
end

"""Exact state at `t` from the closed form of the linear system."""
_truth(X0, t, m) = exp(t * _spring_jac(X0, nothing, nothing, 0.0, m)) * X0

"""Noise-free observations of the truth trajectory."""
_obs(X0, times, m) = [_range_rate(_truth(X0, t, m), nothing, nothing, t, m) for t in times]

"""A fresh oscillator plus the two variables the estimator solves for."""
function _setup(; guess = (4.0, 0.2), cov = (1e14, 1e14), k = 3.0, m = 1.5, h = 5.4)
    osc = _Oscillator(guess[1], guess[2], k, m, h)
    svs = [Vary(_position, osc; guess = guess[1], covariance = cov[1]),
           Vary(_velocity, osc; guess = guess[2], covariance = cov[2])]
    return osc, svs
end

const _TIMES = collect(0.0:0.5:20.0)

# Two properties of the problem that set the tolerances below, neither of them a defect.
#
# A batch fit with a finite a priori returns the MAP estimate, not the maximum-likelihood one,
# so it is biased toward the guess in proportion to the a priori weight. Recovering the truth
# exactly needs a prior weak enough that the bias falls under the tolerance, which is why the
# default covariance here is 1e14 rather than a nominal 1e6. At 1e6 the converged answer sits
# about 7e-7 from the truth, which is correct behaviour and would look like a failure.
#
# A starting guess must not sit at x = 0. The measurement is ρ = √(x² + h²), whose derivative
# with respect to x vanishes there, so H̃ is identically zero and the problem is unobservable at
# that single point. This is a property of the observation geometry rather than of the
# estimator, and it is the reason every guess below is offset from the origin.
const _GUESS = (3.0, 1.0)

@testset "SpringMass — recovers the truth from noise-free data" begin
    osc0, _ = _setup()
    X_true  = [1.2, -0.45]
    obs     = _obs(X_true, _TIMES, osc0)

    # A deliberately wrong starting guess, with loose a priori so the data decides.
    osc, svs  = _setup(guess = _GUESS)
    dyn, meas = _bundles()
    res = _SM.solve_spring_mass_batch!(svs, _TIMES, obs, dyn, meas; model = osc, n_iters = 8)

    # Exact data on a linear system: the estimate is the truth to numerical precision.
    @test res.X_hat ≈ X_true atol = 1e-9
    @test res.iters == 8

    # Noise-free data leaves no residual, which is what shows the reference trajectory and the
    # measurement model agree with how the data was generated.
    @test maximum(maximum(abs, r) for r in res.residuals) < 1e-8
    @test length(res.residuals) == length(_TIMES)

    # The formal covariance is symmetric positive definite and sigma is read off its diagonal.
    @test res.P_hat ≈ res.P_hat' atol = 1e-12
    @test isposdef(Symmetric(res.P_hat))
    @test res.sigma ≈ sqrt.(diag(res.P_hat)) atol = 1e-14
    @test all(res.sigma .> 0)
    @test abs(res.corr) <= 1.0

    # The estimate is written back onto the subject through the variable contract, which is
    # where a user reads it rather than from the result struct.
    @test osc.position ≈ X_true[1] atol = 1e-9
    @test osc.velocity ≈ X_true[2] atol = 1e-9
    @test AstroSolve.current_value(svs[1]) ≈ X_true[1] atol = 1e-9
end

@testset "SpringMass — the solution is unique up to the observation's symmetry" begin
    # This geometry cannot distinguish (x, v) from (−x, −v). ρ = √(x² + h²) is even in x and
    # ρ̇ = x·ẋ/ρ is odd in x and v together, so the two states produce identical observations.
    # The fit is therefore only locally identifiable: it converges to whichever branch the
    # starting guess is nearer, and both branches fit the data exactly.
    #
    # This is a property of range-and-range-rate from a station offset perpendicular to the
    # motion, not a defect. Breaking the symmetry needs a station offset along the line of
    # motion, where ρ = |x − x_s| is not even.
    osc0, _ = _setup()
    X_true  = [-0.8, 0.3]
    obs     = _obs(X_true, _TIMES, osc0)

    answers = map(((0.5, 0.1), (5.0, 5.0), (-3.0, 2.0), (-0.2, -0.1))) do guess
        osc, svs  = _setup(guess = guess)
        dyn, meas = _bundles()
        _SM.solve_spring_mass_batch!(svs, _TIMES, obs, dyn, meas;
                                     model = osc, n_iters = 8).X_hat
    end

    # Every answer is the truth or its reflection, to numerical precision.
    for a in answers
        @test isapprox(a, X_true; atol = 1e-8) || isapprox(a, -X_true; atol = 1e-8)
    end

    # Guesses on the same side of the origin land on the same branch, so the outcome is
    # determined by the guess rather than arbitrary.
    @test answers[1] ≈ answers[2] atol = 1e-9        # both started positive
    @test answers[3] ≈ answers[4] atol = 1e-9        # both started negative
    @test answers[1] ≈ -answers[3] atol = 1e-8       # and the two branches are reflections

    # Both branches fit the data exactly, which is what makes them indistinguishable rather
    # than one being a worse fit the estimator settled for.
    for guess in ((5.0, 5.0), (-3.0, 2.0))
        osc, svs  = _setup(guess = guess)
        dyn, meas = _bundles()
        res = _SM.solve_spring_mass_batch!(svs, _TIMES, obs, dyn, meas;
                                           model = osc, n_iters = 8)
        @test maximum(maximum(abs, r) for r in res.residuals) < 1e-8
    end
end

@testset "SpringMass — the a priori decides when data cannot" begin
    # One observation cannot determine two states, so the a priori carries the answer. Tight
    # about the guess must hold the estimate there; loose must not.
    osc0, _ = _setup()
    X_true  = [1.0, -0.5]
    times   = [2.0]
    obs     = _obs(X_true, times, osc0)

    g = [_GUESS[1], _GUESS[2]]

    osc_t, svs_t = _setup(guess = _GUESS, cov = (1e-10, 1e-10))
    dyn, meas    = _bundles()
    tight = _SM.solve_spring_mass_batch!(svs_t, times, obs, dyn, meas;
                                         model = osc_t, n_iters = 4)
    @test maximum(abs, tight.X_hat .- g) < 1e-3          # held at the guess

    osc_l, svs_l = _setup(guess = _GUESS, cov = (1e6, 1e6))
    dyn2, meas2  = _bundles()
    loose = _SM.solve_spring_mass_batch!(svs_l, times, obs, dyn2, meas2;
                                         model = osc_l, n_iters = 4)
    @test maximum(abs, loose.X_hat .- g) > maximum(abs, tight.X_hat .- g)

    # A tighter a priori also reports a smaller formal uncertainty, which is the point of
    # carrying one at all.
    @test all(tight.sigma .< loose.sigma)
end

@testset "SpringMass — more data tightens the covariance" begin
    osc0, _ = _setup()
    X_true  = [0.7, 0.9]

    few_t = collect(0.0:2.0:8.0)
    osc_f, svs_f = _setup(guess = _GUESS)
    dyn_f, meas_f = _bundles()
    few = _SM.solve_spring_mass_batch!(svs_f, few_t, _obs(X_true, few_t, osc0),
                                       dyn_f, meas_f; model = osc_f, n_iters = 8)

    many_t = collect(0.0:0.1:20.0)
    osc_m, svs_m = _setup(guess = _GUESS)
    dyn_m, meas_m = _bundles()
    many = _SM.solve_spring_mass_batch!(svs_m, many_t, _obs(X_true, many_t, osc0),
                                        dyn_m, meas_m; model = osc_m, n_iters = 8)

    @test all(many.sigma .< few.sigma)
    @test few.X_hat  ≈ X_true atol = 1e-8
    @test many.X_hat ≈ X_true atol = 1e-8
end

@testset "SpringMass — input validation" begin
    osc0, _ = _setup()
    obs = _obs([1.0, 0.0], _TIMES, osc0)
    dyn, meas = _bundles()

    # The exception type and its message are both checked: the expected count is stated and
    # the received count named.
    osc, svs = _setup()
    @test_throws ArgumentError _SM.solve_spring_mass_batch!(svs[1:1], _TIMES, obs, dyn, meas;
                                                            model = osc)
    msg = try
        _SM.solve_spring_mass_batch!(svs[1:1], _TIMES, obs, dyn, meas; model = osc)
    catch e
        sprint(showerror, e)
    end
    @test occursin("length 2", msg)
    @test occursin("got 1", msg)

    osc3, svs3 = _setup()
    three = vcat(svs3, [Vary(_position, osc3; guess = 1.0, covariance = 1.0)])
    @test_throws ArgumentError _SM.solve_spring_mass_batch!(three, _TIMES, obs, dyn, meas;
                                                            model = osc3)
end

@testset "SpringMass — Vary rejects contradictory declarations" begin
    # Vary is the user-facing verb, so its refusals are part of the interface being tested. A
    # covariance and a box bound are two different claims about one variable and no solver reads
    # both; process noise without a covariance has nothing to grow.
    osc = _Oscillator(4.0, 0.2, 3.0, 1.5, 5.4)
    @test_throws ArgumentError Vary(_position, osc; guess = 1.0, covariance = 1.0,
                                    lower_bound = 0.0)
    @test_throws ArgumentError Vary(_position, osc; guess = 1.0, process_noise = 1e-9)

    box_msg = try
        Vary(_position, osc; guess = 1.0, covariance = 1.0, lower_bound = 0.0)
    catch e
        sprint(showerror, e)
    end
    @test occursin("covariance", box_msg)

    # And a quantity with no setter for this subject cannot be varied, which is checked at
    # declaration rather than discovered when the solver tries to write.
    @test_throws ArgumentError Vary(_ω², osc; guess = 1.0, covariance = 1.0)
end
