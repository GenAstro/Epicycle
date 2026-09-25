# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Variational equations, force-model Jacobians, and OrbitODEProblem.
#
# Truth, by section:
#   propagate_with_stm, _to_times   analytic. A forced harmonic oscillator has a closed-form state
#                                   and state transition matrix.
#   propagate_with_sensitivities_   analytic. ẋ = -k x + g u is solved in closed form, and ∂x/∂k
#                                   is ForwardDiff through that closed form.
#   eval_jacobian                   automatic differentiation of the total right-hand side, summed
#                                   the way `propagate!` sums forces.
#   OrbitODEProblem                 `propagate!` for the state; central finite differences of the
#                                   propagated state for Φ and ∂y/∂μ.

using Test
using LinearAlgebra
using ForwardDiff

using EpicycleBase
using AstroEpochs
using AstroStates
using AstroUniverse
using AstroModels
using AstroProp
using AstroProp: Vern9, Tsit5
using AstroProp: propagate_with_stm, propagate_with_stm_to_times,
                 propagate_with_sensitivities_to_times

# ─────────────────────────────────────────────────────────────────────────────
# Generic variational propagation
# ─────────────────────────────────────────────────────────────────────────────

# ẏ₁ = y₂,  ẏ₂ = -ω² y₁ + u. `model` carries ω, `u` is a constant forcing.
_osc_f(y, u, p, t, model)    = [y[2], -model.ω^2 * y[1] + u]
_osc_dfdy(y, u, p, t, model) = [0.0 1.0; -model.ω^2 0.0]

function _osc_truth(y0, u, ω, dt)
    c, s = cos(ω * dt), sin(ω * dt)
    Φ  = [c s/ω; -ω*s c]
    yp = [u / ω^2, 0.0]                  # particular solution
    return yp .+ Φ * (y0 .- yp), Φ
end

@testset "propagate_with_stm — forced harmonic oscillator" begin
    model = (ω = 2.0,)
    y0, u = [0.4, -0.3], 0.25
    t0, t1 = 0.3, 2.1

    y, Φ = propagate_with_stm(_osc_f, _osc_dfdy, y0, t0, t1, u, nothing, model;
                              solver = Vern9())
    y_true, Φ_true = _osc_truth(y0, u, model.ω, t1 - t0)
    @test y ≈ y_true atol = 1e-10
    @test Φ ≈ Φ_true atol = 1e-10

    # Backward integration gives the transition matrix for a negative interval.
    yb, Φb = propagate_with_stm(_osc_f, _osc_dfdy, y0, t1, t0, u, nothing, model;
                                solver = Vern9())
    yb_true, Φb_true = _osc_truth(y0, u, model.ω, t0 - t1)
    @test yb ≈ yb_true atol = 1e-10
    @test Φb ≈ Φb_true atol = 1e-10

    # A zero-length interval returns the initial state and the identity without integrating.
    y_same, Φ_same = propagate_with_stm(_osc_f, _osc_dfdy, y0, t0, t0, u, nothing, model)
    @test y_same == y0
    @test Φ_same == Matrix(1.0I, 2, 2)
end

@testset "propagate_with_stm_to_times — samples along one solve" begin
    model = (ω = 1.5,)
    y0, u, t0 = [1.0, 0.2], -0.1, 0.5
    ts = [t0, t0 + 0.4, t0 + 1.7, t0 + 3.2]

    ys, Φs = propagate_with_stm_to_times(_osc_f, _osc_dfdy, y0, t0, ts, u, nothing, model;
                                         solver = Vern9())
    @test length(ys) == length(ts) && length(Φs) == length(ts)
    for (k, t) in enumerate(ts)
        y_true, Φ_true = _osc_truth(y0, u, model.ω, t - t0)
        @test ys[k] ≈ y_true atol = 1e-8
        @test Φs[k] ≈ Φ_true atol = 1e-8
    end
    @test Φs[1] == Matrix(1.0I, 2, 2)

    # The solve only runs forward, so a sample before t0 is refused rather than extrapolated.
    @test_throws ArgumentError propagate_with_stm_to_times(
        _osc_f, _osc_dfdy, y0, t0, [t0 - 0.1, t0 + 1.0], u, nothing, model)
end

# ẋ = -k x + g u, with k a parameter and g a model gain.
_dec_f(y, u, p, t, model)    = [-p.k * y[1] + model.g * u]
_dec_dfdy(y, u, p, t, model) = fill(-p.k, 1, 1)
_dec_dfdp(y, u, p, t, model) = fill(-y[1], 1, 1)

# Closed form, generic in k so ForwardDiff gives ∂x/∂k.
_dec_truth(x0, k, g, u, τ) = g * u / k + (x0 - g * u / k) * exp(-k * τ)

@testset "propagate_with_sensitivities_to_times — decay with a parameter" begin
    k, x0, t0, u = 0.7, 1.3, 0.2, 0.5
    model = (g = 2.0,)
    ts = [t0, 0.8, 1.5, 2.7]

    ys, Φys, Φps = propagate_with_sensitivities_to_times(
        _dec_f, _dec_dfdy, _dec_dfdp, [x0], t0, ts, u, (k = k,), model; solver = Vern9())

    for (i, t) in enumerate(ts)
        τ = t - t0
        dxdk_true = ForwardDiff.derivative(kk -> _dec_truth(x0, kk, model.g, u, τ), k)
        @test ys[i][1]     ≈ _dec_truth(x0, k, model.g, u, τ) atol = 1e-9
        @test Φys[i][1, 1] ≈ exp(-k * τ)                       atol = 1e-9
        @test Φps[i][1, 1] ≈ dxdk_true                         atol = 1e-8
    end
    @test Φps[1] == zeros(1, 1)                 # no sensitivity at t0

    # With no parameters the sensitivity block is empty, not absent, and B is never called.
    dec_f_fixed    = (y, u, p, t, model) -> [-k * y[1] + model.g * u]
    dec_dfdy_fixed = (y, u, p, t, model) -> fill(-k, 1, 1)
    never          = (y, u, p, t, model) -> error("dfdp called with no parameters")
    ys0, _, Φps0 = propagate_with_sensitivities_to_times(
        dec_f_fixed, dec_dfdy_fixed, never, [x0], t0, [2.7], u, (;), model; solver = Vern9())
    @test size(Φps0[1]) == (1, 0)
    @test ys0[1] ≈ ys[end] atol = 1e-9

    # The solve only runs forward, so a sample before t0 is refused rather than extrapolated.
    @test_throws ArgumentError propagate_with_sensitivities_to_times(
        _dec_f, _dec_dfdy, _dec_dfdp, [x0], t0, [t0 - 0.1], u, (k = k,), model)
end

# ─────────────────────────────────────────────────────────────────────────────
# Force-model Jacobians
# ─────────────────────────────────────────────────────────────────────────────

# A force with no analytic Jacobian methods, so eval_jacobian takes the AD and FD fallbacks.
struct _OpaqueForce{F} <: AstroProp.OrbitODE
    inner::F
end
AstroProp.accel_eval!(f::_OpaqueForce, t::Time, y::AbstractVector, dy::AbstractVector,
                      sc::Spacecraft, params) = accel_eval!(f.inner, t, y, dy, sc, params)

# The right-hand side as `propagate!` assembles it: kinematics once, accelerations summed.
function _total_rhs(fm::ForceModel, t, y, sc)
    dy = zeros(eltype(y), 6)
    dy[1:3] .= y[4:6]
    for force in fm.forces
        acc = zeros(eltype(y), 6)
        accel_eval!(force, t, y, acc, sc, [])
        dy[4:6] .+= acc[4:6]
    end
    return dy
end

_jac_epoch() = Time("2020-10-20T12:00:00", TDB(), ISOT())
_jac_state() = [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]
_jac_sc()    = Spacecraft(state = CartesianState(_jac_state()), time = _jac_epoch(),
                          mass = 1000.0, srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))

@testset "eval_jacobian — state Jacobian against AD of the total right-hand side" begin
    t, y, sc = _jac_epoch(), _jac_state(), _jac_sc()
    cfg = JacobianConfig(partial_y = true)

    models = [
        "point mass with third bodies" =>
            ForceModel(PointMassGravity(earth, (moon, sun))),
        "point mass and SRP" =>
            ForceModel(PointMassGravity(earth, ()), SolarRadiationPressure(earth)),
        "central term and third bodies as separate forces" =>
            ForceModel(PointMassGravity(earth, ()),
                       PointMassGravity(earth, (moon, sun); include_center = false)),
    ]
    for (label, fm) in models
        A_truth = ForwardDiff.jacobian(yy -> _total_rhs(fm, t, yy, sc), y)
        A = eval_jacobian(cfg, fm, y, sc, t).partial_y
        @testset "$label" begin
            @test A ≈ A_truth rtol = 1e-10
        end
    end
end

@testset "eval_jacobian — analytic methods and fallbacks agree" begin
    t, y, sc = _jac_epoch(), _jac_state(), _jac_sc()
    body = CelestialBody("Earth", 398600.4415, 6378.137, 1/298.257223563, 399)
    mv   = ModelVariable(body, Mu())
    cfg  = JacobianConfig(partial_y = true, partial_p = [mv])

    analytic = ForceModel(PointMassGravity(body, ()))
    opaque   = ForceModel(_OpaqueForce(PointMassGravity(body, ())))

    r_analytic = eval_jacobian(cfg, analytic, y, sc, t)
    r_opaque   = eval_jacobian(cfg, opaque, y, sc, t)

    # The plan records which force types have analytic methods.
    @test PointMassGravity in r_analytic.plan.analytic_state
    @test (PointMassGravity, Mu) in r_analytic.plan.analytic_param
    @test isempty(r_opaque.plan.analytic_state)
    @test isempty(r_opaque.plan.analytic_param)

    @test r_analytic.partial_y ≈ r_opaque.partial_y rtol = 1e-10
    @test r_analytic.partial_p[mv] ≈ r_opaque.partial_p[mv] rtol = 1e-6
    @test r_analytic.partial_p[mv] ≈ fd_differentiate_wrt(analytic, body, Mu(), y, sc, t) rtol = 1e-6
    @test body.mu == 398600.4415              # the finite difference restored the field

    # eval_jacobian! reuses its storage, so a second state must not see the first one's values.
    y2 = [-6500.0, 2100.0, 400.0, -2.0, -6.9, 1.1]
    eval_jacobian!(r_analytic, analytic, y2, sc, t)
    @test r_analytic.partial_y ≈ ForwardDiff.jacobian(yy -> _total_rhs(analytic, t, yy, sc), y2) rtol = 1e-10
end

@testset "fd_differentiate_wrt restores the field when the evaluation throws" begin
    body = CelestialBody("Earth", 398600.4415, 6378.137, 1/298.257223563, 399)
    @test_throws ErrorException fd_differentiate_wrt(() -> error("model failed"), body, Mu())
    @test body.mu == 398600.4415
end

# ─────────────────────────────────────────────────────────────────────────────
# OrbitODEProblem
# ─────────────────────────────────────────────────────────────────────────────

_ode_epoch() = Time("2020-10-20T12:00:00", UTC(), ISOT())
_ode_sc(state = CartesianState(_jac_state())) =
    Spacecraft(state = state, time = _ode_epoch(), mass = 1000.0,
               srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
_ode_integ() = IntegratorConfig(Vern9(); abstol = 1e-12, reltol = 1e-12, dt = 60.0)

@testset "OrbitODEProblem — state matches propagate!" begin
    duration = 3600.0
    for (label, fm) in [
            "point mass with third bodies" => ForceModel(PointMassGravity(earth, (moon, sun))),
            "point mass and SRP" =>
                ForceModel(PointMassGravity(earth, ()), SolarRadiationPressure(earth))]
        prop = OrbitPropagator(fm, _ode_integ())
        sc_ref = _ode_sc()
        propagate!(prop, sc_ref, StopAt(sc_ref, PropDurationSeconds(), duration))
        result = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(); duration_s = duration))
        @testset "$label" begin
            @test result.y_final ≈ to_posvel(sc_ref) rtol = 1e-9
            @test result.Φ === nothing
            @test isempty(result.S_p)
            @test result.sol === nothing
        end
    end
end

@testset "OrbitODEProblem — Keplerian input and dense output" begin
    prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), _ode_integ())
    duration = 1800.0
    cart = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(); duration_s = duration))

    kep = KeplerianState(CartesianState(_jac_state()), earth.mu)
    from_kep = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(kep); duration_s = duration))
    @test from_kep.y_final ≈ cart.y_final rtol = 1e-9

    dense = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(); duration_s = duration, dense = true))
    @test dense.sol !== nothing
    @test dense.sol(duration) ≈ cart.y_final rtol = 1e-12
    @test dense.sol(0.0) ≈ _jac_state() rtol = 1e-12
end

@testset "OrbitODEProblem — STM and ∂y/∂μ against finite differences" begin
    body = CelestialBody("Earth", 398600.4415, 6378.137, 1/298.257223563, 399)
    prop = OrbitPropagator(ForceModel(PointMassGravity(body, (moon, sun))), _ode_integ())
    duration = 2700.0
    mv  = ModelVariable(body, Mu())
    y0  = _jac_state()

    result = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(); duration_s = duration,
                                             stm = STMConfig(Φ = true, S_p = [mv])))
    final(state) = AstroProp.solve(
        OrbitODEProblem(prop, _ode_sc(CartesianState(state)); duration_s = duration)).y_final

    @test result.y_final ≈ final(y0) rtol = 1e-10

    Φ_fd = zeros(6, 6)
    for j in 1:6
        h = j <= 3 ? 1e-3 : 1e-6                                  # km, km/s
        e = zeros(6); e[j] = h
        Φ_fd[:, j] = (final(y0 .+ e) .- final(y0 .- e)) ./ (2h)
    end
    @test result.Φ ≈ Φ_fd rtol = 1e-6

    h  = 1e-5 * body.mu          # small enough that the central difference truncation is below 1e-6
    μ0 = body.mu
    set_field!(body, Mu(), μ0 + h); y_plus  = final(y0)
    set_field!(body, Mu(), μ0 - h); y_minus = final(y0)
    set_field!(body, Mu(), μ0)
    @test result.S_p[mv] ≈ (y_plus .- y_minus) ./ (2h) rtol = 1e-6

    # Parameter sensitivity alone, without the STM.
    only_p = AstroProp.solve(OrbitODEProblem(prop, _ode_sc(); duration_s = duration,
                                             stm = STMConfig(S_p = [mv])))
    @test only_p.Φ === nothing
    @test only_p.S_p[mv] ≈ result.S_p[mv] rtol = 1e-9
end

nothing
