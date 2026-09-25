# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The automatic-differentiation generators a transcription calls when a user declared no partial.
#
# Hermite-Simpson reaches only some of them. The rest are the extension interface another
# transcription builds on, LGL in EpicycleEnterprise among them, so they are checked here directly
# rather than through a solve.
#
# Truth: analytic. The dynamics and the path function below are written so that every partial,
# with respect to the state, the control, the parameters, the time and one model field, has a
# closed form.

using Test
using ForwardDiff
using AstroSolve
using AstroSolve: EvalContext, _make_dynamics_jac_y_ad, _make_dynamics_jac_u_ad,
                  _make_dynamics_jac_p_ad, _make_dynamics_jac_t_ad,
                  _make_dynamics_jac_modelfield_ad, _make_path_jac_y_ad, _make_path_jac_u_ad,
                  _make_path_jac_p_ad, _make_path_jac_t_ad, _make_path_jac_modelfield_ad,
                  _ad_field_deriv

# One type parameter per field, which model-field seeding requires: it rebuilds the model with
# the seeded field a Dual and the others unchanged, so the fields cannot share one parameter.
struct _AdModel{K<:Real, C<:Real}
    k::K
    c::C
end

# ẏ₁ = k y₁ u₁ + p₁ t²,   ẏ₂ = sin y₂ + u₂ t
function _ad_dynamics!(dy, y, ctx, t)
    k, u, p = ctx.model.k, ctx.u, ctx.params
    dy[1] = k * y[1] * u[1] + (isempty(p) ? 0.0 : p[1]) * t^2
    dy[2] = sin(y[2]) + u[2] * t
    return nothing
end

# g₁ = y₁² + u₁ p₁,   g₂ = y₂ t³,   g₃ = k y₁
function _ad_path!(g, y, ctx, t)
    k, u, p = ctx.model.k, ctx.u, ctx.params
    g[1] = y[1]^2 + u[1] * (isempty(p) ? 0.0 : p[1])
    g[2] = y[2] * t^3
    g[3] = k * y[1]
    return nothing
end

@testset "AD generators — dynamics partials against closed form" begin
    y, t = [0.7, -1.2], 1.3
    k, u, p = 2.5, [0.4, -0.9], [3.0]
    ctx = EvalContext(_AdModel(k, 0.1), u, p)

    dFy = zeros(2, 2); _make_dynamics_jac_y_ad(_ad_dynamics!, 2)(dFy, y, ctx, t)
    @test dFy ≈ [k*u[1] 0.0; 0.0 cos(y[2])] atol = 1e-14

    dFu = zeros(2, 2); _make_dynamics_jac_u_ad(_ad_dynamics!, 2)(dFu, y, ctx, t)
    @test dFu ≈ [k*y[1] 0.0; 0.0 t] atol = 1e-14

    dFp = zeros(2, 1); _make_dynamics_jac_p_ad(_ad_dynamics!, 2)(dFp, y, ctx, t)
    @test dFp ≈ reshape([t^2, 0.0], 2, 1) atol = 1e-14

    dFt = zeros(2); _make_dynamics_jac_t_ad(_ad_dynamics!, 2)(dFt, y, ctx, t)
    @test dFt ≈ [2p[1]*t, u[2]] atol = 1e-14

    dFk = zeros(2)
    _make_dynamics_jac_modelfield_ad(_ad_dynamics!, 2, AstroSolve.Accessors.@optic(_.k))(dFk, y, ctx, t)
    @test dFk ≈ [y[1]*u[1], 0.0] atol = 1e-14

    # With no parameters the parameter Jacobian is left as it was handed over.
    sentinel = fill(7.0, 2, 1)
    _make_dynamics_jac_p_ad(_ad_dynamics!, 2)(sentinel, y, EvalContext(_AdModel(k, 0.1), u), t)
    @test sentinel == fill(7.0, 2, 1)
end

@testset "AD generators — path-constraint partials against closed form" begin
    y, t = [0.7, -1.2], 1.3
    k, u, p = 2.5, [0.4, -0.9], [3.0]
    ctx = EvalContext(_AdModel(k, 0.1), u, p)

    dgy = zeros(3, 2); _make_path_jac_y_ad(_ad_path!, 3)(dgy, y, ctx, t)
    @test dgy ≈ [2y[1] 0.0; 0.0 t^3; k 0.0] atol = 1e-14

    dgu = zeros(3, 2); _make_path_jac_u_ad(_ad_path!, 3)(dgu, y, ctx, t)
    @test dgu ≈ [p[1] 0.0; 0.0 0.0; 0.0 0.0] atol = 1e-14

    dgp = zeros(3, 1); _make_path_jac_p_ad(_ad_path!, 3)(dgp, y, ctx, t)
    @test dgp ≈ reshape([u[1], 0.0, 0.0], 3, 1) atol = 1e-14

    dgt = zeros(3); _make_path_jac_t_ad(_ad_path!, 3)(dgt, y, ctx, t)
    @test dgt ≈ [0.0, 3y[2]*t^2, 0.0] atol = 1e-14

    dgk = zeros(3)
    _make_path_jac_modelfield_ad(_ad_path!, 3, AstroSolve.Accessors.@optic(_.k))(dgk, y, ctx, t)
    @test dgk ≈ [0.0, 0.0, y[1]] atol = 1e-14

    # Empty controls and parameters leave the blocks untouched.
    no_u = EvalContext(_AdModel(k, 0.1), Float64[], Float64[])
    su = fill(7.0, 3, 0); _make_path_jac_u_ad(_ad_path!, 3)(su, y, no_u, t)
    sp = fill(7.0, 3, 1); _make_path_jac_p_ad(_ad_path!, 3)(sp, y, no_u, t)
    @test size(su) == (3, 0)
    @test sp == fill(7.0, 3, 1)
end

@testset "AD generators — one model field is seeded, its siblings are not" begin
    m = _AdModel(2.0, 5.0)
    seen = Ref{Any}(nothing)
    d = _ad_field_deriv(m, AstroSolve.Accessors.@optic(_.c)) do md
        seen[] = md
        md.k * md.c^2
    end
    @test d ≈ 2 * 2.0 * 5.0 atol = 1e-14          # ∂(k c²)/∂c = 2 k c
    @test seen[].k isa Real && !(seen[].k isa ForwardDiff.Dual)
    @test seen[].c isa ForwardDiff.Dual
    @test m.c == 5.0                               # the carrier itself is not modified
end

nothing
