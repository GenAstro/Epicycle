# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Force partials: the state Jacobian of each force with a state_jac! of its own, against
# ForwardDiff through the force's whole accel_eval!, the fallback it replaces.
#
# HarmonicGravity builds ∂a/∂r = Rᵀ G R from the field's gravity gradient G in its own axes. The
# zonal field has no analytic gradient, so G comes from geopotential_gradient's default, ForwardDiff
# in the three body-fixed position components; that default is checked separately against central
# differences, so the reference is not AD alone. The Jacobian is evaluated at states at low and high
# latitude and at GEO.
#
# Measured on 2026-10-06: zonal J2–J4 Jacobian vs full AD 2.2e-16; the default gradient vs central
# differences about 1e-8, the differences' own truncation and round-off.

using Test
using LinearAlgebra
using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth, sun
using StaticArrays: SVector

@testset "Force partials: HarmonicGravity state Jacobian" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    sc = Spacecraft(; state = CartesianState([6878.137, 0, 0, 0, 7.6, 0.0]), time = t, mass = 1000.0)
    states = ([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820], [-3000.0, 5000.0, 4200.0, -5.0, -3.0, 4.0],
              [0.0, 0.0, 7000.0, 7.5, 0.0, 0.0], [42164.0, 100.0, -50.0, 0.0, 3.07, 0.0])
    f = HarmonicGravity(earth; degree = 4, order = 0, model = Zonal())

    @testset "against ForwardDiff through accel_eval!" begin
        for y in states
            A = zeros(6, 6); AstroProp.state_jac!(A, f, t, y, sc)
            A_ad = zeros(6, 6); AstroProp._state_jac_ad!(A_ad, f, t, y, sc)
            @test norm(A[4:6, 1:3] - A_ad[4:6, 1:3]) / norm(A_ad[4:6, 1:3]) < 1e-13
            @test A[4:6, 4:6] == zeros(3, 3)                 # gravity does not depend on velocity
        end
    end

    @testset "geopotential_gradient's default against central differences" begin
        for r in (SVector(6878137.0, 0.0, 0.0), SVector(-3.0e6, 5.0e6, 4.2e6), SVector(1.0, -2.0, 7.0e6))
            G = geopotential_gradient(f.model, f.data, r, 0.0, f.degree, f.order)
            G_fd = zeros(3, 3)
            for j in 1:3
                h = 1.0                                       # m
                e = SVector(ntuple(k -> k == j ? h : 0.0, 3))
                G_fd[:, j] = (geopotential_accel(f.model, f.data, r + e, 0.0, f.degree, f.order) -
                              geopotential_accel(f.model, f.data, r - e, 0.0, f.degree, f.order)) / 2h
            end
            @test norm(G - G_fd) / norm(G) < 1e-6
        end
    end
end

# AtmosphericDrag: ∂a/∂v and the rotation's part of ∂a/∂r are closed form; the density's gradient
# comes from density_gradient, by default one ForwardDiff evaluation of the density alone, checked
# here against central differences. Measured on 2026-10-06 for the exponential atmosphere: Jacobian
# vs full AD 2.9e-16 (∂a/∂r) and 2.5e-16 (∂a/∂v); gradient vs central differences about 1e-9.
@testset "Force partials: AtmosphericDrag state Jacobian" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    sc = Spacecraft(; state = CartesianState([6578.137, 0, 0, 0, 7.7, 0.0]), time = t, mass = 1000.0,
                    drag = SphericalDrag(c_d = 2.2, drag_area = 10.0))
    f = AtmosphericDrag(earth; model = Exponential())
    states = ([6578.137, 120.0, 300.0, 0.0, 7.7, 0.3], [-3000.0, 4500.0, 3600.0, -5.0, -3.5, 4.0],
              [100.0, -200.0, 6700.0, 7.6, 0.2, 0.0])
    for y in states
        A = zeros(6, 6); AstroProp.state_jac!(A, f, t, y, sc)
        B = zeros(6, 6); AstroProp._state_jac_ad!(B, f, t, y, sc)
        @test norm(A[4:6, 1:3] - B[4:6, 1:3]) / norm(B[4:6, 1:3]) < 1e-13
        @test norm(A[4:6, 4:6] - B[4:6, 4:6]) / norm(B[4:6, 4:6]) < 1e-13
    end
    jd = AstroProp.force_epoch(nothing, t).utc
    R, _ = AstroProp._rotation_blocks(AstroProp.force_rotation(nothing, f.orientation, 399, t))
    r = SVector(6578.137, 120.0, 300.0)
    ρ, ∇ρ = density_gradient(f.model, jd, r, R)
    @test ρ == density(f.model, jd, r, R)
    h = 1.0e-3
    ∇fd = [(density(f.model, jd, r + h * SVector(ntuple(k -> k == j, 3)), R) -
            density(f.model, jd, r - h * SVector(ntuple(k -> k == j, 3)), R)) / 2h for j in 1:3]
    @test norm(∇ρ - ∇fd) / norm(∇ρ) < 1e-6
end

# SolarRadiationPressure: closed form but for the shadow factor's gradient, by ForwardDiff on the
# shadow function alone. Points in sunlight, in the penumbra and in the umbra, with both shadows,
# and about the Sun, where there is no shadow. Measured on 2026-10-06: ≤ 6.3e-16 against full AD.
@testset "Force partials: SolarRadiationPressure state Jacobian" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    sc = Spacecraft(; state = CartesianState([6878.137, 0, 0, 0, 7.6, 0.0]), time = t, mass = 1000.0,
                    srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
    rs = AstroProp.force_position(nothing, earth, sun, t); ŝ = rs / norm(rs)
    n̂ = normalize(cross(ŝ, SVector(0.0, 0.0, 1.0)))
    at(θ) = [6878.137 * (-cos(θ) * ŝ + sin(θ) * n̂); 0.0; 7.6; 0.0]          # θ from the anti-Sun line
    for shadow in (DualCone(), SmoothedConical()), θ in (1.6, 1.188, 0.3)  # sunlight, penumbra, umbra
        f = SolarRadiationPressure(earth; shadow = shadow)
        y = at(θ)
        A = zeros(6, 6); AstroProp.state_jac!(A, f, t, y, sc)
        B = zeros(6, 6); AstroProp._state_jac_ad!(B, f, t, y, sc)
        @test norm(A[4:6, 1:3] - B[4:6, 1:3]) <= 1e-13 * max(norm(B[4:6, 1:3]), eps())
        @test A[4:6, 4:6] == zeros(3, 3)
    end
    fs = SolarRadiationPressure(sun)
    y = [1.5e8, 1.0e6, 0.0, 0.0, 30.0, 0.0]
    A = zeros(6, 6); AstroProp.state_jac!(A, fs, t, y, sc)
    B = zeros(6, 6); AstroProp._state_jac_ad!(B, fs, t, y, sc)
    @test norm(A[4:6, 1:3] - B[4:6, 1:3]) / norm(B[4:6, 1:3]) < 1e-13
end

# eval_accel_and_jacobian!: one pass, a force with accel_and_jac! evaluated once for both. It must
# give accel_eval!'s acceleration and eval_jacobian!'s Jacobian to round-off: a density evaluated on
# dual numbers can differ from the plain evaluation in its last bit (one unit, 4e-20 km/s², here).
@testset "eval_accel_and_jacobian! against separate calls" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    sc = Spacecraft(; state = CartesianState([6578.137, 0, 0, 0, 7.7, 0.0]), time = t, mass = 1000.0,
                    drag = SphericalDrag(c_d = 2.2, drag_area = 10.0), srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
    fm = ForceModel(HarmonicGravity(earth; degree = 4, order = 0, model = Zonal()),
                    PointMassGravity(earth, (sun,); include_center = false),
                    AtmosphericDrag(earth; model = Exponential()), SolarRadiationPressure(earth))
    y = [6578.137, 120.0, 300.0, 0.0, 7.7, 0.3]
    cfg = JacobianConfig(partial_y = true)
    r1 = JacobianResult(cfg; n_state = 6, forces = fm, sc = sc, t_example = t)
    r2 = JacobianResult(cfg; n_state = 6, forces = fm, sc = sc, t_example = t)
    @test r1.plan.combined == Set([typeof(fm.forces[1]), typeof(fm.forces[3])])
    acc = zeros(6); eval_accel_and_jacobian!(acc, r1, fm, y, sc, t)
    ctx = AstroProp.ForceContext(); AstroProp._reset!(ctx, t)
    acc0 = zeros(6); AstroProp._eval_all!(fm, t, y, acc0, sc, (context = ctx,))
    eval_jacobian!(r2, fm, y, sc, t)
    @test acc ≈ acc0 rtol = 1e-14
    @test r1.partial_y ≈ r2.partial_y rtol = 1e-14
end
