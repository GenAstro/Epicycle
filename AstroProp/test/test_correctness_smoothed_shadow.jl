# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# SmoothedConical shadow, against DualCone.
#
# There is no external truth for a smoothed eclipse: GMAT and Orekit model the penumbra exactly. The
# truth is DualCone, which force_srp_spherical.jl checks against GMAT. SmoothedConical is held to it
# in three ways: the lighting factor across the shadow edge at LEO and GEO, the derivative of the
# SRP acceleration through the penumbra, and a one-day LEO propagation.

using Test
using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth, sun
using OrdinaryDiffEqVerner: Vern9
using ForwardDiff
using LinearAlgebra: norm
using StaticArrays: SVector

const _AP = AstroProp

@testset "SmoothedConical shadow vs DualCone" begin

    R_sun = sun.equatorial_radius
    R_occ = earth.equatorial_radius
    r_sun = SVector(1.496e8, 0.0, 0.0)     # Sun along +x, so the shadow lies along −x

    # A point at radius `rad` behind the Earth, `y` km off the shadow axis.
    behind(rad, y) = SVector(-sqrt(rad^2 - y^2), y, 0.0)
    F(model, r) = _AP._shadow_factor(model, r, r_sun, R_sun, R_occ)

    @testset "construction" begin
        @test SmoothedConical().sharpness == 3.25
        @test SmoothedConical(sharpness = 5).sharpness == 5.0
        @test_throws ArgumentError SmoothedConical(sharpness = 0.0)
        @test_throws ArgumentError SmoothedConical(sharpness = -1.0)
    end

    @testset "lighting factor at $(name)" for (name, rad) in (("LEO", 7000.0), ("GEO", 42164.0))
        sm = SmoothedConical()
        worst = 0.0
        for y in range(R_occ - 1500.0, min(R_occ + 1500.0, rad - 1.0); length = 3001)
            r  = behind(rad, y)
            fs = F(sm, r)
            @test 0.0 ≤ fs ≤ 1.0
            worst = max(worst, abs(fs - F(DualCone(), r)))
        end
        # Measured 0.039 at LEO and 0.042 at GEO with the default sharpness.
        @test worst < 0.045

        # Away from the penumbra the two agree: full sun on the lit side, none in the umbra.
        @test F(sm, SVector(rad, 0.0, 0.0)) ≈ 1.0 atol = 1e-12
        @test F(sm, SVector(-rad, 0.0, 0.0)) ≈ 0.0 atol = 1e-12
    end

    @testset "smooth derivative through the penumbra" begin
        sc = Spacecraft(; state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                        time = Time("2020-10-20T12:00:00", UTC(), ISOT()), mass = 100.0,
                        srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
        srp = SolarRadiationPressure(earth; shadow = SmoothedConical())
        t   = sc.time
        r_sun_true = SVector{3}(_AP.translate(earth, sun, t.tdb.jd))
        ŝ = r_sun_true / norm(r_sun_true)
        n̂ = SVector(-ŝ[2], ŝ[1], 0.0) / norm(SVector(-ŝ[2], ŝ[1], 0.0))

        # The SRP acceleration as a function of the offset y from the shadow axis, at 7000 km.
        function accel(y)
            r = -sqrt(7000.0^2 - y^2) * ŝ + y * n̂
            x = [r[1], r[2], r[3], 0.0, 7.5, 0.0]
            dx = zeros(eltype(x), 6)
            accel_eval!(srp, t, x, dx, sc, nothing)
            return dx[4:6]
        end

        # ForwardDiff against a central difference at points across the penumbra and either side.
        for y in R_occ .+ (-120.0, -40.0, -10.0, 0.0, 10.0, 40.0, 120.0)
            ad = ForwardDiff.derivative(accel, y)
            h  = 1e-3
            fd = (accel(y + h) - accel(y - h)) / (2h)
            @test norm(ad - fd) ≤ 1e-6 * norm(ad) + 1e-22
            @test all(isfinite, ad)
        end
        # Deep in the umbra the logistic is far out on its tail; the derivative must stay finite.
        @test all(isfinite, ForwardDiff.derivative(accel, 0.0))
    end

    @testset "one-day LEO propagation" begin
        function propagate_with(shadow)
            sc = Spacecraft(; state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                            time = Time("2020-10-20T12:00:00", UTC(), ISOT()), mass = 100.0,
                            name = "LEO", srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
            forces = ForceModel(PointMassGravity(earth, ()), SolarRadiationPressure(earth; shadow))
            prop = OrbitPropagator(forces,
                                   IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0))
            return propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 86400.0)).u[end]
        end
        y_dual   = propagate_with(DualCone())
        y_smooth = propagate_with(SmoothedConical())
        # SRP moves this orbit 278 m in a day; the two shadows differ by 7.6 cm (measured).
        @test norm(y_dual[1:3] - y_smooth[1:3]) * 1e3 < 0.15
    end
end
