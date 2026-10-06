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
using AstroUniverse: earth
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
