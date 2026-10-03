# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Steps ending at the shadow boundaries (src/discontinuities.jl).
#
# With DualCone the SRP acceleration has a kink at each penumbra edge, and the propagator ends a
# step at each one rather than stepping across it. The tests check the mechanism and its purpose:
#
#   - no step straddles a kink: wherever a kink function changes sign between two saved steps, one
#     of them is within 0.1 s of the kink;
#   - the integration error at an ordinary tolerance is what the tolerance buys on a smooth
#     problem, against a converged run (Vern9, 1e-12, 2 s maximum step);
#   - a force model without kinks gets no callback, and DualCone's lighting factor stays defined
#     on the edges a step now lands on.
#
# Measured on 2026-10-02, one day of LEO under point-mass gravity and SRP, Vern9 at 1e-11:
# 0.25 mm from the converged run, which the 2 mm tolerance is about eight times.

using Test
using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth, sun
using OrdinaryDiffEqVerner: Vern9
using LinearAlgebra: norm
using StaticArrays: SVector

const _APS = AstroProp

_shadow_sc() = Spacecraft(; state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                          time = Time("2020-10-20T12:00:00", UTC(), ISOT()), mass = 1000.0,
                          srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))

function _shadow_fly(srp, tol; dtmax = Inf)
    sc = _shadow_sc()
    forces = ForceModel(PointMassGravity(earth, ()), srp)
    sol = propagate!(OrbitPropagator(forces, IntegratorConfig(Vern9(); reltol = tol, abstol = tol,
                                                              dt = 60.0)),
                     sc, StopAt(sc, PropDurationSeconds(), 86400.0); dtmax = dtmax)
    return sol
end

@testset "Shadow boundaries end a step" begin
    srp = SolarRadiationPressure(earth; shadow = DualCone())
    sol = _shadow_fly(srp, 1e-11)
    epoch0 = _shadow_sc().time.tt

    @testset "no step straddles a kink" begin
        g(i) = (out = zeros(2);
                t = epoch0 + sol.t[i] / 86400.0;
                r = SVector{3}(sol.u[i][1], sol.u[i][2], sol.u[i][3]);
                _APS._kink_values!(out, 0, srp, r, _APS.force_position(nothing, earth, sun, t)); out)
        # A kink function changes by about the orbit rate, 1.1e-3 rad/s; 0.1 s of it is 1.1e-4.
        near = 1.1e-4
        crossings = 0; straddled = 0
        gprev = g(1)
        for i in 2:length(sol.t)
            gi = g(i)
            for j in 1:2
                if sign(gprev[j]) != sign(gi[j])
                    crossings += 1
                    min(abs(gprev[j]), abs(gi[j])) < near || (straddled += 1)
                end
            end
            gprev = gi
        end
        @test crossings ≥ 4 * 14          # four edges an orbit, about 15 orbits
        @test straddled == 0
    end

    @testset "error at an ordinary tolerance" begin
        y  = sol.u[end]
        yc = _shadow_fly(srp, 1e-12; dtmax = 2.0).u[end]
        err = norm(y[1:3] - yc[1:3]) * 1e3                  # m
        println("shadow events: Vern9 1e-11 is ", err * 1e3, " mm from the converged run")
        @test err < 2e-3
    end

    @testset "no callback without kinks" begin
        @test _APS._kink_callback(ForceModel(PointMassGravity(earth, ())), epoch0, (1:6,)) === nothing
        smooth = SolarRadiationPressure(earth; shadow = SmoothedConical())
        @test _APS._kink_callback(ForceModel(PointMassGravity(earth, ()), smooth), epoch0, (1:6,)) === nothing
    end

    @testset "DualCone defined on its edges" begin
        # Sweep the separation through both penumbra edges in steps far finer than round-off
        # resolves near them, where the lighting factor used to take sqrt of a negative number.
        r_sun = SVector(1.496e8, 0.0, 0.0)
        R_sun, R_occ = srp.R_sun, srp.R_occ
        ok = true
        for c in range(1.17, 1.21; length = 200_001)       # the edges are near 1.19 rad from anti-Sun
            r = 6878.137 * SVector(-cos(c), sin(c), 0.0)
            F = _APS._shadow_factor(DualCone(), r, r_sun, R_sun, R_occ)
            ok &= isfinite(F) && 0 ≤ F ≤ 1
        end
        @test ok
    end
end
