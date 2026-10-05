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

    @testset "no step straddles a kink on steps of an hour (GEO)" begin
        # GEO near the equinox, entering the shadow about three hours in. At 1e-9 Vern9 takes steps
        # of over an hour; a look-ahead capped at 900 s let two of them straddle both edges.
        t0 = Time("2024-03-20T00:00:00", UTC(), ISOT())
        ŝ = let s = _APS.force_position(nothing, earth, sun, t0.tt); s / norm(s) end
        z = SVector(0.0, 0.0, 1.0)
        e1 = let c = SVector(z[2]*ŝ[3] - z[3]*ŝ[2], z[3]*ŝ[1] - z[1]*ŝ[3], z[1]*ŝ[2] - z[2]*ŝ[1]); c / norm(c) end
        a = 42164.0; θ = -π/2 - 0.8
        r = a * (cos(θ) * (-ŝ) + sin(θ) * e1)
        v = sqrt(earth.mu / a) * let c = SVector(z[2]*r[3] - z[3]*r[2], z[3]*r[1] - z[1]*r[3], z[1]*r[2] - z[2]*r[1]); c / norm(c) end
        sc = Spacecraft(; state = CartesianState(vcat(r, v)), time = t0, mass = 1000.0,
                        srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
        geo = propagate!(OrbitPropagator(ForceModel(PointMassGravity(earth, (sun,)), srp),
                                         IntegratorConfig(Vern9(); reltol = 1e-9, abstol = 1e-9, dt = 60.0)),
                         sc, StopAt(sc, PropDurationSeconds(), 30 * 3600.0))
        e0 = t0.tt
        gk(i) = (out = zeros(2);
                 rr = SVector{3}(geo.u[i][1], geo.u[i][2], geo.u[i][3]);
                 _APS._kink_values!(out, 0, srp, rr, _APS.force_position(nothing, earth, sun, e0 + geo.t[i] / 86400.0)); out)
        G = [gk(i) for i in eachindex(geo.t)]
        # A kink function changes at about the orbit rate, 7.3e-5 rad/s; within 0.1 s of a kink it
        # is under 1e-5, so a crossing with both ends beyond that was stepped across.
        straddled = count(i -> any(j -> sign(G[i][j]) != sign(G[i+1][j]) &&
                                        min(abs(G[i][j]), abs(G[i+1][j])) > 1e-5, 1:2),
                          1:length(G)-1)
        @test maximum(diff(geo.t)) > 900.0          # the steps the old cap could not see across
        @test any(i -> any(j -> sign(G[i][j]) != sign(G[i+1][j]), 1:2), 1:length(G)-1)   # it crosses
        @test straddled == 0
    end

    @testset "a trial step inside the Earth is rejected, not an error" begin
        # The default initial step is 5000 s; its first trial puts a stage inside the Earth, where
        # the apparent radius of the Earth was asin of a ratio above 1, a DomainError.
        for shadow in (DualCone(), SmoothedConical())
            sc = _shadow_sc()
            f = ForceModel(PointMassGravity(earth, ()), SolarRadiationPressure(earth; shadow = shadow))
            propagate!(OrbitPropagator(f, IntegratorConfig(Vern9(); reltol = 1e-10, abstol = 1e-10)),
                       sc, StopAt(sc, PropDurationSeconds(), 6 * 3600.0))
            @test norm(to_vector(sc.state)[1:3]) > 6378.0
        end
        @test _APS._shadow_factor(DualCone(), SVector(3000.0, 0.0, 0.0), SVector(-1.496e8, 0.0, 0.0),
                                  srp.R_sun, srp.R_occ) == 0.0
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
