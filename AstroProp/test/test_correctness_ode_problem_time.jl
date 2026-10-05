# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# OrbitODEProblem integrates in the time scale propagate! does: TT about the Earth, TDB about any
# other body. It used to integrate in the spacecraft's own scale, usually UTC, which is not uniform
# across a leap second and is not the dynamics' time about another body.

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroUniverse
using LinearAlgebra: norm
using Test

# The two paths step differently, so they agree only to integration noise, about 1e-12 of the
# radius (1.4e-9 km over a day about the Moon). Reading the ephemeris at UTC instead, 69.184 s
# off, moves the Earth case by 4.4 mm and the Moon case by 26 cm, so 1e-7 km sits well between.
const _SAME = 1e-7                                     # km

function _both(sc_of, forces, duration)
    prop = OrbitPropagator(ForceModel(forces...),
                           IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0))
    s1 = sc_of()
    y1 = Vector(propagate!(prop, s1, StopAt(s1, PropDurationSeconds(), duration)).u[end])
    y2 = solve(OrbitODEProblem(prop, sc_of(); duration_s = duration)).y_final
    return y1, y2
end

@testset "OrbitODEProblem — the propagator's time scale" begin
    @testset "about the Earth, across a leap second (TT)" begin
        # 2016-12-31T23:59:60 UTC: a UTC clock would lose a second here.
        sc() = Spacecraft(state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                          time = Time("2016-12-31T22:00:00", UTC(), ISOT()))
        y1, y2 = _both(sc, (PointMassGravity(earth, (moon, sun)),), 4 * 3600.0)
        @test norm(y1[1:3] - y2[1:3]) < _SAME
    end

    @testset "about the Moon (TDB)" begin
        sc() = Spacecraft(state = CartesianState([1837.4, 0.0, 0.0, 0.0, 1.6335, 0.0]),
                          time = Time("2020-10-20T12:00:00", UTC(), ISOT()),
                          coord_sys = CoordinateSystem(moon, ICRF()))
        y1, y2 = _both(sc, (PointMassGravity(moon, (earth, sun)),), 86400.0)
        @test norm(y1[1:3] - y2[1:3]) < _SAME
    end
end
