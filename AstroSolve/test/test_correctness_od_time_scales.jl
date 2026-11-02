# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Orbit determination with the spacecraft's epoch in any time scale, over an arc that crosses a leap
# second.
#
# The truth is the definition: the same instant written in TT, TDB or UTC is the same problem, so a
# fit must not depend on which one the user wrote. The arc crosses the 2016-12-31 leap second on
# purpose. Counting elapsed time in UTC there would put the later records a full second out, about
# 7 km along track, and the fits would not agree.

using Test
using AstroEpochs
using AstroUniverse
using AstroStates
using AstroFrames
using AstroModels
using AstroProp
using AstroSolve
using LinearAlgebra
using OrdinaryDiffEqHighOrderRK: DP8

const _SCALES_Y = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]
const _SCALES_GUESS = _SCALES_Y .+ [0.5, -0.5, 0.25, 5e-4, -5e-4, 2.5e-4]
const _SCALES_EPOCH_TT = Time("2016-12-31T23:30:00.000", UTC(), ISOT()).tt

_scales_station() = GroundStation(name = "DSS-14", body = earth, latitude = 35.4267, longitude = -116.89,
                                  altitude = 1.0, min_elevation = -90.0)
_scales_propagator() = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                                       IntegratorConfig(DP8(); dt = 60.0, reltol = 1e-11, abstol = 1e-11))
_scales_tracking(station, sc) = AbstractMeasurement[
    TwoWayRange(SignalPath(station, sc, station);   noise = MeasurementNoise(15.0e-3)),
    TwoWayDoppler(SignalPath(station, sc, station); noise = MeasurementNoise(2.0e-5))]

function _scales_fit(epoch, records, station)
    sat = Spacecraft(state = CartesianState(copy(_SCALES_GUESS)), time = epoch, name = "estimate")
    y0  = Vary(state, sat; guess = _SCALES_GUESS,
               covariance = Diagonal([1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]))
    problem = ODProblem(spacecraft = sat, propagator = _scales_propagator(),
                        measurements = _scales_tracking(station, sat), solve_for = [y0])
    return solve!(problem, copy(records); method = Batch(n_iters = 6, tol = 1e-9, verbose = false))
end

@testset "OD — the epoch's time scale does not change the fit" begin
    station = _scales_station()
    truth   = Spacecraft(state = CartesianState(copy(_SCALES_Y)), time = _SCALES_EPOCH_TT, name = "truth")
    records = simulate(_scales_tracking(station, truth), truth, _scales_propagator(), 120.0:120.0:3600.0;
                       seed = 42)
    @test records[1].t_receive.tt.isot < "2017-01-01T00:00:00.000" < records[end].t_receive.tt.isot

    base = _scales_fit(_SCALES_EPOCH_TT, records, station)
    @test norm(base.X_hat[1:3] - _SCALES_Y[1:3]) < 1.0

    for (name, epoch) in (("TDB", _SCALES_EPOCH_TT.tdb), ("UTC", _SCALES_EPOCH_TT.utc), ("TAI", _SCALES_EPOCH_TT.tai))
        fit = _scales_fit(epoch, records, station)
        @test norm(fit.X_hat[1:3] - base.X_hat[1:3]) < 1e-6       # km
        @test norm(fit.X_hat[4:6] - base.X_hat[4:6]) < 1e-9       # km/s
    end
end

@testset "phase time — seconds from a UTC epoch are elapsed seconds, across a leap second" begin
    a = Time("2016-12-31T23:59:00.000", UTC(), ISOT())
    @test AstroSolve._epoch_after(a, 121.0).isot == "2017-01-01T00:01:00.000"
    @test AstroSolve._epoch_after(a, 60.5).isot == "2016-12-31T23:59:60.500"
    @test AstroSolve._epoch_after(a, 0.0).scale === UTC()
    tt = a.tt
    @test AstroSolve._epoch_after(tt, 86400.0) == tt + 1.0
end

@testset "OD — simulate takes epochs in any scale, and any epoch scale for the spacecraft" begin
    station = _scales_station()
    for epoch in (_SCALES_EPOCH_TT, _SCALES_EPOCH_TT.tdb, _SCALES_EPOCH_TT.utc)
        truth = Spacecraft(state = CartesianState(copy(_SCALES_Y)), time = epoch, name = "truth")
        by_seconds = simulate(_scales_tracking(station, truth), truth, _scales_propagator(),
                              [600.0, 3000.0]; noise = false)
        by_epoch   = simulate(_scales_tracking(station, truth), truth, _scales_propagator(),
                              [(_SCALES_EPOCH_TT + 600.0 / 86400).utc, (_SCALES_EPOCH_TT + 3000.0 / 86400).tdb];
                              noise = false)
        @test length(by_seconds) == length(by_epoch) == 4
        for (a, b) in zip(by_seconds, by_epoch)
            @test a.observed ≈ b.observed rtol = 1e-9
            @test abs(a.t_receive.tt - b.t_receive.tt) * 86400 < 1e-5
        end
    end
end

nothing
