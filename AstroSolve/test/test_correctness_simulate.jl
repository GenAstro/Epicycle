# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Simulated tracking data.
#
# `simulate` produces the observations a known trajectory would have produced, which is what a
# study, a test or a documentation page needs and what the estimators read. The checks here are
# against things that hold by construction rather than against recorded numbers: a noiseless
# measurement equals the predictor the estimator itself uses, a seeded run repeats, a station that
# cannot see the spacecraft records nothing, and a batch fit of noiseless data returns the state
# the data was made from.
#
# The visibility check is the one worth stating. `is_visible` was two different functions of the
# same name — AstroModels' station cutoff and a permissive fallback inside Measurements — so a
# station in a signal path was never asked whether the spacecraft was above its horizon. The
# elevation test below fails against that.

using AstroSolve
using AstroSolve: simulate
using AstroModels
using AstroProp
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase
using LinearAlgebra
using OrdinaryDiffEqHighOrderRK: DP8
using Test

const _SIM_EPOCH = AstroEpochs.Time("2020-03-01T00:00:00.000", TT(), ISOT())
const _SIM_TRUTH = [6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]
const _SIM_SIG_R = 15.0e-3      # km
const _SIM_SIG_D = 2.0e-5       # km/s

_sim_station(; min_elevation = -90.0) =
    GroundStation(name = "DSS-14", body = earth, latitude = 35.4267, longitude = -116.89,
                  altitude = 1.0, min_elevation = min_elevation)

_sim_prop() = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                              IntegratorConfig(DP8(); dt = 60.0, reltol = 1e-12, abstol = 1e-12))

_sim_sc(y = _SIM_TRUTH; name = "Sat") =
    Spacecraft(state = CartesianState(copy(y)), time = _SIM_EPOCH,
               coord_sys = CoordinateSystem(earth, ICRF()), name = name)

_sim_tracking(gs, sc) = AbstractMeasurement[
    TwoWayRange(SignalPath(gs, sc, gs);   noise = MeasurementNoise(_SIM_SIG_R)),
    TwoWayDoppler(SignalPath(gs, sc, gs); noise = MeasurementNoise(_SIM_SIG_D))]

@testset "simulate — a record per measurement per epoch, in time order" begin
    gs, sc = _sim_station(), _sim_sc()
    times  = 60.0:60.0:600.0
    recs   = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), times)

    # With no elevation cutoff every epoch yields both measurements.
    @test length(recs) == 2 * length(times)
    @test count(r -> r.measurement_type === :RANGE, recs) == length(times)
    @test count(r -> r.measurement_type === :DOPPLER, recs) == length(times)
    @test all(r -> r.participant_1 == "DSS-14", recs)
    @test issorted([Float64(r.t_receive - _SIM_EPOCH) for r in recs])

    # The epochs are the ones asked for. The tolerance is microseconds because an epoch is a
    # Julian date in days, which resolves to about that in double precision.
    range_times = [Float64(r.t_receive - _SIM_EPOCH) * 86400.0
                   for r in recs if r.measurement_type === :RANGE]
    @test range_times ≈ collect(times) atol = 1e-4

    # The spacecraft handed in is left at its epoch: simulation flies a copy.
    @test Float64(sc.time - _SIM_EPOCH) == 0.0
    @test to_posvel(sc) == _SIM_TRUTH
end

@testset "simulate — a noiseless record is the predictor the estimator uses" begin
    # Truth: the observed value with `noise = false` is the same number the estimator computes for
    # the same state and epoch, so a fit of it has identically zero residual.
    gs, sc = _sim_station(), _sim_sc()
    recs = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), [300.0, 600.0]; noise = false)

    truth = _sim_sc()
    propagate!(_sim_prop(), truth, StopAt(truth, PropDurationSeconds(), 300.0))
    y = to_posvel(truth)
    specs = _sim_tracking(gs, truth)
    @test recs[1].observed ≈ AstroSolve.BatchLeastSquares._predict(specs[1], y, truth.time, truth) atol = 1e-9
    @test recs[2].observed ≈ AstroSolve.BatchLeastSquares._predict(specs[2], y, truth.time, truth) atol = 1e-12

    # Noise is added when asked for, and it is the size the measurement declares.
    noisy = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), [300.0, 600.0]; seed = 1)
    @test noisy[1].observed != recs[1].observed
    @test abs(noisy[1].observed - recs[1].observed) < 10 * _SIM_SIG_R
end

@testset "simulate — a seed repeats a run, and a different seed does not" begin
    gs, sc = _sim_station(), _sim_sc()
    a = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), 60.0:60.0:300.0; seed = 42)
    b = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), 60.0:60.0:300.0; seed = 42)
    c = simulate(_sim_tracking(gs, sc), sc, _sim_prop(), 60.0:60.0:300.0; seed = 43)
    @test [r.observed for r in a] == [r.observed for r in b]
    @test [r.observed for r in a] != [r.observed for r in c]
end

@testset "simulate — a station records only what it can see" begin
    # One revolution from one station: a cutoff of 5° leaves a single pass, and a cutoff no orbit
    # clears leaves nothing. Against the permissive fallback both would record every epoch.
    sc    = _sim_sc()
    times = 60.0:60.0:5400.0
    open_sky = simulate(_sim_tracking(_sim_station(), sc), sc, _sim_prop(), times)
    masked   = simulate(_sim_tracking(_sim_station(min_elevation = 5.0), sc), sc, _sim_prop(), times)
    closed   = simulate(_sim_tracking(_sim_station(min_elevation = 89.0), sc), sc, _sim_prop(), times)

    @test length(open_sky) == 2 * length(times)
    @test 0 < length(masked) < length(open_sky)
    @test isempty(closed)

    # What survives the cutoff is contiguous in time: a pass, not scattered epochs.
    t_pass = sort!(unique([Float64(r.t_receive - _SIM_EPOCH) * 86400.0 for r in masked]))
    @test all(isapprox.(diff(t_pass), 60.0; atol = 1e-4))
end

@testset "simulate — noiseless data fits back to the state it was made from" begin
    # End to end, against the estimator: a batch fit of noiseless single-station data returns the
    # truth state, from a guess a kilometre and a metre per second away.
    gs    = _sim_station(min_elevation = 5.0)
    truth = _sim_sc(; name = "Sat (truth)")
    recs  = simulate(_sim_tracking(gs, truth), truth, _sim_prop(), 30.0:30.0:21600.0; noise = false)
    @test length(recs) > 50

    guess = _SIM_TRUTH .+ [1.0, -1.0, 0.5, 1e-3, -1e-3, 5e-4]
    sat   = _sim_sc(guess; name = "Sat (estimate)")
    y0    = Vary(state, sat; guess = guess,
                 covariance = Diagonal([1e2, 1e2, 1e2, 1e-2, 1e-2, 1e-2]))
    problem = ODProblem(spacecraft = sat, propagator = _sim_prop(),
                        measurements = _sim_tracking(gs, sat), solve_for = [y0])

    result = solve!(problem, recs; method = Batch(n_iters = 10, tol = 1e-9))
    @test maximum(abs, result.X_hat[1:3] .- _SIM_TRUTH[1:3]) < 1e-3     # a metre
    @test maximum(abs, result.X_hat[4:6] .- _SIM_TRUTH[4:6]) < 1e-6     # a millimetre per second

    # The ODProblem carries the same three arguments, so it can be simulated from directly.
    from_problem = simulate(problem, 30.0:30.0:600.0; noise = false)
    direct = simulate(problem.measurements, problem.spacecraft, problem.propagator,
                      30.0:30.0:600.0; noise = false)
    @test [r.observed for r in from_problem] == [r.observed for r in direct]
end

@testset "simulate — bad input is refused where it is found" begin
    gs, sc = _sim_station(), _sim_sc()
    specs  = _sim_tracking(gs, sc)

    @test_throws ArgumentError simulate(AbstractMeasurement[], sc, _sim_prop(), [60.0])
    @test_throws ArgumentError simulate(specs, sc, _sim_prop(), Float64[])
    @test_throws ArgumentError simulate(specs, sc, _sim_prop(), [120.0, 60.0])
    @test_throws ArgumentError simulate(specs, sc, _sim_prop(), [-60.0, 60.0])

    msg = try
        simulate(specs, sc, _sim_prop(), [120.0, 60.0])
    catch e
        sprint(showerror, e)
    end
    @test occursin("ascend", msg)
end
