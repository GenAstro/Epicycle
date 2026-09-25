# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Simulated tracking data.
#
# An estimator needs observations, and the observations a test, a study or a documentation page
# needs are the ones a known trajectory would have produced. Everything the calculation takes was
# already here — the measurement predictors, the noise draw, the visibility predicate and the
# record type — and assembling them was left to the caller, who had to reach into the estimator's
# internals to do it. `simulate` is that assembly, and its records are what the estimators and
# `write_records` already read.

"""
    simulate(measurements, spacecraft, propagator, times; noise = true, seed = nothing)
    simulate(problem::ODProblem, times; noise = true, seed = nothing)

Tracking data a known trajectory would have produced.

`spacecraft` carries the truth state and its epoch; it is propagated on a copy, so the caller's
spacecraft is left where it was. At each requested epoch every measurement whose participants can
see the spacecraft is evaluated, and the value is recorded with a draw from that measurement's
noise added.

# Arguments
- `measurements`: the measurements to take, each carrying its own `SignalPath` and noise, such as
  `TwoWayRange(SignalPath(station, sat, station); noise = MeasurementNoise(0.015))`.
- `spacecraft`: the truth spacecraft, whose epoch is the start of the arc.
- `propagator`: the propagator flown between epochs.
- `times`: when to observe, either seconds past the spacecraft's epoch or epochs, ascending.
- `noise`: `false` records the computed value, with no draw added.
- `seed`: an integer seeds the draws, so a run repeats. Omitted, the global generator is used.

# Notes
A measurement is skipped, rather than recorded, when any participant reports the spacecraft out of
sight: a `GroundStation` applies its `min_elevation` cutoff, and a participant type that says
nothing about visibility sees everything. A pass therefore yields fewer records than epochs asked
for, and an arc that is never in view yields none.

An epoch is recorded as the propagation stopped on it, in the spacecraft's own time scale. It
carries a Julian date where a record read from a tracking data file carries a calendar string; the
two mix freely, because an estimator works in elapsed time from its own epoch.

Truth and estimate stay separate. Simulate from the truth spacecraft, then estimate with another
one built from the guess; handing the same spacecraft to both fits data to itself.

# Returns
A `Vector{ObservationRecord}` in time order, each carrying its measurement type, epoch, observed
value and the first participant's name. Pass it to [`solve!`](@ref) with `Batch` or `Sequential`,
or to [`write_records`](@ref) to write a tracking data file.

Throws `ArgumentError` when no measurement is given, when `times` is not ascending, or when a
measurement has no tracking-data type to record it under.

# Example
```julia
using Epicycle

station = GroundStation(name = "DSS-14", body = earth, latitude = 35.4267,
                        longitude = -116.89, altitude = 1.0, min_elevation = 5.0)

truth = Spacecraft(state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
                   time  = Time("2020-03-01T00:00:00.000", TT(), ISOT()),
                   coord_sys = CoordinateSystem(earth, ICRF()), name = "Sat")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); dt = 60.0, reltol = 1e-12, abstol = 1e-12))

measurements = [TwoWayRange(SignalPath(station, truth, station);
                            noise = MeasurementNoise(0.015)),
                TwoWayDoppler(SignalPath(station, truth, station);
                              noise = MeasurementNoise(2.0e-5))]

# One orbit at a one-minute cadence, the passes Goldstone can see.
records = simulate(measurements, truth, prop, 60.0:60.0:5400.0; seed = 42)
```
"""
function simulate(measurements::AbstractVector, sc::Spacecraft, prop::OrbitPropagator,
                  times::AbstractVector; noise::Bool = true, seed = nothing)
    isempty(measurements) && throw(ArgumentError(
        "simulate: at least one measurement is required; got none."))
    isempty(times) && throw(ArgumentError(
        "simulate: at least one observation epoch is required; got none."))

    truth   = deepcopy(sc)
    offsets = _simulate_offsets(times, truth.time)
    rng     = seed === nothing ? Random.default_rng() : Random.Xoshiro(seed)

    records = ObservationRecord[]
    flown   = 0.0
    for dt in offsets
        dt > flown && propagate!(prop, truth, StopAt(truth, PropDurationSeconds(), dt - flown))
        flown = dt
        t = truth.time
        y = to_posvel(truth)
        r = @view y[1:3]
        for m in measurements
            all(p -> is_visible(p, r, t), m.path.participants) || continue
            value = BatchLeastSquares._predict(m, y, truth.time, truth)
            noise && (value += draw(m.noise, rng))
            push!(records, ObservationRecord(_record_type(m), t, value, m.path.names[1]))
        end
    end
    return records
end

simulate(problem::ODProblem, times::AbstractVector; kwargs...) =
    simulate(problem.measurements, problem.spacecraft, problem.propagator, times; kwargs...)

# Seconds past the epoch, whichever way the caller asked. A `Time` difference is in days.
_simulate_offsets(times::AbstractVector{<:Real}, ::AstroEpochs.Time) =
    _ascending(collect(float.(times)))
_simulate_offsets(times::AbstractVector{<:AstroEpochs.Time}, epoch::AstroEpochs.Time) =
    _ascending([Float64(t - epoch) * 86400.0 for t in times])

function _ascending(offsets::Vector{Float64})
    issorted(offsets) || throw(ArgumentError(
        "simulate: observation epochs must ascend; they are propagated through in the order given."))
    offsets[1] >= 0.0 || throw(ArgumentError(
        "simulate: observation epochs start $(-offsets[1]) s before the spacecraft's epoch; " *
        "the arc is flown forward from it."))
    return offsets
end

# The tracking-data type each measurement is recorded under. `read_records` maps these back the
# other way, in BatchLeastSquares.
_record_type(::TwoWayRange)   = :RANGE
_record_type(::TwoWayDoppler) = :DOPPLER
_record_type(m::AbstractMeasurement) = throw(ArgumentError(
    "simulate: $(nameof(typeof(m))) has no tracking-data type to record it under; " *
    "two-way range and two-way Doppler have."))
