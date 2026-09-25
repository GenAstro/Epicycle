# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# Measurements are data-only definitions. Propagation and measurement prediction
# remain with the estimation problem that owns the trajectory.

module Measurements

using Random
# `is_visible` is AstroModels': a station's elevation cutoff is its own. Importing it is what makes
# the fallback below a method rather than a second function of the same name, which is what it was
# until 2026-09-19 — so a `GroundStation` in a signal path was never asked whether it could see.
import AstroModels: is_visible

# =============================================================================
# 1. SignalPath
# =============================================================================

"""
    SignalPath(participants...)

An ordered route through the participants in a tracking signal.

# Fields
- `participants::NTuple{N,Any}`: Participants in signal-hop order.
- `names::NTuple{N,String}`: Participant names used as registry and tracking-data keys.

# Notes
Each participant must have a `name` field convertible to `String`, or be a string.
The path must contain at least one participant, and adjacent participants must have
different names. A participant may appear again later in the path, as in the
round trip `SignalPath(station, spacecraft, station)`.

Throws an `ArgumentError` for an empty path or adjacent participants with the
same name.

# Example
```julia
using AstroSolve
station = (name = "DSS-14",)
sat = (name = "Sat",)
SignalPath(station, sat, station)
```
"""
struct SignalPath{N}
    participants::NTuple{N,Any}
    names::NTuple{N,String}
end

function SignalPath(participants...)
    n = length(participants)
    n >= 1 || throw(ArgumentError(
        "SignalPath: a path needs at least one participant; got none"))
    names = ntuple(i -> _participant_name(participants[i]), n)
    for i in 2:n
        names[i] == names[i-1] &&
            throw(ArgumentError(
                "SignalPath: adjacent participants must differ, since a participant " *
                "cannot relay to itself; positions $(i-1) and $(i) are both " *
                "$(repr(names[i]))"))
    end
    return SignalPath{n}(participants, names)
end

# Bare strings support paths reconstructed from tracking-data files.
_participant_name(obj) = String(getfield(obj, :name))
_participant_name(s::AbstractString) = String(s)

Base.length(::SignalPath{N}) where {N} = N
Base.getindex(p::SignalPath, i::Integer) = p.participants[i]
Base.iterate(p::SignalPath, args...) = iterate(p.participants, args...)

"""
    participant_names(path::SignalPath)

Return the participant names in signal-hop order.

# Arguments
- `path::SignalPath`: Signal route to inspect.

# Returns
An `NTuple{N,String}` containing the registry name of each participant.

# Example
```jldoctest
using AstroSolve
path = SignalPath("DSS-14", "Sat", "DSS-14")
AstroSolve.participant_names(path)

# output
("DSS-14", "Sat", "DSS-14")
```
"""
participant_names(p::SignalPath) = p.names

# =============================================================================
# 2. Measurement noise (defined before observable tags so specs can
#    carry it as a field — F4: noise rides on the function)
# =============================================================================

"""
    AbstractMeasurementNoise

Extension interface for noise models attached to measurements.

# Notes
A concrete subtype implements `draw(noise, rng)` and `variance(noise)`.
Scalar models return one residual sample and one variance. Future vector-valued
models may return a covariance matrix.

# Example
```jldoctest
using AstroSolve
MeasurementNoise(0.015) isa AstroSolve.AbstractMeasurementNoise

# output
true
```
"""
abstract type AbstractMeasurementNoise end

"""
    MeasurementNoise(sigma)

Zero-mean Gaussian noise for a scalar measurement.

# Fields
- `sigma::Float64`: Standard deviation in the measurement's units.

# Notes
`sigma` must be nonnegative. A value of zero produces deterministic,
noise-free measurements. For a range reported in kilometres, 15 m of noise is
specified as `0.015`.

Throws a `DomainError` when `sigma` is negative.

# Example
```julia
using AstroSolve
MeasurementNoise(15.0e-3)
```
"""
struct MeasurementNoise <: AbstractMeasurementNoise
    sigma::Float64

    function MeasurementNoise(sigma::Real)
        σ = Float64(sigma)
        σ >= 0.0 || throw(DomainError(σ,
            "MeasurementNoise: sigma is a standard deviation and must be " *
            "non-negative; got $(σ)"))
        return new(σ)
    end
end

"""
    draw(noise::MeasurementNoise, rng::AbstractRNG) -> Float64

Draw one Gaussian measurement residual.

# Arguments
- `noise::MeasurementNoise`: Noise model in the measurement's units.
- `rng::AbstractRNG`: Random-number generator used for the draw.

# Notes
Returns exactly `0.0` when `noise.sigma == 0.0`.

# Returns
A `Float64` residual in the measurement's units.

# Example
```jldoctest
using AstroSolve, Random
AstroSolve.draw(MeasurementNoise(0.0), MersenneTwister(1))

# output
0.0
```
"""
@inline function draw(noise::MeasurementNoise, rng::Random.AbstractRNG)
    return noise.sigma == 0.0 ? 0.0 : noise.sigma * randn(rng)
end

"""
    variance(noise::MeasurementNoise) -> Float64

Return the measurement variance used by an estimator.

# Arguments
- `noise::MeasurementNoise`: Scalar Gaussian noise model.

# Returns
The `Float64` value `noise.sigma^2`, in squared measurement units.

# Example
```jldoctest
using AstroSolve
AstroSolve.variance(MeasurementNoise(0.015))

# output
0.000225
```
"""
variance(noise::MeasurementNoise) = noise.sigma^2

# =============================================================================
# 3. Observable tags
# =============================================================================

"""
    AbstractMeasurement

Extension interface for tracking measurements used by an estimation problem.

# Notes
A concrete subtype carries its `SignalPath`, noise model, and any
measurement-specific parameters. Measurement prediction uses the trajectory
owned by the estimation problem.

# Example
```julia
using AstroSolve
station = (name = "DSS-14",)
sat = (name = "Sat",)
path = SignalPath(station, sat, station)
tracking = AbstractMeasurement[
    TwoWayRange(path; noise = MeasurementNoise(15.0e-3)),
    TwoWayDoppler(path; noise = MeasurementNoise(2.0e-5))]
```
"""
abstract type AbstractMeasurement end

"""
    TwoWayRange(path::SignalPath; noise = MeasurementNoise(0.0))

A geometric two-way range measurement in kilometres.

# Fields
- `path::P`: Signal route, where `P <: SignalPath`.
- `noise::N`: Noise model, where `N <: AbstractMeasurementNoise`.

# Notes
Prediction sums the geometric distance of every leg in `path`; it does not
implicitly double a one-way path. The canonical round trip is
`SignalPath(station, spacecraft, station)`. Noise is specified as a standard
deviation in kilometres and defaults to zero.

# Example
```julia
using AstroSolve
station = (name = "DSS-14",)
sat = (name = "Sat",)
TwoWayRange(SignalPath(station, sat, station); noise = MeasurementNoise(15.0e-3))
```
"""
struct TwoWayRange{P<:SignalPath, N<:AbstractMeasurementNoise} <: AbstractMeasurement
    path::P
    noise::N
end
TwoWayRange(path::SignalPath;
            noise::AbstractMeasurementNoise = MeasurementNoise(0.0)) =
    TwoWayRange{typeof(path), typeof(noise)}(path, noise)

"""
    TwoWayDoppler(path::SignalPath; noise = MeasurementNoise(0.0))

A geometric two-way Doppler measurement reported as range rate in kilometres per second.

# Fields
- `path::P`: Signal route, where `P <: SignalPath`.
- `noise::N`: Noise model, where `N <: AbstractMeasurementNoise`.

# Notes
Prediction sums the range rate of every leg in `path`. For
`SignalPath(station, spacecraft, station)`, the result is twice the one-way
range rate of the spacecraft relative to the station. The model is geometric
and does not apply light-time or relativistic corrections. Noise is specified
as a standard deviation in kilometres per second and defaults to zero.

# Example
```julia
using AstroSolve
station = (name = "DSS-14",)
sat = (name = "Sat",)
TwoWayDoppler(SignalPath(station, sat, station); noise = MeasurementNoise(2.0e-5))
```
"""
struct TwoWayDoppler{P<:SignalPath, N<:AbstractMeasurementNoise} <: AbstractMeasurement
    path::P
    noise::N
end
TwoWayDoppler(path::SignalPath;
              noise::AbstractMeasurementNoise = MeasurementNoise(0.0)) =
    TwoWayDoppler{typeof(path), typeof(noise)}(path, noise)

# =============================================================================
# 4. Visibility — open generic
# =============================================================================

"""
    is_visible(participant, r_sat, t) -> Bool

Return whether a signal-path participant can observe a spacecraft position at an epoch.

# Arguments
- `participant`: Signal-path participant.
- `r_sat`: Spacecraft position passed to the participant's visibility method.
- `t`: Observation epoch.

# Notes
The fallback method returns `true` for participants with no visibility
restriction. `GroundStation` uses the method defined by AstroModels, which
expects `r_sat` in kilometres in GCRF axes and applies the station's
minimum-elevation limit.

# Returns
`true` when the participant permits the observation.

# Example
```jldoctest
using AstroSolve, AstroEpochs
participant = (name = "Relay",)
epoch = Time("2020-03-01T00:00:00.000", TT(), ISOT())
AstroSolve.is_visible(participant, [7000.0, 0.0, 0.0], epoch)

# output
true
```
"""
is_visible(::Any, r_sat, t) = true

# =============================================================================
# 5. Exports
# =============================================================================

export SignalPath, participant_names,
       AbstractMeasurement, TwoWayRange, TwoWayDoppler,
       AbstractMeasurementNoise, MeasurementNoise,
       draw, variance, is_visible

end # module Measurements
