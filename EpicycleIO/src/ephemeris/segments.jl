# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# What a file holds, in the file's own terms: names as written, units km and km/s, epochs in the
# file's time system. The readers produce these and the writers consume them; turning them into
# Epicycle frames and trajectories is a separate step.

"""
    CovarianceRecord

One 6×6 position–velocity covariance from an ephemeris file.

# Fields
- `epoch`: The epoch of the covariance.
- `frame`: The reference frame name, as an ephemeris file gives it (e.g. `"EME2000"`).
- `matrix`: The symmetric covariance: km² for position–position, km²/s for position–velocity,
  km²/s² for velocity–velocity, ordered x, y, z, vx, vy, vz.

# Examples
```jldoctest
using AstroEpochs, StaticArrays
c = EpicycleIO.CovarianceRecord(Time("2026-01-01T00:00:00", UTC(), ISOT()), "EME2000",
                                SMatrix{6,6}(zeros(6, 6)))
c.frame

# output

"EME2000"
```
"""
struct CovarianceRecord
    epoch::Time
    frame::String
    matrix::SMatrix{6, 6, Float64, 36}
end

"""
    EphemerisSegmentData

One segment of an ephemeris file: metadata and the states that follow it. Interpolation never
crosses from one segment to the next.

# Fields
- `object_name`, `object_id`: The object the states describe.
- `center`: The origin name, upper case (e.g. `"EARTH"`).
- `frame`: The reference frame name (e.g. `"EME2000"`).
- `frame_epoch`: The frame's epoch where the frame needs one, else `nothing`.
- `time_system`: The time system name, upper case (e.g. `"UTC"`).
- `start_time`, `stop_time`: The span the segment declares.
- `useable_start`, `useable_stop`: The span recommended for use, or `nothing`.
- `interpolation`: The declared method, or `nothing` when the file declares none.
- `times`: Epochs, strictly increasing.
- `states`: Position and velocity at each epoch, km and km/s.
- `accelerations`: Acceleration at each epoch, km/s², or `nothing`.
- `covariances`: Covariances in this segment, by increasing epoch.
- `comments`: Comment text, without the `COMMENT` keyword.

# Examples
```jldoctest
fieldcount(EpicycleIO.EphemerisSegmentData)

# output

16
```
"""
Base.@kwdef mutable struct EphemerisSegmentData
    object_name::String = "UNKNOWN"
    object_id::String = "UNKNOWN"
    center::String = ""
    frame::String = ""
    frame_epoch::Union{Nothing, Time} = nothing
    time_system::String = ""
    start_time::Union{Nothing, Time} = nothing
    stop_time::Union{Nothing, Time} = nothing
    useable_start::Union{Nothing, Time} = nothing
    useable_stop::Union{Nothing, Time} = nothing
    interpolation::Union{Nothing, AbstractInterpolationMethod} = nothing
    times::Vector{Time} = Time[]
    states::Vector{SVector{6, Float64}} = SVector{6, Float64}[]
    accelerations::Union{Nothing, Vector{SVector{3, Float64}}} = nothing
    covariances::Vector{CovarianceRecord} = CovarianceRecord[]
    comments::Vector{String} = String[]
end

"""
    EphemerisFileData

The contents of an ephemeris file: its header and its segments, plus anything the reader
corrected.

# Fields
- `format`: `:oem`, `:stk` or `:spk`.
- `version`: The format version, e.g. `"3.0"`.
- `creation_date`: As written in the file, or empty.
- `originator`: As written in the file, or empty.
- `message_id`: OEM 3.0 only; empty when absent.
- `classification`: OEM 3.0 only; empty when absent.
- `comments`: Header comments.
- `segments`: The [`EphemerisSegmentData`](@ref) in file order.
- `diagnostics`: Problems the reader corrected, as [`Diagnostic`](@ref)s.

# Examples
```jldoctest
EpicycleIO.EphemerisFileData(format = :oem).version

# output

"3.0"
```
"""
Base.@kwdef mutable struct EphemerisFileData
    format::Symbol = :oem
    version::String = "3.0"
    creation_date::String = ""
    originator::String = ""
    message_id::String = ""
    classification::String = ""
    comments::Vector{String} = String[]
    segments::Vector{EphemerisSegmentData} = EphemerisSegmentData[]
    diagnostics::Vector{Diagnostic} = Diagnostic[]
end

# Seconds from a to b, measured in a's time scale; b may be in another. One scale throughout
# matters: elapsed TDB and elapsed UTC differ by about 1 µs an hour, so a table of epochs and the
# time it is interpolated at must be measured the same way.
_seconds(a::Time, b::Time) = (getproperty(b, a.scale) - a) * 86400.0
