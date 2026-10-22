# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The public interface: ephemeris files read into, and written from, the trajectory a spacecraft
# carries, a SpacecraftHistory. What a file holds beyond times, states and a coordinate system
# (interpolation, accelerations, covariance, the object's name and ID, the usable span, and where
# the data came from) travels in each segment's `metadata`, under the keys in `_META`.

"""
    OEM()

The CCSDS Orbit Ephemeris Message format, `.oem`. Passed as `format =` when the extension does not
say which format a file is in.

# Examples
```jldoctest
OEM()

# output

OEM()
```
"""
struct OEM end

"""
    STKEphemeris()

The STK ephemeris format, `.e`. Passed as `format =` when the extension does not say which format
a file is in.

# Examples
```jldoctest
STKEphemeris()

# output

STKEphemeris()
```
"""
struct STKEphemeris end

"""
    SPK()

The JPL SPK format, `.bsp`. Passed as `format =` when the extension does not say which format a
file is in.

# Examples
```jldoctest
SPK()

# output

SPK()
```
"""
struct SPK end

# A time this close to a segment's first or last epoch is taken as that epoch. Files store epochs
# to finite precision (SPK to about 0.1 µs), so a history read back can end a few nanoseconds
# before the time it was written for; this keeps such a time evaluable without extrapolating.
const _EDGE = 1e-6                                                     # s

const _META = (interpolation = "interpolation", accelerations = "accelerations",
               covariances = "covariances", object_name = "object_name", object_id = "object_id",
               useable_start = "useable_start", useable_stop = "useable_stop", source = "source")

function _format_of(path::AbstractString, format)
    format === nothing || return format
    ext = lowercase(splitext(path)[2])
    ext == ".oem" && return OEM()
    ext == ".e" && return STKEphemeris()
    ext in (".bsp", ".spk") && return SPK()
    throw(ArgumentError("the extension $(repr(ext)) does not say which ephemeris format $(basename(path)) is; " *
                        "give format = OEM(), STKEphemeris() or SPK()"))
end

# ============================================================================================
# Reading
# ============================================================================================

"""
    read_ephemeris(path; format, kwargs...) -> SpacecraftHistory

Read an ephemeris file into a `SpacecraftHistory`: CCSDS OEM, STK `.e` or JPL SPK.

# Arguments
- `path`: The file. Its extension, `.oem`, `.e` or `.bsp`, chooses the format.
- `format`: [`OEM`](@ref)`()`, [`STKEphemeris`](@ref)`()` or [`SPK`](@ref)`()`, when the extension
  does not say.
- `strict`: Treat every correction as an error. Default `false`.
- `skip_bad_lines`: OEM and STK. Drop an unreadable data line, with a warning, and split the
  segment there. Default `false`.
- `ref_frame`, `center`: Names to read in place of the file's frame and centre, such as
  `"EME2000"` and `"EARTH"`. `time_system`: OEM only. Each override is recorded in the
  diagnostics with the file's original value.
- `naif_id`: SPK. The object to read, when the file holds more than one.
- `step`: SPK. Seconds between samples for SPK types that store none, such as Chebyshev.

# Notes
- Each segment of the file is a segment of the history, in the coordinate system the file names.
  Interpolation never crosses from one segment to the next.
- A file's problems are reported all at once. What can be corrected without changing the data is
  corrected and listed in one warning; the rest stops the read with an [`EphemerisReadError`](@ref)
  giving every problem with its line. [`diagnostics`](@ref) returns what was corrected or assumed.
- Beyond times and states, each segment's `metadata` holds the declared interpolation,
  accelerations and covariance where the file has them (km, km/s, km/s², km² and so on), the
  object's name and ID, the usable span, and the file's details.

# Returns
A `SpacecraftHistory`, one segment per file segment.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0, name = "SAT")
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 1800.0))
path = tempname() * ".oem"
write_ephemeris(path, sat; step = 60.0)
eph = read_ephemeris(path)
length(epochs(eph))

# output

31
```
"""
function read_ephemeris(path::AbstractString; format = nothing, naif_id = nothing, step = nothing,
                        strict::Bool = false, skip_bad_lines::Bool = false,
                        ref_frame = nothing, center = nothing, time_system = nothing)
    fmt = _format_of(path, format)
    data = if fmt isa OEM
        read_oem(path; strict, skip_bad_lines, ref_frame = _name(ref_frame), center = _name(center),
                 time_system = _name(time_system))
    elseif fmt isa STKEphemeris
        time_system === nothing || throw(ArgumentError("read_ephemeris: an STK file's times are UTC-based by definition; time_system is for OEM"))
        read_stk(path; strict, skip_bad_lines, ref_frame = _name(ref_frame), center = _name(center))
    else
        (ref_frame === nothing && center === nothing && time_system === nothing) ||
            throw(ArgumentError("read_ephemeris: an SPK file names its frame and centre by NAIF code; overrides are for OEM and STK"))
        read_spk(path; naif_id, step, strict)
    end
    return _to_history(data, path)
end

_name(x::Nothing) = nothing
_name(x::AbstractString) = String(x)
_name(x::Symbol) = String(x)

function _to_history(data::EphemerisFileData, path)
    source = (path = String(path), format = data.format, version = data.version, originator = data.originator,
              creation_date = data.creation_date, message_id = data.message_id, comments = data.comments,
              diagnostics = data.diagnostics)
    segs = HistorySegment[]
    for s in data.segments
        axes, _, problem = _frame_axes(s.frame, s.center)
        body, _ = _center_body(s.center)
        cs = CoordinateSystem(body, axes)
        md = Dict{String, Any}(
            _META.interpolation => s.interpolation, _META.object_name => s.object_name,
            _META.object_id => s.object_id, _META.source => source)
        s.accelerations === nothing || (md[_META.accelerations] = s.accelerations)
        isempty(s.covariances) || (md[_META.covariances] = s.covariances)
        s.useable_start === nothing || (md[_META.useable_start] = s.useable_start)
        s.useable_stop === nothing || (md[_META.useable_stop] = s.useable_stop)
        md["time_system"] = s.time_system
        md["frame_name"] = s.frame
        states = [CartesianState(Vector(x)) for x in s.states]
        push!(segs, HistorySegment(s.times, states, cs; name = s.object_name, metadata = md))
    end
    return SpacecraftHistory(segs)
end

# ============================================================================================
# Inspecting
# ============================================================================================

_segments(h::SpacecraftHistory) = filter(!isempty, h.segments)
_segments(sat::Spacecraft) = _segments(sat.history)

"""
    epochs(eph) -> Vector{Time}

Every epoch in a history or a spacecraft's history, in segment order. An epoch shared by two
segments, such as a maneuver time, appears once for each.

# Arguments
- `eph`: A `SpacecraftHistory`, or a `Spacecraft`.

# Returns
A vector of `Time`.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
path = tempname() * ".oem"
write_ephemeris(path, sat; step = 60.0)
length(epochs(read_ephemeris(path)))

# output

11
```
"""
epochs(eph::Union{SpacecraftHistory, Spacecraft}) = reduce(vcat, (s.times for s in _segments(eph)); init = Time{Float64}[])

"""
    span(eph) -> Tuple{Time, Time}

The first and last epochs of a history. A history may have gaps inside its span; see
[`coverage`](@ref).

# Arguments
- `eph`: A `SpacecraftHistory`, or a `Spacecraft`.

# Returns
A tuple `(first, last)`.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
t_start, t_stop = span(sat)
round((t_stop - t_start) * 86400; digits = 6)

# output

600.0
```
"""
function span(eph::Union{SpacecraftHistory, Spacecraft})
    c = coverage(eph)
    isempty(c) && throw(ArgumentError("span: the history is empty"))
    return first(c[1]), last(c[end])
end

"""
    coverage(eph) -> Vector{Tuple{Time, Time}}

The intervals a history covers, in time order: segments that touch or overlap are merged, so a
gap between segments shows as a gap between intervals. A segment's usable span, where its file
gives one, limits what it covers.

# Arguments
- `eph`: A `SpacecraftHistory`, or a `Spacecraft`.

# Returns
A vector of `(start, stop)` tuples.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
length(coverage(sat))

# output

1
```
"""
function coverage(eph::Union{SpacecraftHistory, Spacecraft})
    iv = sort([_usable(s) for s in _segments(eph)]; by = x -> _seconds(first(_usable(_segments(eph)[1])), x[1]))
    out = Tuple{Time, Time}[]
    for (a, b) in iv
        if !isempty(out) && _seconds(out[end][2], a) <= 0
            _seconds(out[end][2], b) > 0 && (out[end] = (out[end][1], b))
        else
            push!(out, (a, b))
        end
    end
    return out
end

function _usable(s::HistorySegment)
    a = get(s.metadata, _META.useable_start, first(s.times))
    b = get(s.metadata, _META.useable_stop, last(s.times))
    return (a, b)
end

"""
    diagnostics(eph) -> Vector{Diagnostic}

What the reader corrected or assumed when it read the file a history came from, as
[`Diagnostic`](@ref)s. Empty for a history that was not read from a file, or for a file read
without corrections.

# Arguments
- `eph`: A `SpacecraftHistory` from [`read_ephemeris`](@ref).

# Returns
A vector of `Diagnostic`, each with its line, the field, and what was done.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
path = tempname() * ".bsp"
write_ephemeris(path, sat; step = 60.0, naif_id = -900001)
first(diagnostics(read_ephemeris(path))).message

# output

"SPICE's J2000 frame read as ICRF, as NAIF defines it"
```
"""
function diagnostics(eph::SpacecraftHistory)
    for s in _segments(eph)
        haskey(s.metadata, _META.source) && return s.metadata[_META.source].diagnostics
    end
    return Diagnostic[]
end

# ============================================================================================
# Evaluating
# ============================================================================================

"""
    state_at(eph, t; overlap = :error) -> Coordinate

The interpolated state at `t`, with its coordinate system and epoch.

# Arguments
- `eph`: A `SpacecraftHistory`, or a `Spacecraft`.
- `t`: The epoch, in any time scale.
- `overlap`: Where two segments' usable spans overlap: `:error` (default), `:earlier` or `:later`.

# Notes
- The segment that covers `t` is interpolated with its declared method, Lagrange degree 7 when
  none is declared, and never with points from another segment.
- At an epoch two segments share, such as a maneuver, the later segment applies: the state after
  the event.
- `t` outside every segment, inside a gap, is an `ArgumentError` naming the nearest covered times.
  A time within 1 µs of a segment's first or last epoch is taken as that epoch, since files store
  epochs to finite precision.
- Hermite interpolation uses the tabulated velocities. Lagrange interpolates positions and
  velocities independently, as OEM and STK readers do; where a file has positions only, the
  velocity is the derivative of the position polynomial.

# Returns
A `Coordinate` in the segment's coordinate system.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
c = state_at(sat, Time("2026-10-21T00:05:00", UTC(), ISOT()))
round(to_vector(CartesianState(c))[1]; digits = 3)

# output

6637.168
```
"""
function state_at(eph::Union{SpacecraftHistory, Spacecraft}, t::Time; overlap::Symbol = :error)
    overlap in (:error, :earlier, :later) || throw(ArgumentError("state_at: overlap must be :error, :earlier or :later"))
    segs = _segments(eph)
    isempty(segs) && throw(ArgumentError("state_at: the history is empty"))
    covering = Int[]
    for (k, s) in enumerate(segs)
        a, b = _usable(s)
        _seconds(a, t) >= -_EDGE && _seconds(t, b) >= -_EDGE && push!(covering, k)
    end
    if isempty(covering)
        c = coverage(eph)
        before = filter(iv -> _seconds(iv[2], t) > 0, c)
        after = filter(iv -> _seconds(t, iv[1]) > 0, c)
        near = String[]
        isempty(before) || push!(near, "the history before it ends at $(_label(last(before)[2]))")
        isempty(after) || push!(near, "the history after it starts at $(_label(first(after)[1]))")
        throw(ArgumentError("state_at: $(_label(t)) is not covered" *
                            (isempty(before) || isempty(after) ? "; there is no extrapolation" : ", it is in a gap") *
                            ". " * join(near, "; ")))
    end
    k = if length(covering) == 1
        only(covering)
    else
        # Two segments touching at t: the later one, the state after the event.
        a, b = segs[covering[1]], segs[covering[end]]
        shared = length(covering) == 2 && abs(_seconds(t, _usable(a)[2])) <= _EDGE && abs(_seconds(_usable(b)[1], t)) <= _EDGE
        if shared || overlap == :later
            covering[end]
        elseif overlap == :earlier
            covering[1]
        else
            throw(ArgumentError("state_at: segments $(join(covering, " and ")) both cover $(_label(t)); " *
                                "choose one with overlap = :earlier or :later"))
        end
    end
    return _interpolate_segment(segs[k], t)
end

_label(t::Time) = _format_epoch(t, uppercase(String(t.scale)) in keys(_TIME_SYSTEMS) ? uppercase(String(t.scale)) : "UTC", 3) *
                  " " * uppercase(String(t.scale))

function _interpolate_segment(s::HistorySegment, t::Time)
    # One state, as a maneuver records: exact at its epoch, which is all it covers.
    length(s.times) == 1 && return Coordinate(s.states[1], s.coordinate_system, t)
    t0 = first(s.times)
    x = [_seconds(t0, ti) for ti in s.times]
    xt = clamp(_seconds(t0, t), first(x), last(x))             # within _EDGE of an end: that end
    n = length(x)
    m = get(s.metadata, _META.interpolation, nothing)
    m === nothing && (m = Lagrange(min(7, n - 1)))
    n >= required_points(m) || (m = m isa Hermite ? Hermite(2n - 1) : Lagrange(n - 1))
    vals = [to_vector(st) for st in s.states]
    r = [SVector{3}(v[1:3]) for v in vals]
    v = [SVector{3}(v[4:6]) for v in vals]
    if all(w -> all(isnan, w), v)                              # positions only
        pos, vel = interpolate_with_derivative(m isa Hermite ? Lagrange(min(7, n - 1)) : m, x, r, xt)
    elseif m isa Hermite
        pos, vel = interpolate_with_derivative(m, x, r, v, xt)
    else
        pos, vel = interpolate(m, x, r, xt), interpolate(m, x, v, xt)
    end
    return Coordinate(CartesianState(vcat(pos, vel)), s.coordinate_system, t)
end

"""
    position_difference(a, b) -> Vector{Float64}

The distance between two trajectories at each epoch of `a` that `b` covers, in km.

# Arguments
- `a`: The reference, a `SpacecraftHistory` or a `Spacecraft`.
- `b`: The trajectory compared with it, a `SpacecraftHistory` or a `Spacecraft`.

# Notes
- `b` is interpolated at `a`'s epochs and expressed in `a`'s coordinate system before the
  difference is taken. Epochs of `a` that `b` does not cover are skipped.
- At an epoch two segments of `b` share, the later segment applies, as in [`state_at`](@ref).

# Returns
A vector of distances in km, one per epoch of `a` covered by `b`.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 1800.0))
path = tempname() * ".oem"
write_ephemeris(path, sat; step = 60.0)
maximum(position_difference(read_ephemeris(path), sat)) < 1e-6

# output

true
```
"""
function position_difference(a::Union{SpacecraftHistory, Spacecraft}, b::Union{SpacecraftHistory, Spacecraft})
    out = Float64[]
    cov = coverage(b)
    for s in _segments(a), (t, st) in zip(s.times, s.states)
        any(iv -> _seconds(iv[1], t) >= -_EDGE && _seconds(t, iv[2]) >= -_EDGE, cov) || continue
        cb = state_at(b, t; overlap = :later)
        xb = cb.coord_sys == s.coordinate_system ? to_vector(CartesianState(cb)) :
             to_vector(CartesianState(cb, s.coordinate_system))
        xa = to_vector(st)
        push!(out, sqrt(sum(abs2, xa[1:3] - xb[1:3])))
    end
    return out
end

# ============================================================================================
# Writing
# ============================================================================================

"""
    write_ephemeris(path, source; format, coordinate_system, step, interpolation, time_system,
                    naif_id, object_name, object_id, originator, version, digits,
                    discard_unsupported) -> String

Write a spacecraft's trajectory, or any `SpacecraftHistory`, to an ephemeris file: CCSDS OEM,
STK `.e` or JPL SPK.

# Arguments
- `path`: The file to write. Its extension, `.oem`, `.e` or `.bsp`, chooses the format; it is
  replaced if it exists.
- `source`: A `Spacecraft`, whose history is written, or a `SpacecraftHistory`, such as one from
  [`read_ephemeris`](@ref).
- `format`: [`OEM`](@ref)`()`, [`STKEphemeris`](@ref)`()` or [`SPK`](@ref)`()`, when the extension
  does not say.
- `coordinate_system`: Write in this coordinate system, converting every state. Default: each
  segment's own, except for SPK, which is written in ICRF axes about the segment's origin.
- `step`: Resample every segment every `step` seconds, each segment's first and last epochs
  included, by interpolating it. Default `nothing`: the epochs as they are.
- `interpolation`: The method the file declares, such as `Hermite(7)`. Default: what the source
  declares, else `Lagrange(7)`.
- `time_system`: OEM. The time system epochs are written in, e.g. `UTC()` (default) or `TDB()`, or
  `"GPS"`.
- `naif_id`: SPK. The object's NAIF ID; required for SPK.
- `object_name`, `object_id`: The object. Default: the source's, else the spacecraft's name, else
  `"UNKNOWN"`.
- `originator`: OEM. Who made the file. Default `"EPICYCLE"`.
- `version`: OEM. `"3.0"` (default) or `"2.0"`.
- `digits`: OEM and STK. Decimal places of seconds in written epochs. Default 9.
- `discard_unsupported`: SPK. Write without accelerations and covariance, which SPK cannot hold,
  rather than refuse. Default `false`.

# Notes
- Every rule of the format is checked before the file is opened, and every problem is reported
  in one `ArgumentError`. Nothing is written on failure, and an existing file is left as it was.
- Each segment of the history is a segment of the file, so a maneuver stays a discontinuity. The
  one-state segment `maneuver!` records repeats the first state of the propagation after it and
  is not written.
- What the source carries beyond states (accelerations, covariance) is written where the format
  can hold it. A conversion to another coordinate system rotates states only: a source with
  accelerations or covariance in another frame is refused rather than written half-converted.

# Returns
`path`.

# Examples
```jldoctest
using Epicycle
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
sat = Spacecraft(state = CartesianState([7000.0, 0, 0, 0, 7.5, 1.0]), time = t0)
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())), IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12))
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 600.0))
path = tempname() * ".e"
write_ephemeris(path, sat; step = 60.0, interpolation = Hermite(7))
filter(startswith("InterpolationSamplesM1"), readlines(path))

# output

1-element Vector{String}:
 "InterpolationSamplesM1  3"
```
"""
function write_ephemeris(path::AbstractString, source::Union{SpacecraftHistory, Spacecraft};
                         format = nothing, coordinate_system = nothing, step = nothing, interpolation = nothing,
                         time_system = UTC(), naif_id = nothing, object_name = nothing, object_id = nothing,
                         originator::AbstractString = "EPICYCLE", version::AbstractString = "3.0",
                         digits::Integer = 9, discard_unsupported::Bool = false)
    fmt = _format_of(path, format)
    segs = _writable(_segments(source))
    isempty(segs) && throw(ArgumentError("write_ephemeris: nothing to write; the history is empty" *
        (source isa Spacecraft ? " (was the spacecraft propagated?)" : "")))
    default_name = source isa Spacecraft ? source.name : ""
    sys = fmt isa SPK ? "TDB" : fmt isa STKEphemeris ? "UTC" : _time_system_name(time_system)
    out = EphemerisSegmentData[]
    for s in segs
        target = coordinate_system !== nothing ? coordinate_system :
                 fmt isa SPK && !(typeof(s.coordinate_system.axes) in (ICRF, MJ2000Ec)) ?
                     CoordinateSystem(s.coordinate_system.origin, ICRF()) : s.coordinate_system
        push!(out, _to_segment(s, target, step, interpolation, sys, object_name, object_id, default_name))
    end
    data = EphemerisFileData(format = fmt isa OEM ? :oem : fmt isa STKEphemeris ? :stk : :spk,
                             originator = originator, segments = out)
    if fmt isa OEM
        write_oem(path, data; digits, version)
    elseif fmt isa STKEphemeris
        write_stk(path, data; digits)
    else
        naif_id === nothing && throw(ArgumentError("write_ephemeris: SPK needs the object's NAIF ID; give naif_id = …, " *
                                                  "negative for a spacecraft"))
        write_spk(path, data; naif_id, discard_unsupported)
    end
    return path
end

# A segment of one state, as a maneuver records, repeats the first state of the propagation after
# it: dropped, since nothing is lost. Any other one-state segment cannot be interpolated.
function _writable(segs)
    out = HistorySegment[]
    for (k, s) in enumerate(segs)
        if length(s.times) == 1
            # The same epoch and state, compared in the neighbour's coordinate system: a maneuver
            # records in the spacecraft's, a propagation in its own.
            function same(o, j)
                abs(_seconds(s.times[1], o.times[j])) <= _EDGE || return false
                x = s.coordinate_system == o.coordinate_system ? to_vector(s.states[1]) :
                    to_vector(CartesianState(Coordinate(s.states[1], s.coordinate_system, s.times[1]), o.coordinate_system))
                y = to_vector(o.states[j])
                return maximum(abs.(x[1:3] - y[1:3])) <= 1e-9 * max(1.0, maximum(abs, y[1:3])) &&
                       maximum(abs.(x[4:6] - y[4:6])) <= 1e-9 * max(1.0, maximum(abs, y[4:6]))
            end
            (k < length(segs) && same(segs[k + 1], 1)) || (k > 1 && same(segs[k - 1], length(segs[k - 1].times))) ||
                throw(ArgumentError("write_ephemeris: segment $k$(isempty(s.name) ? "" : " ($(s.name))") holds one state, " *
                                    "which no format can interpolate"))
            continue
        end
        push!(out, s)
    end
    return out
end

function _to_segment(s::HistorySegment, target, step, interpolation, sys, object_name, object_id, default_name)
    md = s.metadata
    cs = s.coordinate_system
    times = s.times
    states = [SVector{6, Float64}(to_vector(x)) for x in s.states]
    acc = get(md, _META.accelerations, nothing)
    covs = get(md, _META.covariances, CovarianceRecord[])
    m = interpolation !== nothing ? interpolation : get(md, _META.interpolation, nothing)
    if step !== nothing
        step > 0 || throw(ArgumentError("write_ephemeris: step must be positive; got $step"))
        (acc === nothing && isempty(covs)) || throw(ArgumentError("write_ephemeris: resampling with step " *
            "would drop the accelerations and covariance this source carries; write it without step"))
        t0 = first(times)
        total = _seconds(t0, last(times))
        offs = collect(0.0:Float64(step):total)
        total - last(offs) > 1e-9 && push!(offs, total)
        newtimes = [Time(t0.jd1, t0.jd2 + o / 86400, t0.scale, :jd) for o in offs]
        newtimes[end] = last(times)
        states = [SVector{6, Float64}(to_vector(CartesianState(_interpolate_segment(s, t)))) for t in newtimes]
        times = newtimes
    end
    if target != cs
        (acc === nothing && isempty(covs)) || throw(ArgumentError("write_ephemeris: converting to another " *
            "coordinate system rotates states only; this source carries accelerations or covariance in its own frame"))
        states = [SVector{6, Float64}(to_vector(CartesianState(Coordinate(CartesianState(Vector(x)), cs, t), target)))
                  for (x, t) in zip(states, times)]
    end
    center = _center_name(target.origin)
    frame = _frame_name(target.axes, center)
    name = something(object_name, get(md, _META.object_name, nothing), isempty(default_name) ? nothing : default_name, "UNKNOWN")
    id = something(object_id, get(md, _META.object_id, nothing), "UNKNOWN")
    return EphemerisSegmentData(object_name = name, object_id = id, center = center, frame = frame,
        time_system = sys, interpolation = something(m, Lagrange(7)), times = collect(times), states = states,
        accelerations = acc, covariances = covs,
        useable_start = step === nothing ? get(md, _META.useable_start, nothing) : nothing,
        useable_stop = step === nothing ? get(md, _META.useable_stop, nothing) : nothing)
end
