# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# JPL SPK files (.bsp), through CSPICE (SPICE.jl). Writing uses SPK type 9 (Lagrange) or type 13
# (Hermite), one SPK segment per Epicycle segment. Reading takes the stored samples of types 9
# and 13 directly, and samples any other type through SPICE at a step the caller gives.
#
# Layout of a type 9 or 13 segment, as SPICE writes it: N states (6N values), N epochs, an epoch
# directory, then two values: for type 9 the polynomial degree and N; for type 13 the window
# size less one and N, the Hermite degree being 2·(window size) − 1.

# CSPICE keeps global state and is not safe to call from several threads at once.
const _SPICE_LOCK = ReentrantLock()

const _SPK_FRAMES = Dict(1 => "ICRF", 17 => "J2000_ECLIPTIC")      # SPICE frame code → name
const _SPK_FRAME_NAMES = Dict("ICRF" => "J2000", "J2000_ECLIPTIC" => "ECLIPJ2000")
const _SPK_MAX_DEGREE = 27                                         # types 9 and 13, CSPICE MAXDEG
const _SPK_SEGID_MAX = 40
const _J2000_TDB_JD = 2451545.0

# Ephemeris time ↔ Time. Whole days and the fraction are kept apart, so the conversion adds no
# error to what a 64-bit ephemeris time can hold (about 0.1 µs this century).
function _et(t::Time)
    u = t.tdb
    days = (u.jd1 - _J2000_TDB_JD) + floor(u.jd2)
    return days * 86400.0 + (u.jd2 - floor(u.jd2)) * 86400.0
end
function _from_et(et::Real)
    days = floor(et / 86400.0)
    return Time(_J2000_TDB_JD + days, (et - days * 86400.0) / 86400.0, :tdb, :jd)
end

function _naif_center(code::Integer)
    for (name, body) in _BODIES
        body.naifid == code && return name
    end
    return nothing
end

# ============================================================================================
# Reading
# ============================================================================================

"""
    read_spk(path; naif_id = nothing, step = nothing, strict = false) -> EphemerisFileData

Read the trajectory of one object from a JPL SPK file.

# Arguments
- `path`: The SPK file.
- `naif_id`: The NAIF ID of the object to read. Needed when the file holds more than one.
- `step`: Seconds between samples for SPK types that store no samples, such as Chebyshev
  (types 2, 3, 21). Not used for types 9 and 13.
- `strict`: Treat every correction as an error.

# Notes
- Types 9 (Lagrange) and 13 (Hermite) are read as stored: the samples, epochs and degree, with
  no interpolation. A segment's coverage, where narrower than its samples, becomes its usable span.
- Other types are evaluated by SPICE at every `step` seconds. The reader then compares Lagrange
  degree 7 between the samples with SPICE at the midpoints and records the largest difference
  as a note, so the caller can judge the step.
- SPICE's J2000 frame is read as ICRF, as NAIF defines it; ECLIPJ2000 as J2000_ECLIPTIC.
- Epochs are SPICE ephemeris time (TDB seconds past J2000) and keep its 64-bit precision, about
  0.1 µs in this century.

# Returns
An [`EphemerisFileData`](@ref) with `format = :spk`, one segment per SPK segment for the object.

# Examples
```jldoctest
using AstroEpochs, StaticArrays
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
times = [Time(t0.jd1, t0.jd2 + k / 1440, :utc, :jd) for k in 0:9]
states = [SVector(7000.0, 450.0k, 0.0, 0.0, 7.5, 0.0) for k in 0:9]
seg = EpicycleIO.EphemerisSegmentData(center = "EARTH", frame = "ICRF", time_system = "UTC",
          times = times, states = states, interpolation = EpicycleIO.Hermite(7))
path = tempname() * ".bsp"
EpicycleIO.write_spk(path, EpicycleIO.EphemerisFileData(segments = [seg]); naif_id = -900001)
f = EpicycleIO.read_spk(path)
(f.segments[1].interpolation, length(f.segments[1].times))

# output

(AstroRoutines.Hermite(7), 10)
```
"""
function read_spk(path::AbstractString; naif_id = nothing, step = nothing, strict::Bool = false)
    isfile(path) || throw(ArgumentError("read_spk: no file at $path"))
    d = _Diagnostics(path; strict)
    file = EphemerisFileData(format = :spk, version = "")
    segs = lock(_SPICE_LOCK) do
        _spk_summaries(d, path)
    end
    if !_has_errors(d)
        targets = unique(s.ic[1] for s in segs)
        if isempty(targets)
            _error!(d, 0, "file", "the file holds no SPK segments")
        elseif naif_id === nothing && length(targets) > 1
            _error!(d, 0, "naif_id", "the file holds objects $(join(targets, ", ")); give naif_id to choose one")
        else
            target = naif_id === nothing ? only(targets) : Int(naif_id)
            target in targets || _error!(d, 0, "naif_id", "object $target is not in the file, which holds $(join(targets, ", "))")
            for (k, s) in enumerate(segs)
                s.ic[1] == target || continue
                seg = _spk_segment!(d, path, k, s, step)
                seg === nothing || push!(file.segments, seg)
            end
        end
    end
    file.diagnostics = _finish(d)
    return file
end

# Every segment's descriptor and data, read under the SPICE lock.
function _spk_summaries(d, path)
    out = NamedTuple[]
    handle = nothing
    try
        handle = SPICE.dafopr(path)
        SPICE.dafbfs(handle)
        while SPICE.daffna()
            dc, ic = SPICE.dafus(SPICE.dafgs(), 2, 6)
            name = String(strip(SPICE.dafgn()))
            typ = ic[4]
            data = nothing
            if typ in (9, 13)
                tail = SPICE.dafgda(handle, ic[6] - 1, ic[6])
                n = round(Int, tail[2])
                raw = SPICE.dafgda(handle, ic[5], ic[5] + 7n - 1)
                data = (param = round(Int, tail[1]), n = n, raw = raw)
            end
            push!(out, (dc = collect(dc), ic = collect(ic), name = name, data = data))
        end
    catch e
        _error!(d, 0, "file", "SPICE could not read it as an SPK file: $(sprint(showerror, e))")
    finally
        handle === nothing || SPICE.dafcls(handle)
    end
    return out
end

function _spk_segment!(d, path, k, s, step)
    target, centre_code, frame_code, typ = s.ic[1], s.ic[2], s.ic[3], s.ic[4]
    where_ = "SPK segment $k ($(isempty(s.name) ? "unnamed" : repr(s.name)))"
    center = _naif_center(centre_code)
    center === nothing && return (_error!(d, 0, "center", "$where_ is relative to NAIF body $centre_code, which is not an Epicycle centre"); nothing)
    frame = get(_SPK_FRAMES, frame_code, nothing)
    seg = EphemerisSegmentData(object_name = isempty(s.name) ? "UNKNOWN" : s.name, object_id = string(target),
                               center = center, time_system = "TDB")
    if s.data !== nothing
        frame === nothing && return (_error!(d, 0, "frame", "$where_ is in SPICE frame $frame_code; Epicycle reads J2000 and ECLIPJ2000 segments"); nothing)
        frame == "ICRF" && _note!(d, 0, "frame", "SPICE's J2000 frame read as ICRF, as NAIF defines it")
        n, raw = s.data.n, s.data.raw
        seg.frame = frame
        seg.states = [SVector{6, Float64}(raw[6j - 5:6j]) for j in 1:n]
        seg.times = [_from_et(raw[6n + j]) for j in 1:n]
        seg.interpolation = typ == 9 ? Lagrange(s.data.param) : Hermite(2 * (s.data.param + 1) - 1)
        first_et, last_et = raw[6n + 1], raw[7n]
        s.dc[1] > first_et && (seg.useable_start = _from_et(s.dc[1]))
        s.dc[2] < last_et && (seg.useable_stop = _from_et(s.dc[2]))
        (s.dc[1] < first_et || s.dc[2] > last_et) && _fixed!(d, 0, "coverage",
            "$where_ declares coverage beyond its samples; limited to the samples")
    else
        step === nothing && return (_error!(d, 0, "step", "$where_ is SPK type $typ, which stores no samples; " *
                                            "read with step = seconds to sample it through SPICE"); nothing)
        step > 0 || return (_error!(d, 0, "step", "step must be positive; got $step"); nothing)
        _spk_sample!(d, path, s, seg, Float64(step), where_) || return nothing
    end
    seg.start_time, seg.stop_time = first(seg.times), last(seg.times)
    return seg
end

# Evaluate a segment that stores no samples, through SPICE, and measure the sampling error.
function _spk_sample!(d, path, s, seg, step, where_)
    a, b = s.dc
    ets = collect(a:step:b)
    # The end, unless the last sample is already there to within rounding of the coverage.
    b - last(ets) > 1e-6 * step ? push!(ets, b) : (ets[end] = b)
    length(ets) >= 8 || (_error!(d, 0, "step", "$where_ spans $(b - a) s; a step of $step s gives $(length(ets)) samples and at least 8 are needed"); return false)
    target, centre = s.ic[1], s.ic[2]
    sample(et) = SPICE.spkgeo(target, et, "J2000", centre)[1]
    states, mids = Vector{Float64}[], Float64[]
    ok = lock(_SPICE_LOCK) do
        try
            SPICE.furnsh(path)
            for et in ets
                push!(states, collect(sample(et)))
            end
            mids = [sample((ets[j] + ets[j + 1]) / 2)[1:3] for j in 1:length(ets) - 1]
            true
        catch e
            _error!(d, 0, "file", "SPICE could not evaluate $where_: $(sprint(showerror, e))")
            false
        finally
            try SPICE.unload(path) catch end
        end
    end
    ok || return false
    seg.frame = "ICRF"
    seg.states = [SVector{6, Float64}(x) for x in states]
    seg.times = _from_et.(ets)
    seg.interpolation = Lagrange(7)
    r = [x[1:3] for x in seg.states]
    err = maximum(sqrt(sum(abs2, interpolate(Lagrange(7), ets, r, (ets[j] + ets[j + 1]) / 2) - mids[j])) for j in 1:length(mids))
    _note!(d, 0, "step", "$where_ (SPK type $(s.ic[4])) sampled every $step s through SPICE, in its J2000 frame read " *
                         "as ICRF; Lagrange 7 between samples is within $(@sprintf("%.3g", err)) km of SPICE")
    return true
end

# ============================================================================================
# Writing
# ============================================================================================

"""
    write_spk(path, data::EphemerisFileData; naif_id, discard_unsupported = false)

Write `data` as a JPL SPK file: one type 13 segment for each Hermite segment, and one type 9
segment for each Lagrange segment.

# Arguments
- `path`: The file to write. It is replaced if it exists.
- `data`: The segments to write.
- `naif_id`: The NAIF ID of the object, an integer. Spacecraft conventionally use negative IDs.
- `discard_unsupported`: SPK stores no accelerations or covariance; with `false`, the default,
  data that carry them are refused rather than written without them.

# Notes
- Everything is checked before the file is opened, and every problem is reported in one
  `ArgumentError`; nothing is written on failure.
- Segments must be in ICRF (written as SPICE's J2000 frame, as NAIF defines it) or
  J2000_ECLIPTIC (ECLIPJ2000). Data in EME2000 must be rotated to ICRF first: SPICE's J2000 is
  ICRF, and EME2000 differs from it by the 23 mas frame bias.
- A segment's usable span, where given, becomes its SPK coverage; otherwise its first to last
  epoch. Usable spans of consecutive segments may not overlap, as SPICE would take the later
  segment wherever two cover a time.
- Degree at most 27, CSPICE's limit; Hermite degree odd. A segment without an interpolation
  is written as Lagrange degree 7.
- Epochs are converted to SPICE ephemeris time by Epicycle, so no leapseconds kernel is needed.

# Examples
See [`read_spk`](@ref).
"""
function write_spk(path::AbstractString, data::EphemerisFileData; naif_id, discard_unsupported::Bool = false)
    problems = String[]
    naif_id isa Integer || push!(problems, "naif_id must be an integer; got $(repr(naif_id))")
    naif_id isa Integer && _naif_center(naif_id) !== nothing &&
        push!(problems, "naif_id $naif_id is a natural body's; use the object's own ID (negative for a spacecraft)")
    isempty(data.segments) && push!(problems, "there are no segments to write")
    for (k, s) in enumerate(data.segments)
        what = "segment $k"
        body, prob = _center_body(s.center)
        body === nothing && push!(problems, "$what CENTER: $prob")
        haskey(_SPK_FRAME_NAMES, s.frame) || push!(problems, "$what is in $(s.frame); SPK segments are written in ICRF " *
            "or J2000_ECLIPTIC" * (s.frame == "EME2000" ? ". Rotate EME2000 to ICRF first: SPICE's J2000 is ICRF" : ""))
        n = length(s.times)
        n == length(s.states) || push!(problems, "$what has $n epochs but $(length(s.states)) states")
        n == 0 && (push!(problems, "$what has no states"); continue)
        any(j -> _seconds(s.times[j - 1], s.times[j]) <= 0, 2:n) && push!(problems, "$what epochs do not increase strictly")
        any(x -> !all(isfinite, x), s.states) && push!(problems, "$what has a state that is not finite " *
                                                        "(SPK stores velocities; a file read with positions only cannot be written)")
        m = something(s.interpolation, Lagrange(7))
        m.degree <= _SPK_MAX_DEGREE || push!(problems, "$what degree $(m.degree) is over $_SPK_MAX_DEGREE, CSPICE's limit")
        n >= required_points(m) || push!(problems, "$what has $n states; $m needs $(required_points(m))")
        if !discard_unsupported
            s.accelerations === nothing || push!(problems, "$what has accelerations, which SPK cannot store; " *
                                                  "write with discard_unsupported = true to drop them")
            isempty(s.covariances) || push!(problems, "$what has covariance, which SPK cannot store; write OEM or STK " *
                                                "to keep it, or write with discard_unsupported = true")
        end
        name = s.object_name
        (length(name) <= _SPK_SEGID_MAX && all(c -> isascii(c) && isprint(c), name)) ||
            push!(problems, "$what object name is used as the SPK segment ID, which is at most $_SPK_SEGID_MAX printable ASCII characters")
    end
    segs = data.segments
    for k in 2:length(segs)
        (isempty(segs[k - 1].times) || isempty(segs[k].times)) && continue
        a = something(segs[k - 1].useable_stop, segs[k - 1].stop_time, last(segs[k - 1].times))
        b = something(segs[k].useable_start, segs[k].start_time, first(segs[k].times))
        _seconds(a, b) < 0 && push!(problems, "segments $(k - 1) and $k overlap; SPICE would take segment $k " *
                                    "wherever both cover a time. Set their usable spans to separate them")
    end
    isempty(problems) || throw(ArgumentError("write_spk: nothing written; $(length(problems)) problem$(length(problems) == 1 ? "" : "s"):\n  " *
                                             join(problems, "\n  ")))
    dir = dirname(abspath(path))
    isdir(dir) || throw(ArgumentError("write_spk: folder $dir does not exist"))
    tmp = joinpath(dir, ".spk_" * string(rand(UInt64); base = 16) * ".bsp")   # SPICE will not open an existing file
    lock(_SPICE_LOCK) do
        handle = nothing
        try
            handle = SPICE.spkopn(tmp, "Epicycle", 0)
            for s in segs
                m = something(s.interpolation, Lagrange(7))
                ets = _et.(s.times)
                first_et = s.useable_start === nothing ? ets[1] : _et(s.useable_start)
                last_et = s.useable_stop === nothing ? ets[end] : _et(s.useable_stop)
                states = [collect(x) for x in s.states]
                centre = _BODIES[uppercase(s.center)].naifid
                frame = _SPK_FRAME_NAMES[s.frame]
                if m isa Hermite
                    SPICE.spkw13(handle, naif_id, centre, frame, first_et, last_et, s.object_name, m.degree, states, ets)
                else
                    SPICE.spkw09(handle, naif_id, centre, frame, first_et, last_et, s.object_name, m.degree, states, ets)
                end
            end
            SPICE.spkcls(handle)
            handle = nothing
            mv(tmp, path; force = true)
        catch
            handle === nothing || try SPICE.spkcls(handle) catch end
            rm(tmp; force = true)
            rethrow()
        end
    end
    return path
end
