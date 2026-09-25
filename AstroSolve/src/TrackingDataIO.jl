# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# Read and write the CCSDS Tracking Data Messages used by estimators.
#
# Supported subset:
#
#   · Format: CCSDS 503.0-B-2 KVN (Keyword=Value Notation)
#   · Observables: RANGE only
#   · Single segment per file
#   · TIME_SYSTEM ∈ {UTC, TAI, TT, TDB, TCG, TCB}
#   · MODE = SEQUENTIAL only
#
# References (CCSDS 503.0-B-2):
#   §3.3.1  Header keywords    : CCSDS_TDM_VERS, CREATION_DATE, ORIGINATOR
#   §3.3.2  Metadata keywords  : TIME_SYSTEM, START_TIME, STOP_TIME,
#                                PARTICIPANT_n, MODE, PATH, RANGE_UNITS
#   §3.3.3  Data keywords      : RANGE
#   Annex B Examples           : KVN structure (META_START/META_STOP,
#                                DATA_START/DATA_STOP)

module TrackingDataIO

using Printf
using AstroEpochs
using AstroEpochs: Time, UTC, TAI, TT, TDB, TCG, TCB, ISOT
using AstroEpochs.Tempo: jd2cal

# =============================================================================
# 1. Format tag
# =============================================================================

abstract type AbstractTrackingDataFormat end

"""
    CCSDS_KVN()

Format tag for a CCSDS 503.0-B-2 Tracking Data Message written in Keyword-Value
Notation (KVN).

# Example
```julia
using AstroSolve

TrackingDataFile("tracking.tdm", CCSDS_KVN())
```
"""
struct CCSDS_KVN <: AbstractTrackingDataFormat end

# =============================================================================
# 2. ObservationRecord
# =============================================================================

"""
    ObservationRecord(measurement_type, t_receive, observed, computed,
                      participant_1)

A tracking observation and the participant that supplied it.

# Fields
- `measurement_type::Symbol`: Observable type, such as `:RANGE` or `:DOPPLER`.
- `t_receive::Time`: Receive epoch in the segment's time scale.
- `observed::Float64`: Measured value.
- `computed::Float64`: Model-predicted value, or `NaN` when it has not been set.
- `participant_1::String`: First participant, typically the receiving ground station.

# Example
```julia
using AstroEpochs: Time, TT, ISOT
using AstroSolve: ObservationRecord

epoch = Time("2020-03-01T00:00:00.000", TT(), ISOT())

ObservationRecord(:RANGE, epoch, 2015.3274, "DSS-14")
```
"""
struct ObservationRecord
    measurement_type::Symbol
    t_receive::Time
    observed::Float64
    computed::Float64
    participant_1::String
end

# 4-arg shortcut: computed = NaN.
ObservationRecord(measurement_type::Symbol, t_receive::Time,
                  observed::Real, participant_1::AbstractString) =
    ObservationRecord(measurement_type, t_receive, Float64(observed),
                      NaN, String(participant_1))

# 5-arg explicit form (positional, allows Real-typed observed/computed).
ObservationRecord(measurement_type::Symbol, t_receive::Time,
                  observed::Real, computed::Real,
                  participant_1::AbstractString) =
    ObservationRecord(measurement_type, t_receive, Float64(observed),
                      Float64(computed), String(participant_1))

# =============================================================================
# 3. TrackingDataFile
# =============================================================================

"""
    TrackingDataFile(path, format = CCSDS_KVN())

The path and format of a tracking data file. It carries no open I/O state; `format` is what
[`read_records`](@ref) and [`write_records`](@ref) use to choose a reader and a writer.

# Fields
- `path::String`: the file's path on disk.
- `format::AbstractTrackingDataFormat`: the file format, which selects the reader and writer.
  [`CCSDS_KVN`](@ref) unless another is given.

# Example
```julia
using AstroSolve

TrackingDataFile("tracking.tdm", CCSDS_KVN())
```
"""
struct TrackingDataFile{F<:AbstractTrackingDataFormat}
    path::String
    format::F
end

TrackingDataFile(path::AbstractString) = TrackingDataFile(String(path), CCSDS_KVN())
TrackingDataFile(path::AbstractString, fmt::AbstractTrackingDataFormat) =
    TrackingDataFile{typeof(fmt)}(String(path), fmt)

# =============================================================================
# 4. Header / segment-metadata structs
# =============================================================================

"""
    TDMHeader(creation_date::String, originator::String)

CCSDS TDM header keywords (`CREATION_DATE`, `ORIGINATOR`).
`creation_date` is an ISO8601 string; the writer copies it verbatim.

# Fields
- `creation_date::String`: the `CREATION_DATE` keyword, written out verbatim.
- `originator::String`: the `ORIGINATOR` keyword, naming who produced the file.

# Example
```julia
TDMHeader("2020-03-01T00:00:00", "GenAstro")
```
"""
struct TDMHeader
    creation_date::String
    originator::String
end

"""
    TDMSegmentMeta(time_system, participant_1, participant_2,
                    mode, path, range_units)

Metadata for one CCSDS TDM segment.

# Fields
- `time_system::String`: Time system used by the segment epochs.
- `participant_1::String`: First participant, typically a ground station.
- `participant_2::String`: Second participant, typically the tracked spacecraft.
- `mode::String`: Tracking mode. The supported value is `"SEQUENTIAL"`.
- `path::String`: Ordered signal path through the participants.
- `range_units::String`: Units used for range observations.

# Notes
The writer derives `START_TIME` and `STOP_TIME` from the records. Defaults of
`mode = "SEQUENTIAL"`, `path = "1,2,1"`, and `range_units = "km"` match the
orbit-determination simulator.

# Example
```julia
using AstroSolve: TDMSegmentMeta

TDMSegmentMeta(time_system = "TT", participant_1 = "DSS-14", participant_2 = "Sat")
```
"""
struct TDMSegmentMeta
    time_system::String
    participant_1::String
    participant_2::String
    mode::String
    path::String
    range_units::String
end

function TDMSegmentMeta(; time_system::AbstractString,
                          participant_1::AbstractString,
                          participant_2::AbstractString,
                          mode::AbstractString = "SEQUENTIAL",
                          path::AbstractString = "1,2,1",
                          range_units::AbstractString = "km")
    return TDMSegmentMeta(String(time_system),
                          String(participant_1),
                          String(participant_2),
                          String(mode),
                          String(path),
                          String(range_units))
end

# =============================================================================
# 5. Constants — the supported subset of CCSDS 503.0-B-2
# =============================================================================

const TDM_VERSION = "2.0"

# CCSDS time systems we understand AND AstroEpochs models.  GPS is
# allowed by the spec but not modeled in AstroEpochs yet, so it's not
# in this set.
const ALLOWED_TIME_SYSTEMS = ("UTC", "TAI", "TT", "TDB", "TCG", "TCB")

const ALLOWED_MODES        = ("SEQUENTIAL",)            # v1 subset
const ALLOWED_RANGE_UNITS  = ("km", "s", "RU")          # spec-allowed; we just store
const SUPPORTED_OBSERVABLES = (:RANGE, :DOPPLER)        # v1 subset

# Required metadata keys for the v1 subset.  Reader checks presence;
# writer always emits these.
const REQUIRED_META_KEYS = ("TIME_SYSTEM", "START_TIME", "STOP_TIME",
                            "PARTICIPANT_1", "PARTICIPANT_2",
                            "MODE", "PATH", "RANGE_UNITS")

# =============================================================================
# 6. Internal helpers
# =============================================================================

@inline _split_kvn(line::AbstractString) = begin
    parts = split(line, "="; limit = 2)
    length(parts) == 2 || throw(ArgumentError(
        "CCSDS_KVN: every metadata and header line must be KEY = VALUE; " *
        "got $(repr(line))"))
    return strip(parts[1]), strip(parts[2])
end

function _scale_tag(time_system::AbstractString)
    s = uppercase(strip(time_system))
    if !(s in ALLOWED_TIME_SYSTEMS)
        allowed = join(ALLOWED_TIME_SYSTEMS, ", ")
        throw(ArgumentError(
            "TIME_SYSTEM must be one of $allowed; got $(repr(time_system))"))
    end
    s == "UTC" ? UTC() :
    s == "TAI" ? TAI() :
    s == "TT"  ? TT()  :
    s == "TDB" ? TDB() :
    s == "TCG" ? TCG() :
    s == "TCB" ? TCB() :
    # Unreachable: the membership test above admits only tags this chain maps. It is here so
    # that adding a tag to ALLOWED_TIME_SYSTEMS without adding it here fails loudly (§9.7)
    # rather than falling off the end of the chain.
    throw(ArgumentError(                                              # COV_EXCL_LINE
        "TIME_SYSTEM $(repr(s)) is listed in ALLOWED_TIME_SYSTEMS " *  # COV_EXCL_LINE
        "but has no AstroEpochs scale tag; add one to _scale_tag"))    # COV_EXCL_LINE
end

# Compare Times by absolute JD (works across scales for sorting; OK
# because all records in one segment share a scale by construction).
_jd_total(t::Time) = t.jd1 + t.jd2

# Format `t` as a millisecond-precision ISOT string with proper carry
# from 60.000 s into the next minute (avoids a bug in AstroEpochs'
# `_to_isot` which can emit `:60.000` for boundary fractions).
function _format_isot(t::Time)
    # Calendar date + civil-day fraction (fd ∈ [0, 1)).
    y, m, d, fd = jd2cal(t.jd1, t.jd2)

    # Snap fd to the millisecond grid.
    ms = round(Int, fd * 86_400_000)
    if ms >= 86_400_000
        # Boundary case: carry into next civil day.
        ms -= 86_400_000
        y, m, d, _ = jd2cal(t.jd1 + 1.0, t.jd2)
    elseif ms < 0
        # Unreachable through jd2cal, which normalizes the day fraction into [0, 1) so `ms` is
        # never negative. Kept as the mirror of the forward carry above: if a future epoch
        # backend returns a signed fraction, this carries rather than emitting a negative field.
        ms += 86_400_000                                  # COV_EXCL_LINE
        y, m, d, _ = jd2cal(t.jd1 - 1.0, t.jd2)           # COV_EXCL_LINE
    end

    h,  rem1 = divrem(ms,   3_600_000)
    mi, rem2 = divrem(rem1, 60_000)
    s,  msr  = divrem(rem2, 1_000)

    return @sprintf("%04d-%02d-%02dT%02d:%02d:%02d.%03d",
                    y, m, d, h, mi, s, msr)
end

# =============================================================================
# 7. Writer — CCSDS_KVN
# =============================================================================

"""
    write_records(file::TrackingDataFile{CCSDS_KVN}, records, header,
                  metas::AbstractVector{TDMSegmentMeta})
    write_records(file, records, header, meta::TDMSegmentMeta)

Write `records` to `file.path` as a CCSDS 503.0-B-2 KVN tracking data message.

# Arguments
- `file`    — the destination, whose `format` selects the KVN writer.
- `records` — the observations to write. Each contributes one data line.
- `header`  — the `CREATION_DATE` and `ORIGINATOR` keywords.
- `metas`   — one segment per entry. A single `meta` writes one segment.

# Notes
One segment is emitted per meta. Records are routed to a segment by matching their
`participant_1` against `meta.participant_1`, and are sorted by `t_receive` within each segment.
A record with an empty `participant_1` belongs to the only meta when there is one, which keeps
callers working that predate the field; with two or more metas it has no segment and is refused.

A meta that no record matches is dropped rather than written, because a TDM segment with no data
lines is malformed. Only `observed` is written: `computed` is carried through the OD pipeline
rather than by the file.

`START_TIME` and `STOP_TIME` are derived from the records in each segment, not taken from the
meta, so they always describe what was actually written.

# Returns
`file`, unchanged, so a write can be piped straight into `read_records`.

Throws `ArgumentError` when `records` or `metas` is empty, when a meta names a `TIME_SYSTEM` or
`MODE` outside the supported subset, when two metas share a `participant_1`, when a record's
`measurement_type` is not in `SUPPORTED_OBSERVABLES`, when a record's `participant_1` matches no
meta, or when no record matches any meta.

# Example
```julia
using AstroEpochs: Time, TT, ISOT
using AstroSolve

t  = Time("2020-03-01T00:00:00.000", TT(), ISOT())
recs = [ObservationRecord(:RANGE, t, 42164.0, "DSS-14")]
meta = TDMSegmentMeta(time_system = "TT", participant_1 = "DSS-14",
                      participant_2 = "SAT-1")
file = TrackingDataFile(tempname())
write_records(file, recs, TDMHeader("2020-03-01T00:00:00", "GEN ASTRO"), meta)
back, header, metas = read_records(file)
```
"""
function write_records(file::TrackingDataFile{CCSDS_KVN},
                       records::AbstractVector{ObservationRecord},
                       header::TDMHeader,
                       metas::AbstractVector{TDMSegmentMeta})
    isempty(records) && throw(ArgumentError(
        "write_records: records must hold at least one observation; got an empty vector"))
    isempty(metas) && throw(ArgumentError(
        "write_records: metas must hold at least one segment; got an empty vector"))

    for meta in metas
        meta.time_system in ALLOWED_TIME_SYSTEMS ||
            throw(ArgumentError(
                "write_records: meta TIME_SYSTEM must be one of " *
                "$(join(ALLOWED_TIME_SYSTEMS, ", ")); got $(repr(meta.time_system))"))
        meta.mode in ALLOWED_MODES ||
            throw(ArgumentError(
                "write_records: meta MODE must be one of " *
                "$(join(ALLOWED_MODES, ", ")); got $(repr(meta.mode))"))
    end

    # No two metas may share participant_1 (segment routing would be ambiguous).
    let seen = Set{String}()
        for meta in metas
            meta.participant_1 in seen &&
                throw(ArgumentError(
                    "write_records: each segment meta must carry a distinct participant_1, " *
                    "since records are routed to segments by that name; " *
                    "got $(repr(meta.participant_1)) twice"))
            push!(seen, meta.participant_1)
        end
    end

    for (i, r) in pairs(records)
        r.measurement_type in SUPPORTED_OBSERVABLES ||
            throw(ArgumentError(
                "write_records: record $i measurement_type must be one of " *
                "$(join(SUPPORTED_OBSERVABLES, ", ")); got $(repr(r.measurement_type))"))
    end

    # Group records by participant_1.  Empty participant_1 is only
    # acceptable when there is exactly one meta (back-compat shim).
    groups = Dict{String, Vector{ObservationRecord}}()
    for meta in metas
        groups[meta.participant_1] = ObservationRecord[]
    end
    for (i, r) in pairs(records)
        key = r.participant_1
        if isempty(key)
            length(metas) == 1 ||
                throw(ArgumentError(
                    "write_records: record $i must name its participant_1 when more than one " *
                    "segment meta is given, since there is no single segment to route it to; " *
                    "got an empty name with $(length(metas)) metas"))
            key = metas[1].participant_1
        end
        haskey(groups, key) ||
            throw(ArgumentError(
                "write_records: record $i participant_1 must match a segment meta; " *
                "got $(repr(key)), and the metas name " *
                "$(join((repr(m.participant_1) for m in metas), ", "))"))
        push!(groups[key], r)
    end

    # Drop metas whose group is empty — a TDM segment with no data
    # lines is malformed.  If *all* groups are empty, that's an error
    # because the user asked us to write nothing.
    nonempty_metas = TDMSegmentMeta[m for m in metas
                                    if !isempty(groups[m.participant_1])]
    isempty(nonempty_metas) &&
        throw(ArgumentError(
            "write_records: at least one record must match a segment meta, or the file would " *
            "have no data lines; none of the $(length(records)) records matched"))

    open(file.path, "w") do io
        # Header
        println(io, "CCSDS_TDM_VERS = ", TDM_VERSION)
        println(io, "CREATION_DATE = ", header.creation_date)
        println(io, "ORIGINATOR    = ", header.originator)
        println(io)

        for (seg_idx, meta) in pairs(nonempty_metas)
            sorted  = sort(groups[meta.participant_1];
                           by = r -> _jd_total(r.t_receive))
            t_start = sorted[1].t_receive
            t_stop  = sorted[end].t_receive

            # Metadata
            println(io, "META_START")
            println(io, "TIME_SYSTEM   = ", meta.time_system)
            println(io, "START_TIME    = ", _format_isot(t_start))
            println(io, "STOP_TIME     = ", _format_isot(t_stop))
            println(io, "PARTICIPANT_1 = ", meta.participant_1)
            println(io, "PARTICIPANT_2 = ", meta.participant_2)
            println(io, "MODE          = ", meta.mode)
            println(io, "PATH          = ", meta.path)
            println(io, "RANGE_UNITS   = ", meta.range_units)
            println(io, "META_STOP")
            println(io)

            # Data
            println(io, "DATA_START")
            for r in sorted
                # Always write `observed`.  `computed` is round-trip-only
                # via OD pipeline metadata, not part of TDM data lines.
                println(io, String(r.measurement_type), " = ",
                        _format_isot(r.t_receive), " ", r.observed)
            end
            println(io, "DATA_STOP")
            seg_idx == length(nonempty_metas) || println(io)
        end
    end

    return file
end

# Single-segment convenience overload.  All records (if tagged) must
# match meta.participant_1.
write_records(file::TrackingDataFile{CCSDS_KVN},
              records::AbstractVector{ObservationRecord},
              header::TDMHeader,
              meta::TDMSegmentMeta) =
    write_records(file, records, header, TDMSegmentMeta[meta])

# =============================================================================
# 8. Reader — CCSDS_KVN
# =============================================================================

"""
    read_records(file::TrackingDataFile{CCSDS_KVN})
        -> (records::Vector{ObservationRecord},
            header::TDMHeader,
            metas::Vector{TDMSegmentMeta})

Read `file.path` as a CCSDS 503.0-B-2 KVN tracking data message.

# Arguments
- `file` — the source, whose `format` selects the KVN reader.

# Notes
A file may carry several META/DATA segments, which is how one file holds passes from more than
one ground station. Each record is tagged with its enclosing segment's `PARTICIPANT_1`, so the
station a measurement came from survives the read and an estimator can pick the matching
measurement specification.

Records come back in file order, segment by segment and line by line within a segment. That is
not necessarily time order across segments, since stations are usually written one after another
rather than interleaved; `build_od_closures` sorts by epoch for this reason.

Each record's epoch is parsed in the time scale its segment declares, so a file mixing scales
across segments reads correctly.

`COMMENT` lines are skipped wherever they appear, which CCSDS allows anywhere in a message. Blank
lines are skipped too.

# Returns
`(records, header, metas)` — the observations, the header keywords, and one `TDMSegmentMeta` per
segment in file order.

Throws `ArgumentError` when the file does not exist, when `CCSDS_TDM_VERS` is missing or is not
the supported version, when a META segment omits a required key, when `TIME_SYSTEM` or `MODE` is
outside the supported subset, when the data section names an unsupported observable, when a
key=value or data line is malformed, when a META or DATA segment is unterminated or nested, or
when the file carries no META segment or no observations.

# Example
<!-- doc-fragment -->
```julia
file = TrackingDataFile("pass.tdm")
records, header, metas = read_records(file)
epochs  = [r.t_receive for r in records]
station = metas[1].participant_1
```
"""
function read_records(file::TrackingDataFile{CCSDS_KVN})
    isfile(file.path) || throw(ArgumentError(
        "read_records: file.path must name a readable file; got $(repr(file.path))"))

    creation_date = ""
    originator    = ""
    saw_version   = false

    metas    = TDMSegmentMeta[]
    records  = ObservationRecord[]

    # Per-segment scratch state
    in_meta        = false
    in_data        = false
    cur_meta_dict  = Dict{String,String}()
    cur_meta       = nothing            # ::Union{Nothing, TDMSegmentMeta}
    cur_scale_tag  = UTC()
    cur_part1      = ""

    open(file.path, "r") do io
    for (lineno, raw) in enumerate(eachline(io))
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "COMMENT") && continue

        if line == "META_START"
            in_meta && throw(ArgumentError(
                "read_records: $(file.path) line $lineno: META_START must follow a META_STOP; " *
                "the previous META segment is still open"))
            in_data && throw(ArgumentError(
                "read_records: $(file.path) line $lineno: META_START must not appear inside a " *
                "DATA segment; close it with DATA_STOP first"))
            in_meta = true
            cur_meta_dict = Dict{String,String}()
            continue
        elseif line == "META_STOP"
            in_meta || throw(ArgumentError(
                "read_records: $(file.path) line $lineno: META_STOP must close a segment " *
                "opened by META_START; no META segment is open"))
            in_meta = false

            for k in REQUIRED_META_KEYS
                haskey(cur_meta_dict, k) ||
                    throw(ArgumentError(
                        "read_records: $(file.path): the META segment must define $k; " *
                        "required keys are $(join(REQUIRED_META_KEYS, ", "))"))
            end
            cur_meta_dict["MODE"] in ALLOWED_MODES ||
                throw(ArgumentError(
                    "read_records: $(file.path): MODE must be one of " *
                    "$(join(ALLOWED_MODES, ", ")); got $(repr(cur_meta_dict["MODE"]))"))
            cur_scale_tag = _scale_tag(cur_meta_dict["TIME_SYSTEM"])
            cur_part1     = cur_meta_dict["PARTICIPANT_1"]
            cur_meta = TDMSegmentMeta(
                time_system   = cur_meta_dict["TIME_SYSTEM"],
                participant_1 = cur_meta_dict["PARTICIPANT_1"],
                participant_2 = cur_meta_dict["PARTICIPANT_2"],
                mode          = cur_meta_dict["MODE"],
                path          = cur_meta_dict["PATH"],
                range_units   = cur_meta_dict["RANGE_UNITS"],
            )
            push!(metas, cur_meta)
            continue
        elseif line == "DATA_START"
            in_meta && throw(ArgumentError(
                "read_records: $(file.path) line $lineno: DATA_START must follow META_STOP; " *
                "the META segment is still open"))
            in_data && throw(ArgumentError(
                "read_records: $(file.path) line $lineno: DATA_START must follow a DATA_STOP; " *
                "the previous DATA segment is still open"))
            cur_meta === nothing && throw(ArgumentError(
                "read_records: $(file.path) line $lineno: DATA_START must be preceded by a " *
                "META segment, which says what time system and participants the data uses"))
            in_data = true
            continue
        elseif line == "DATA_STOP"
            in_data || throw(ArgumentError(
                "read_records: $(file.path) line $lineno: DATA_STOP must close a segment " *
                "opened by DATA_START; no DATA segment is open"))
            in_data = false
            continue
        end

        if in_meta
            k, v = _split_kvn(line)
            cur_meta_dict[k] = v
        elseif in_data
            k, v = _split_kvn(line)
            sym = Symbol(k)
            sym in SUPPORTED_OBSERVABLES ||
                throw(ArgumentError(
                    "read_records: $(file.path) line $lineno: the observable must be one of " *
                    "$(join(SUPPORTED_OBSERVABLES, ", ")); got $k"))

            parts = split(v)
            length(parts) == 2 ||
                throw(ArgumentError(
                    "read_records: $(file.path) line $lineno: a $k line must be " *
                    "$k = <epoch> <value>; got $(repr(line))"))
            t_iso, val_str = parts
            t   = Time(String(t_iso), cur_scale_tag, ISOT())
            val = parse(Float64, val_str)
            push!(records,
                  ObservationRecord(sym, t, val, NaN, cur_part1))
        else
            # Header lines
            k, v = _split_kvn(line)
            if k == "CCSDS_TDM_VERS"
                saw_version = true
                v == TDM_VERSION ||
                    throw(ArgumentError(
                        "read_records: $(file.path): CCSDS_TDM_VERS must be $TDM_VERSION; " *
                        "got $(repr(v))"))
            elseif k == "CREATION_DATE"
                creation_date = v
            elseif k == "ORIGINATOR"
                originator = v
            else
                throw(ArgumentError(
                    "read_records: $(file.path) line $lineno: a line outside META and DATA " *
                    "segments must be one of CCSDS_TDM_VERS, CREATION_DATE, ORIGINATOR; " *
                    "got $k"))
            end
        end
    end
    end  # open

    in_meta && throw(ArgumentError(
        "read_records: $(file.path): every META_START must be closed by META_STOP; " *
        "the file ends inside a META segment"))
    in_data && throw(ArgumentError(
        "read_records: $(file.path): every DATA_START must be closed by DATA_STOP; " *
        "the file ends inside a DATA segment"))
    saw_version || throw(ArgumentError(
        "read_records: $(file.path): the file must open with CCSDS_TDM_VERS = $TDM_VERSION; " *
        "no version line was found"))
    isempty(metas) && throw(ArgumentError(
        "read_records: $(file.path): the file must carry at least one META segment, " *
        "which says what time system and participants the data uses; none was found"))
    isempty(records) && throw(ArgumentError(
        "read_records: $(file.path): the file must carry at least one observation; " *
        "supported observables are $(join(SUPPORTED_OBSERVABLES, ", "))"))

    return records, TDMHeader(creation_date, originator), metas
end

# =============================================================================
# 9. Exports
# =============================================================================

export ObservationRecord,
       TrackingDataFile,
       AbstractTrackingDataFormat, CCSDS_KVN,
       TDMHeader, TDMSegmentMeta,
       write_records, read_records,
       TDM_VERSION,
       ALLOWED_TIME_SYSTEMS, ALLOWED_MODES, ALLOWED_RANGE_UNITS,
       SUPPORTED_OBSERVABLES, REQUIRED_META_KEYS

end # module TrackingDataIO
