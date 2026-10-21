# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# CCSDS Orbit Ephemeris Message, keyword = value form (CCSDS 502.0-B-3, §5 and §7). Section
# numbers in comments refer to that document.

const _OEM_VERSIONS = ("1.0", "2.0", "3.0")             # §7.9.1
const _OEM_WRITE_VERSIONS = ("2.0", "3.0")
const _OEM_HEADER_KEYS = ("CCSDS_OEM_VERS", "CLASSIFICATION", "CREATION_DATE", "ORIGINATOR", "MESSAGE_ID")
const _OEM_META_KEYS = ("OBJECT_NAME", "OBJECT_ID", "CENTER_NAME", "REF_FRAME", "REF_FRAME_EPOCH",
                        "TIME_SYSTEM", "START_TIME", "USEABLE_START_TIME", "USEABLE_STOP_TIME",
                        "STOP_TIME", "INTERPOLATION", "INTERPOLATION_DEGREE")
const _OEM_COV_KEYS = ("EPOCH", "COV_REF_FRAME")
# Orbit-relative frames a covariance may be given in, besides the ephemeris frames.
const _COV_RELATIVE_FRAMES = ("RTN", "RSW", "RIC", "QSW", "TNW", "VNB", "LVLH")
const _MAX_LINE = 254                                         # §7.3.2
const _DEFAULT_INTERPOLATION = Lagrange(7)
const _KEY_WIDTH = 20

# ============================================================================================
# Reading
# ============================================================================================

"""
    read_oem(path; strict = false, skip_bad_lines = false,
             ref_frame = nothing, center = nothing, time_system = nothing) -> EphemerisFileData

Read a CCSDS Orbit Ephemeris Message in keyword = value form, versions 1.0, 2.0 and 3.0.

# Arguments
- `path`: The file to read.
- `strict`: Treat every correction as an error. Default `false`.
- `skip_bad_lines`: Drop a data line that cannot be read, with a warning, and split its segment
  there so interpolation never spans the gap. Default `false`: such a line is an error.
- `ref_frame`, `center`, `time_system`: Names to use in place of the file's REF_FRAME,
  CENTER_NAME and TIME_SYSTEM, in every segment. Each override is recorded in the diagnostics
  with the file's original value.

# Notes
- Problems the reader can correct without changing what the data mean are corrected and
  recorded: keywords out of order or in lower case, tabs, a byte-order mark, misplaced comments,
  an exact duplicate record, Fortran `D` exponents, day-of-year epochs, or START_TIME and
  STOP_TIME that disagree with the data. They are returned in `diagnostics` and summarized in
  one warning.
- Problems that would change or guess at the data stop the read with an
  [`EphemerisReadError`](@ref) listing every one, with line numbers: a missing or unknown frame,
  centre or time system, a data line that cannot be read, epochs that go backwards, different
  states at one epoch, or a segment too short for its declared interpolation.

# Returns
An [`EphemerisFileData`](@ref).

# Examples
```jldoctest
path = tempname() * ".oem"
write(path, \"\"\"
CCSDS_OEM_VERS = 3.0
CREATION_DATE = 2026-10-21T00:00:00
ORIGINATOR = EXAMPLE
META_START
OBJECT_NAME = SAT
OBJECT_ID = 2026-001A
CENTER_NAME = EARTH
REF_FRAME = EME2000
TIME_SYSTEM = UTC
START_TIME = 2026-10-21T00:00:00
STOP_TIME = 2026-10-21T00:02:00
META_STOP
2026-10-21T00:00:00 7000.0 0.0 0.0 0.0 7.5 0.0
2026-10-21T00:01:00 6998.4 449.9 0.0 -0.05 7.49 0.0
2026-10-21T00:02:00 6993.6 899.3 0.0 -0.1 7.48 0.0
\"\"\")
f = EpicycleIO.read_oem(path)
length(f.segments[1].times)

# output

3
```
"""
function read_oem(path::AbstractString; strict::Bool = false, skip_bad_lines::Bool = false,
                  ref_frame = nothing, center = nothing, time_system = nothing)
    isfile(path) || throw(ArgumentError("read_oem: no file at $path"))
    text = read(path, String)
    d = _Diagnostics(path; strict)
    if startswith(text, "﻿")
        text = text[nextind(text, 1):end]
        _fixed!(d, 1, "file", "byte-order mark at the start ignored")
    end
    p = _OEMParser(d, EphemerisFileData(format = :oem), skip_bad_lines,
                   ref_frame, center, time_system)
    for (i, raw) in enumerate(split(text, r"\r\n|\n\r|\r|\n"))
        _oem_line!(p, i, raw)
    end
    _oem_end!(p, length(split(text, r"\r\n|\n\r|\r|\n")))
    hint = any(x -> x.severity == :error && x.field in ("REF_FRAME", "CENTER_NAME", "TIME_SYSTEM"), d.items) ?
        "Fix the file, or override a field: read_ephemeris(path; ref_frame = \"…\", center = \"…\", time_system = \"…\")" :
        (any(x -> x.severity == :error && x.field == "data", d.items) && !skip_bad_lines ?
         "To drop unreadable data lines and split the segment there, read with skip_bad_lines = true" : "")
    p.file.diagnostics = _finish(d; hint)
    return p.file
end

# One segment under construction: metadata as raw strings with their lines, then data.
mutable struct _OEMSegment
    start_line::Int
    meta::Dict{String, Tuple{String, Int}}
    meta_order::Vector{String}
    seg::EphemerisSegmentData
    ready::Bool                       # metadata finished and valid enough to read data
    data_lines::Vector{Int}
    breaks::Vector{Int}               # indices where a skipped line splits the data
    skipped::Vector{Int}              # lines dropped at those breaks
    usable::Bool                      # metadata complete: the segment can be returned
    has_acc::Union{Nothing, Bool}
    cov_epoch::Union{Nothing, Tuple{Time, Int}}
    cov_frame::String
    cov_values::Vector{Float64}
    cov_line::Int
    cov_shape_noted::Bool
end
_OEMSegment(line) = _OEMSegment(line, Dict{String, Tuple{String, Int}}(), String[],
                                EphemerisSegmentData(), false, Int[], Int[], Int[], false, nothing,
                                nothing, "", Float64[], 0, false)

mutable struct _OEMParser
    d::_Diagnostics
    file::EphemerisFileData
    skip_bad_lines::Bool
    ref_frame::Any
    center::Any
    time_system::Any
    state::Symbol                     # :start, :header, :meta, :data, :cov, :after_cov
    header_order::Vector{String}
    current::Union{Nothing, _OEMSegment}
    noted::Set{Symbol}                # file-wide corrections reported once
    first_system::String
end
_OEMParser(d, f, s, rf, c, ts) = _OEMParser(d, f, s, rf, c, ts, :start, String[], nothing, Set{Symbol}(), "")

function _once!(p::_OEMParser, kind::Symbol, line, field, msg)
    kind in p.noted && return
    push!(p.noted, kind)
    _fixed!(p.d, line, field, msg)
end

function _oem_line!(p::_OEMParser, i::Int, raw::AbstractString)
    line = raw
    if occursin('\t', line)
        _once!(p, :tab, i, "file", "tab characters read as spaces (first on this line); §7.3.4 allows none")
        line = replace(line, '\t' => ' ')
    end
    if any(c -> !(isascii(c) && (isprint(c) || c == ' ')), line)
        _once!(p, :ascii, i, "file", "characters other than printable ASCII (first on this line); §7.3.4")
    end
    length(line) > _MAX_LINE &&
        _once!(p, :length, i, "file", "lines longer than $_MAX_LINE characters (first on this line); §7.3.2")
    s = strip(line)
    isempty(s) && return                                                  # §7.3.5

    if startswith(uppercase(s), "COMMENT") && (length(s) == 7 || s[8] == ' ')
        return _oem_comment!(p, i, length(s) > 8 ? String(strip(s[9:end])) : "")
    end
    u = uppercase(s)
    if u in ("META_START", "META_STOP", "COVARIANCE_START", "COVARIANCE_STOP")
        u == s || _once!(p, Symbol("case_", u), i, u, "$s in lower case; read as $u")
        return _oem_marker!(p, i, u)
    end
    if p.state in (:start, :header, :meta) || (p.state == :cov && occursin('=', s))
        return _oem_kvn!(p, i, s)
    end
    p.state == :cov && return _oem_cov_values!(p, i, s)
    p.state == :data && return _oem_data!(p, i, s)
    p.state == :after_cov && return _error!(p.d, i, "data",
        "ephemeris line after COVARIANCE_STOP; data must come before the segment's covariance (§5.2.1)")
end

function _oem_comment!(p::_OEMParser, i, text)
    if p.state in (:start, :header)
        p.state == :header && !isempty(p.header_order) && p.header_order[end] != "CCSDS_OEM_VERS" &&
            _once!(p, :hcomment, i, "COMMENT", "header comment after other header keywords; " *
                   "§7.8.9 places them right after CCSDS_OEM_VERS")
        push!(p.file.comments, text)
    elseif p.state == :meta
        isempty(p.current.meta_order) ||
            _once!(p, :mcomment, i, "COMMENT", "metadata comment after other keywords; §7.8.9 places them right after META_START")
        push!(p.current.seg.comments, text)
    elseif p.state == :data
        isempty(p.current.seg.times) ||
            _once!(p, :dcomment, i, "COMMENT", "comment between ephemeris lines; kept, but §7.8.9 allows them only before the first")
        push!(p.current.seg.comments, text)
    else
        push!(p.current.seg.comments, text)
    end
end

function _oem_marker!(p::_OEMParser, i, u)
    if u == "META_START"
        p.state == :start && _error!(p.d, i, "CCSDS_OEM_VERS",
            "the file has no header; the first line must be CCSDS_OEM_VERS = x.y (§7.3.6)")
        p.state == :meta && _error!(p.d, i, "META_START", "META_START inside a metadata block; META_STOP missing")
        p.state == :cov && _error!(p.d, i, "COVARIANCE_STOP", "covariance block not closed before META_START")
        p.state in (:header, :start) && _oem_header_done!(p, i)
        p.current === nothing || _oem_close_segment!(p)
        p.current = _OEMSegment(i)
        p.state = :meta
    elseif u == "META_STOP"
        p.state == :meta || return _error!(p.d, i, "META_STOP", "META_STOP without META_START")
        _oem_meta_done!(p, i)
        p.state = :data
    elseif u == "COVARIANCE_START"
        p.state == :data || return _error!(p.d, i, "COVARIANCE_START",
            "COVARIANCE_START outside a segment's data (§5.2.5.1)")
        p.state = :cov
    else
        p.state == :cov || return _error!(p.d, i, "COVARIANCE_STOP", "COVARIANCE_STOP without COVARIANCE_START")
        _oem_flush_cov!(p, i)
        p.state = :after_cov
    end
end

function _oem_kvn!(p::_OEMParser, i, s)
    parts = split(s, "="; limit = 2)
    if length(parts) != 2
        p.state == :start && return _error!(p.d, i, "CCSDS_OEM_VERS",
            "the first line must be CCSDS_OEM_VERS = x.y (§7.3.6); this is not an OEM header")
        return _error!(p.d, i, "file", "expected KEYWORD = value; got $(repr(String(s)))")
    end
    key, value = String(strip(parts[1])), String(strip(parts[2]))
    ukey = uppercase(key)
    ukey == key || _fixed!(p.d, i, ukey, "keyword in lower case; read as $ukey (§7.4.4)")
    if p.state == :start
        ukey == "CCSDS_OEM_VERS" || _error!(p.d, i, "CCSDS_OEM_VERS",
            "the first line must be CCSDS_OEM_VERS = x.y (§7.3.6); got $ukey")
        p.state = :header
    end
    p.state == :header && return _oem_header_kvn!(p, i, ukey, value)
    p.state == :meta && return _oem_meta_kvn!(p, i, ukey, value)
    return _oem_cov_kvn!(p, i, ukey, value)
end

function _oem_header_kvn!(p::_OEMParser, i, key, value)
    if !(key in _OEM_HEADER_KEYS)
        return _fixed!(p.d, i, key, "not an OEM header keyword; ignored.$(_did_you_mean(key, _OEM_HEADER_KEYS))")
    end
    key in p.header_order && return _error!(p.d, i, key, "$key given twice in the header")
    _check_order!(p, i, key, p.header_order, _OEM_HEADER_KEYS, "header")
    push!(p.header_order, key)
    f = p.file
    if key == "CCSDS_OEM_VERS"
        f.version = value
        value in _OEM_VERSIONS || _fixed!(p.d, i, key,
            "version $value is not 1.0, 2.0 or 3.0 (§7.9.1); read as the current format")
    elseif key == "CREATION_DATE"
        t, notes, problem = _parse_epoch(value, "UTC")
        t === nothing && _fixed!(p.d, i, key, "$problem; kept as text")
        f.creation_date = value
    elseif key == "ORIGINATOR"
        f.originator = value
    elseif key == "MESSAGE_ID"
        f.message_id = value
        f.version in ("1.0", "2.0") && _fixed!(p.d, i, key, "MESSAGE_ID is an OEM 3.0 keyword; read anyway")
    else
        f.classification = value
        f.version in ("1.0", "2.0") && _fixed!(p.d, i, key, "CLASSIFICATION is an OEM 3.0 keyword; read anyway")
    end
end

function _check_order!(p::_OEMParser, i, key, seen, order, what)
    isempty(seen) && return
    last_pos = findfirst(==(seen[end]), order)
    pos = findfirst(==(key), order)
    (last_pos === nothing || pos === nothing || pos > last_pos) && return
    _once!(p, Symbol("order_", what), i, key,
           "$what keywords out of the order §7.4.8 fixes ($key after $(seen[end])); read anyway")
end

function _oem_header_done!(p::_OEMParser, i)
    for key in ("CREATION_DATE", "ORIGINATOR")
        key in p.header_order || _fixed!(p.d, i, key, "mandatory header keyword $key missing (§5.2.2); not needed to read the data")
    end
end

function _oem_meta_kvn!(p::_OEMParser, i, key, value)
    c = p.current
    if !(key in _OEM_META_KEYS)
        return _fixed!(p.d, i, key, "not an OEM metadata keyword; ignored.$(_did_you_mean(key, _OEM_META_KEYS))")
    end
    haskey(c.meta, key) && return _error!(p.d, i, key, "$key given twice in one metadata block")
    _check_order!(p, i, key, c.meta_order, _OEM_META_KEYS, "metadata")
    push!(c.meta_order, key)
    c.meta[key] = (value, i)
end

# Resolve one name: the override if given, else the file's value; record overrides and notes.
function _resolved_name!(p::_OEMParser, c::_OEMSegment, key::String, override, mline)
    file_value, line = get(c.meta, key, ("", mline))
    override === nothing && return file_value, line
    ov = override isa AbstractString ? String(override) : String(nameof(typeof(override)))
    _note!(p.d, line, key, "file $(isempty(file_value) ? "has none" : "says $(repr(file_value))"), read as $ov by override")
    return ov, line
end

function _oem_meta_done!(p::_OEMParser, mline)
    c, seg, d = p.current, p.current.seg, p.d
    for key in ("OBJECT_NAME", "OBJECT_ID")
        if haskey(c.meta, key) && !isempty(c.meta[key][1])
            setfield!(seg, key == "OBJECT_NAME" ? :object_name : :object_id, c.meta[key][1])
        else
            _fixed!(d, mline, key, "mandatory $key missing (§5.2.3); set to UNKNOWN")
        end
    end
    ok = true
    ts, tline = _resolved_name!(p, c, "TIME_SYSTEM", p.time_system, mline)
    if isempty(ts)
        _error!(d, mline, "TIME_SYSTEM", "mandatory TIME_SYSTEM missing (§5.2.3)"); ok = false
    else
        sys, note, problem = _time_system(ts)
        if isempty(problem)
            isempty(note) || _note!(d, tline, "TIME_SYSTEM", note)
            seg.time_system = sys
            if isempty(p.first_system)
                p.first_system = sys
            elseif sys != p.first_system
                _fixed!(d, tline, "TIME_SYSTEM", "$sys differs from the first segment's $(p.first_system); " *
                        "§5.2.4.5 requires one time system per file. Each segment is read in its own.")
            end
        else
            _error!(d, tline, "TIME_SYSTEM", problem); ok = false
        end
    end
    cn, cline = _resolved_name!(p, c, "CENTER_NAME", p.center, mline)
    if isempty(cn)
        _error!(d, mline, "CENTER_NAME", "mandatory CENTER_NAME missing (§5.2.3)"); ok = false
    else
        body, problem = _center_body(cn)
        body === nothing ? (_error!(d, cline, "CENTER_NAME", problem); ok = false) : (seg.center = uppercase(cn))
    end
    fr, fline = _resolved_name!(p, c, "REF_FRAME", p.ref_frame, mline)
    if isempty(fr)
        _error!(d, mline, "REF_FRAME", "mandatory REF_FRAME missing (§5.2.3)"); ok = false
    else
        axes, note, problem = _frame_axes(fr, seg.center)
        if axes === nothing
            _error!(d, fline, "REF_FRAME", problem); ok = false
        else
            isempty(note) || _note!(d, fline, "REF_FRAME", note)
            seg.frame = uppercase(fr)
        end
    end
    if !isempty(seg.time_system)
        for (key, field) in (("START_TIME", :start_time), ("STOP_TIME", :stop_time),
                             ("USEABLE_START_TIME", :useable_start), ("USEABLE_STOP_TIME", :useable_stop),
                             ("REF_FRAME_EPOCH", :frame_epoch))
            haskey(c.meta, key) || continue
            value, line = c.meta[key]
            t, notes, problem = _parse_epoch(value, seg.time_system)
            if t === nothing
                _error!(d, line, key, "$(repr(value)): $problem")
            else
                foreach(n -> _once!(p, Symbol("epoch_", n), line, key, n), notes)
                setfield!(seg, field, t)
            end
        end
        for key in ("START_TIME", "STOP_TIME")
            haskey(c.meta, key) || _fixed!(d, mline, key, "mandatory $key missing (§5.2.3); taken from the data")
        end
    end
    seg.interpolation = _oem_interpolation!(p, c, mline)
    # Data lines are read whenever their epochs can be, so one read reports every problem; a
    # segment with a bad frame or centre is still not returned.
    c.ready = !isempty(seg.time_system)
    c.usable = ok
end

function _oem_interpolation!(p::_OEMParser, c::_OEMSegment, mline)
    d = p.d
    if !haskey(c.meta, "INTERPOLATION")
        haskey(c.meta, "INTERPOLATION_DEGREE") && _fixed!(d, c.meta["INTERPOLATION_DEGREE"][2],
            "INTERPOLATION_DEGREE", "given without INTERPOLATION; ignored")
        return nothing
    end
    name, line = c.meta["INTERPOLATION"]
    u = uppercase(name)
    (u == name || lowercase(name) == name) ||
        _once!(p, :interp_case, line, "INTERPOLATION", "$(repr(name)) in mixed case; §7.5.3 asks for all upper or all lower")
    if !(u in ("LAGRANGE", "HERMITE", "LINEAR"))
        _fixed!(d, line, "INTERPOLATION", "$(repr(name)) is not LAGRANGE, HERMITE or LINEAR; " *
                "$(_DEFAULT_INTERPOLATION) used.$(_did_you_mean(u, ("LAGRANGE", "HERMITE", "LINEAR")))")
        return _DEFAULT_INTERPOLATION
    end
    u == "LINEAR" && return Lagrange(1)
    if !haskey(c.meta, "INTERPOLATION_DEGREE")
        m = u == "HERMITE" ? Hermite(7) : Lagrange(7)
        _fixed!(d, line, "INTERPOLATION_DEGREE", "missing, though INTERPOLATION is given (§5.2.3); degree 7 used")
        return m
    end
    dtext, dline = c.meta["INTERPOLATION_DEGREE"]
    deg = tryparse(Int, dtext)
    if deg === nothing || deg < 1
        _fixed!(d, dline, "INTERPOLATION_DEGREE", "$(repr(dtext)) is not a positive integer; degree 7 used")
        deg = 7
    end
    if u == "HERMITE" && iseven(deg)
        _fixed!(d, dline, "INTERPOLATION_DEGREE", "Hermite degree $deg is even; values and " *
                "derivatives at each point give odd degrees, so degree $(deg + 1) used")
        deg += 1
    end
    return u == "HERMITE" ? Hermite(deg) : Lagrange(deg)
end

# A data line that cannot be used: an error, or with skip_bad_lines a warning and a split.
function _bad_line!(p::_OEMParser, i, msg)
    c = p.current
    if p.skip_bad_lines
        _fixed!(p.d, i, "data", "$msg; line skipped and the segment split here")
        isempty(c.seg.times) || push!(c.breaks, length(c.seg.times) + 1)
        push!(c.skipped, i)
    else
        _error!(p.d, i, "data", msg)
    end
end

function _oem_number(p::_OEMParser, i, tok)
    x = tryparse(Float64, tok)
    x === nothing || return x
    if occursin(r"^[+-]?\d*\.?\d*[dD][+-]?\d+$", tok)
        _once!(p, :fortran, i, "data", "Fortran D exponents (as in $tok) read as E; §7.5.7 allows E or e only")
        return tryparse(Float64, replace(tok, r"[dD]" => "e"))
    end
    return nothing
end

function _oem_data!(p::_OEMParser, i, s)
    c = p.current
    c.ready || return                                       # metadata already reported as unusable
    seg = c.seg
    toks = split(s)
    n = length(toks)
    n in (7, 10) || return _bad_line!(p, i, "data line has $n values; expected 7 (epoch and state) or 10 (with accelerations)")
    t, notes, problem = _parse_epoch(toks[1], seg.time_system)
    t === nothing && return _bad_line!(p, i, "epoch $(repr(String(toks[1]))): $problem")
    foreach(note -> _once!(p, Symbol("epoch_", note), i, "data", note), notes)
    vals = Vector{Float64}(undef, n - 1)
    for k in 2:n
        x = _oem_number(p, i, toks[k])
        x === nothing && return _bad_line!(p, i, "value $(k - 1), $(repr(String(toks[k]))), is not a number")
        isfinite(x) || return _bad_line!(p, i, "value $(k - 1), $(toks[k]), is not finite")
        vals[k - 1] = x
    end
    acc = n == 10
    if c.has_acc === nothing
        c.has_acc = acc
        acc && (seg.accelerations = SVector{3, Float64}[])
    elseif c.has_acc != acc
        return _bad_line!(p, i, "line has $(n) values where this segment's earlier lines have $(c.has_acc ? 10 : 7); " *
                                "accelerations are all or none within a segment")
    end
    state = SVector{6, Float64}(vals[1:6])
    if !isempty(seg.times)
        dt = _seconds(seg.times[end], t)
        if dt == 0
            same = state == seg.states[end] && (!acc || SVector{3}(vals[7:9]) == seg.accelerations[end])
            same && return _fixed!(p.d, i, "data", "duplicate of line $(c.data_lines[end]) at $(toks[1]); one copy kept")
            return _error!(p.d, i, "data", "a different state at the same epoch as line $(c.data_lines[end]) ($(toks[1]))")
        elseif dt < 0
            return _bad_line!(p, i, "epoch $(toks[1]) is before the previous line's ($(c.data_lines[end])); epochs must increase")
        end
    end
    push!(seg.times, t)
    push!(seg.states, state)
    acc && push!(seg.accelerations, SVector{3, Float64}(vals[7:9]))
    push!(c.data_lines, i)
end

function _oem_cov_kvn!(p::_OEMParser, i, key, value)
    c = p.current
    if key == "EPOCH"
        _oem_flush_cov!(p, i)
        t, notes, problem = _parse_epoch(value, c.seg.time_system)
        t === nothing && return _error!(p.d, i, "EPOCH", "$(repr(value)): $problem")
        c.cov_epoch, c.cov_frame, c.cov_line = (t, i), "", i
    elseif key == "COV_REF_FRAME"
        c.cov_epoch === nothing && return _error!(p.d, i, "COV_REF_FRAME", "COV_REF_FRAME before EPOCH")
        u = uppercase(value)
        if !(u in _COV_RELATIVE_FRAMES)
            axes, note, problem = _frame_axes(u, c.seg.center)
            axes === nothing && return _error!(p.d, i, "COV_REF_FRAME", problem)
        end
        c.cov_frame = u
    else
        _error!(p.d, i, key, "not a covariance keyword; expected EPOCH or COV_REF_FRAME (§5.2.5)")
    end
end

function _oem_cov_values!(p::_OEMParser, i, s)
    c = p.current
    c.cov_epoch === nothing && return _error!(p.d, i, "EPOCH", "covariance values before an EPOCH line (§5.2.5.3)")
    toks = split(s)
    row = length(c.cov_values)
    expected = findfirst(k -> k * (k - 1) ÷ 2 == row, 1:6)            # the row this line should be
    if !c.cov_shape_noted && (expected === nothing || length(toks) != expected)
        c.cov_shape_noted = true
        _fixed!(p.d, i, "covariance", "lines not in the lower-triangular shape of §5.2.5.4 (one value, then two, …); read in order")
    end
    for tok in toks
        x = _oem_number(p, i, tok)
        x === nothing && return _error!(p.d, i, "covariance", "$(repr(String(tok))) is not a number")
        push!(c.cov_values, x)
    end
    length(c.cov_values) > 21 && _error!(p.d, i, "covariance",
        "more than 21 values after EPOCH at line $(c.cov_line); a 6×6 lower triangle has 21")
end

function _oem_flush_cov!(p::_OEMParser, i)
    c = p.current
    c.cov_epoch === nothing && return
    vals = c.cov_values
    if length(vals) != 21
        length(vals) < 21 && _error!(p.d, c.cov_line, "covariance",
            "$(length(vals)) values after EPOCH; a 6×6 lower triangle has 21")
    else
        m = zeros(6, 6)
        k = 0
        for r in 1:6, col in 1:r
            k += 1
            m[r, col] = m[col, r] = vals[k]
        end
        bad = findfirst(j -> m[j, j] < 0, 1:6)
        if bad !== nothing
            _error!(p.d, c.cov_line, "covariance", "negative variance in element [$bad,$bad]")
        else
            cov = c.seg.covariances
            t = c.cov_epoch[1]
            if !isempty(cov) && _seconds(cov[end].epoch, t) <= 0
                _fixed!(p.d, c.cov_line, "covariance", "covariance epochs not increasing (§5.2.5.7); sorted")
            end
            push!(cov, CovarianceRecord(t, isempty(c.cov_frame) ? c.seg.frame : c.cov_frame, SMatrix{6, 6}(m)))
            sort!(cov; by = r -> _seconds(cov[1].epoch, r.epoch))
        end
    end
    c.cov_epoch, c.cov_values, c.cov_frame = nothing, Float64[], ""
end

function _oem_end!(p::_OEMParser, nlines)
    p.state == :start && return _error!(p.d, 1, "file", "empty file, or no CCSDS_OEM_VERS line")
    p.state == :meta && _error!(p.d, nlines, "META_STOP", "file ends inside a metadata block")
    p.state == :cov && _error!(p.d, nlines, "COVARIANCE_STOP", "file ends inside a covariance block")
    p.state == :header && _error!(p.d, nlines, "META_START", "the file has a header but no segments")
    p.current === nothing || _oem_close_segment!(p)
    _oem_check_segments!(p)
end

function _oem_close_segment!(p::_OEMParser)
    c, d = p.current, p.d
    p.current = nothing
    c.ready || return
    seg = c.seg
    if isempty(seg.times)
        return _error!(d, c.start_line, "data", "the segment starting at line $(c.start_line) has no ephemeris lines")
    end
    pieces = _split_pieces(seg, c.breaks)
    for (k, piece) in enumerate(pieces)
        _check_piece!(p, c, piece, k == 1, k == length(pieces))
        c.usable && push!(p.file.segments, piece)
    end
end

function _split_pieces(seg::EphemerisSegmentData, breaks)
    isempty(breaks) && return [seg]
    bounds = vcat(1, sort(unique(breaks)), length(seg.times) + 1)
    pieces = EphemerisSegmentData[]
    for k in 1:length(bounds) - 1
        r = bounds[k]:bounds[k + 1] - 1
        isempty(r) && continue
        piece = deepcopy(seg)
        piece.times, piece.states = seg.times[r], seg.states[r]
        seg.accelerations === nothing || (piece.accelerations = seg.accelerations[r])
        a, b = first(piece.times), last(piece.times)
        piece.covariances = filter(cv -> _seconds(a, cv.epoch) >= 0 && _seconds(cv.epoch, b) >= 0, seg.covariances)
        k > 1 && (piece.start_time = a; piece.useable_start = nothing)
        k < length(bounds) - 1 && (piece.stop_time = b; piece.useable_stop = nothing)
        push!(pieces, piece)
    end
    return pieces
end

function _check_piece!(p::_OEMParser, c::_OEMSegment, seg::EphemerisSegmentData, first_piece, last_piece)
    d = p.d
    a, b = first(seg.times), last(seg.times)
    line = c.start_line
    if seg.start_time === nothing
        seg.start_time = a
    elseif _seconds(seg.start_time, a) != 0
        _fixed!(d, get(c.meta, "START_TIME", ("", line))[2], "START_TIME",
                "$(_format_epoch(seg.start_time, seg.time_system, 3)) differs from the first epoch " *
                "$(_format_epoch(a, seg.time_system, 3)); coverage starts at the first epoch")
        seg.start_time = a
    end
    if seg.stop_time === nothing
        seg.stop_time = b
    elseif _seconds(seg.stop_time, b) != 0
        _fixed!(d, get(c.meta, "STOP_TIME", ("", line))[2], "STOP_TIME",
                "$(_format_epoch(seg.stop_time, seg.time_system, 3)) differs from the last epoch " *
                "$(_format_epoch(b, seg.time_system, 3)); coverage ends at the last epoch")
        seg.stop_time = b
    end
    for (field, key, bad) in ((:useable_start, "USEABLE_START_TIME", t -> _seconds(a, t) < 0 || _seconds(t, b) < 0),
                              (:useable_stop, "USEABLE_STOP_TIME", t -> _seconds(a, t) < 0 || _seconds(t, b) < 0))
        t = getfield(seg, field)
        t === nothing && continue
        if bad(t)
            _fixed!(d, get(c.meta, key, ("", line))[2], key, "outside the data; ignored")
            setfield!(seg, field, nothing)
        end
    end
    if seg.useable_start !== nothing && seg.useable_stop !== nothing && _seconds(seg.useable_start, seg.useable_stop) < 0
        _fixed!(d, line, "USEABLE_STOP_TIME", "before USEABLE_START_TIME; both ignored")
        seg.useable_start = seg.useable_stop = nothing
    end
    n = length(seg.times)
    m = seg.interpolation
    need = required_points(m === nothing ? _DEFAULT_INTERPOLATION : m)
    if n < need
        if m === nothing
            seg.interpolation = Lagrange(max(1, n - 1))
            n == 1 ? _error!(d, line, "data", "the segment starting at line $line has one ephemeris line; nothing can be interpolated") :
                     _note!(d, line, "INTERPOLATION", "none declared and the segment has $n points; Lagrange($(n - 1)) used")
        else
            where_ = isempty(c.skipped) ? "the segment starting at line $line" :
                     "after skipping line(s) $(join(c.skipped, ", ")), a piece of the segment starting at line $line"
            _error!(d, line, "INTERPOLATION", "$where_ has $n points; $(m) needs $need (§5.2.4.7)")
        end
    end
end

function _oem_check_segments!(p::_OEMParser)
    segs = p.file.segments
    for k in 2:length(segs)
        prev, cur = segs[k - 1], segs[k]
        pe = something(prev.useable_stop, prev.stop_time)
        cs = something(cur.useable_start, cur.start_time)
        if _seconds(pe, cs) < 0
            _fixed!(p.d, 0, "segments", "segments $(k - 1) and $k overlap in their usable spans; " *
                    "§5.2.4.4 allows only a shared endpoint. Where both cover a time, choose with overlap = :earlier or :later.")
        end
    end
end

# ============================================================================================
# Writing
# ============================================================================================

"""
    write_oem(path, data::EphemerisFileData; digits = 9)

Write `data` as a CCSDS Orbit Ephemeris Message in keyword = value form, version
`data.version` (`"3.0"`, or `"2.0"` for older tools).

# Arguments
- `path`: The file to write. It is replaced if it exists.
- `data`: The header and segments to write.
- `digits`: Decimal places of seconds in epochs, 0 to 9. Default 9.

# Notes
- Everything is checked before the file is opened, and every problem is reported in one
  `ArgumentError`. A file that fails a check is not written, and an existing file at `path` is
  left as it was.
- Checked: the version; mandatory header and metadata values; frames, centres and time systems
  a reader can map back; one time system for the file (§5.2.4.5); epochs strictly increasing;
  finite values; accelerations on every line of a segment or none; enough points for the
  declared interpolation (§5.2.4.7); declared and usable spans consistent with the data, and
  usable spans of consecutive segments not overlapping (§5.2.4.4); covariance epochs increasing,
  values finite and symmetric; printable ASCII and lines of at most 254 characters (§7.3).
- Numbers are written with 16 significant digits, the most §7.5.7 allows. A value read back
  agrees to within one part in 10¹⁶.
- CREATION_DATE is the time of writing, in UTC, unless `data.creation_date` is set.

# Examples
```jldoctest
using AstroEpochs, StaticArrays
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
times = [Time(t0.jd1, t0.jd2 + k / 1440, :utc, :jd) for k in 0:7]
states = [SVector(7000.0, 450.0k, 0.0, 0.0, 7.5, 0.0) for k in 0:7]
seg = EpicycleIO.EphemerisSegmentData(object_name = "SAT", object_id = "2026-001A",
          center = "EARTH", frame = "EME2000", time_system = "UTC", times = times, states = states)
data = EpicycleIO.EphemerisFileData(originator = "EXAMPLE", segments = [seg])
path = tempname() * ".oem"
EpicycleIO.write_oem(path, data)
readline(path)

# output

"CCSDS_OEM_VERS       = 3.0"
```
"""
function write_oem(path::AbstractString, data::EphemerisFileData; digits::Integer = 9)
    0 <= digits <= 9 || throw(ArgumentError("write_oem: digits must be 0 to 9; got $digits"))
    lines = _oem_lines(data, digits)                         # throws, listing every problem
    dir = dirname(abspath(path))
    isdir(dir) || throw(ArgumentError("write_oem: folder $dir does not exist"))
    tmp, io = mktemp(dir)
    try
        for l in lines
            write(io, l, '\n')
        end
        close(io)
        mv(tmp, path; force = true)
    catch
        isopen(io) && close(io)
        rm(tmp; force = true)
        rethrow()
    end
    return path
end

_kv(key, value) = rpad(key, _KEY_WIDTH) * " = " * value
_num(x) = @sprintf("%.15e", x)

function _check_text!(problems, what, s; mandatory = false)
    mandatory && isempty(strip(s)) && push!(problems, "$what is empty; it is mandatory")
    all(c -> isascii(c) && isprint(c), s) || push!(problems, "$what has characters other than printable ASCII (§7.3.4)")
end

function _oem_lines(data::EphemerisFileData, digits)
    problems = String[]
    v = data.version
    v in _OEM_WRITE_VERSIONS || push!(problems, "version $(repr(v)) cannot be written; use \"3.0\" or \"2.0\"")
    if v == "2.0"
        isempty(data.message_id) || push!(problems, "MESSAGE_ID is an OEM 3.0 keyword; it cannot be written in version 2.0")
        isempty(data.classification) || push!(problems, "CLASSIFICATION is an OEM 3.0 keyword; it cannot be written in version 2.0")
    end
    _check_text!(problems, "ORIGINATOR", data.originator; mandatory = true)
    _check_text!(problems, "MESSAGE_ID", data.message_id)
    _check_text!(problems, "CLASSIFICATION", data.classification)
    foreach(c -> _check_text!(problems, "a header comment", c), data.comments)
    isempty(data.segments) && push!(problems, "there are no segments to write")
    systems = unique(s.time_system for s in data.segments)
    length(systems) > 1 && push!(problems, "segments use time systems $(join(systems, ", ")); " *
                                 "§5.2.4.5 requires one per file")
    for (k, s) in enumerate(data.segments)
        _check_oem_segment!(problems, k, s)
    end
    for k in 2:length(data.segments)
        a, b = data.segments[k - 1], data.segments[k]
        (isempty(a.times) || isempty(b.times)) && continue
        ae = something(a.useable_stop, a.stop_time, last(a.times))
        bs = something(b.useable_start, b.start_time, first(b.times))
        _seconds(ae, bs) < 0 && push!(problems, "segments $(k - 1) and $k overlap in their usable spans; " *
            "§5.2.4.4 allows only a shared endpoint. Set USEABLE_STOP_TIME and USEABLE_START_TIME to separate them.")
    end
    isempty(problems) || throw(ArgumentError("write_oem: nothing written; " *
        "$(length(problems)) problem$(length(problems) == 1 ? "" : "s"):\n  " * join(problems, "\n  ")))

    out = String[]
    push!(out, _kv("CCSDS_OEM_VERS", v))
    foreach(c -> push!(out, "COMMENT " * c), data.comments)
    v == "3.0" && !isempty(data.classification) && push!(out, _kv("CLASSIFICATION", data.classification))
    created = isempty(data.creation_date) ? Dates.format(Dates.now(Dates.UTC), "yyyy-mm-ddTHH:MM:SS") : data.creation_date
    push!(out, _kv("CREATION_DATE", created))
    push!(out, _kv("ORIGINATOR", data.originator))
    v == "3.0" && !isempty(data.message_id) && push!(out, _kv("MESSAGE_ID", data.message_id))
    for s in data.segments
        sys = s.time_system
        ep(t) = _format_epoch(t, sys, digits)
        push!(out, "")
        push!(out, "META_START")
        foreach(c -> push!(out, "COMMENT " * c), s.comments)
        push!(out, _kv("OBJECT_NAME", s.object_name))
        push!(out, _kv("OBJECT_ID", s.object_id))
        push!(out, _kv("CENTER_NAME", s.center))
        push!(out, _kv("REF_FRAME", s.frame))
        s.frame_epoch === nothing || push!(out, _kv("REF_FRAME_EPOCH", ep(s.frame_epoch)))
        push!(out, _kv("TIME_SYSTEM", sys))
        push!(out, _kv("START_TIME", ep(something(s.start_time, first(s.times)))))
        s.useable_start === nothing || push!(out, _kv("USEABLE_START_TIME", ep(s.useable_start)))
        s.useable_stop === nothing || push!(out, _kv("USEABLE_STOP_TIME", ep(s.useable_stop)))
        push!(out, _kv("STOP_TIME", ep(something(s.stop_time, last(s.times)))))
        m = s.interpolation
        if m !== nothing
            push!(out, _kv("INTERPOLATION", m isa Hermite ? "HERMITE" : (m.degree == 1 ? "LINEAR" : "LAGRANGE")))
            push!(out, _kv("INTERPOLATION_DEGREE", string(m.degree)))
        end
        push!(out, "META_STOP")
        push!(out, "")
        for (j, t) in enumerate(s.times)
            vals = s.accelerations === nothing ? s.states[j] : vcat(s.states[j], s.accelerations[j])
            push!(out, ep(t) * " " * join(_num.(vals), " "))
        end
        if !isempty(s.covariances)
            push!(out, "")
            push!(out, "COVARIANCE_START")
            for cv in s.covariances
                push!(out, _kv("EPOCH", ep(cv.epoch)))
                cv.frame == s.frame || push!(out, _kv("COV_REF_FRAME", cv.frame))
                for r in 1:6
                    push!(out, join((_num(cv.matrix[r, c]) for c in 1:r), " "))
                end
            end
            push!(out, "COVARIANCE_STOP")
        end
    end
    long = findfirst(l -> length(l) > _MAX_LINE, out)
    long === nothing || throw(ArgumentError("write_oem: nothing written; line $long would be " *
        "$(length(out[long])) characters, over the $_MAX_LINE §7.3.2 allows: $(first(out[long], 60))…"))
    return out
end

function _check_oem_segment!(problems, k, s::EphemerisSegmentData)
    what = "segment $k"
    _check_text!(problems, "$what OBJECT_NAME", s.object_name; mandatory = true)
    _check_text!(problems, "$what OBJECT_ID", s.object_id; mandatory = true)
    foreach(c -> _check_text!(problems, "a $what comment", c), s.comments)
    body, problem = _center_body(s.center)
    body === nothing && push!(problems, "$what CENTER_NAME: $problem")
    axes, note, problem = _frame_axes(s.frame, s.center)
    axes === nothing && push!(problems, "$what REF_FRAME: $problem")
    sys, note, problem = _time_system(s.time_system)
    (isempty(problem) && sys == s.time_system) ||
        push!(problems, "$what TIME_SYSTEM: $(isempty(problem) ? "write $sys rather than $(s.time_system)" : problem)")
    n = length(s.times)
    n == length(s.states) || push!(problems, "$what has $n epochs but $(length(s.states)) states")
    n == 0 && return push!(problems, "$what has no states")
    isempty(problem) || return                               # epochs cannot be checked without a time system
    for j in 2:n
        if _seconds(s.times[j - 1], s.times[j]) <= 0
            push!(problems, "$what epochs do not increase strictly at state $j")
            break
        end
    end
    bad = findfirst(x -> !all(isfinite, x), s.states)
    bad === nothing || push!(problems, "$what state $bad is not finite")
    if s.accelerations !== nothing
        length(s.accelerations) == n || push!(problems, "$what has $n states but $(length(s.accelerations)) accelerations")
        bad = findfirst(x -> !all(isfinite, x), s.accelerations)
        bad === nothing || push!(problems, "$what acceleration $bad is not finite")
    end
    m = s.interpolation
    m === nothing || n >= required_points(m) ||
        push!(problems, "$what has $n states; $(m) needs $(required_points(m)) (§5.2.4.7)")
    a, b = first(s.times), last(s.times)
    s.start_time === nothing || _seconds(s.start_time, a) >= 0 ||
        push!(problems, "$what START_TIME is after the first epoch")
    s.stop_time === nothing || _seconds(b, s.stop_time) >= 0 ||
        push!(problems, "$what STOP_TIME is before the last epoch")
    for (name, t) in (("USEABLE_START_TIME", s.useable_start), ("USEABLE_STOP_TIME", s.useable_stop))
        t === nothing && continue
        (_seconds(a, t) >= 0 && _seconds(t, b) >= 0) || push!(problems, "$what $name is outside the data")
    end
    if s.useable_start !== nothing && s.useable_stop !== nothing
        _seconds(s.useable_start, s.useable_stop) >= 0 || push!(problems, "$what USEABLE_STOP_TIME is before USEABLE_START_TIME")
    end
    for (j, cv) in enumerate(s.covariances)
        all(isfinite, cv.matrix) || push!(problems, "$what covariance $j is not finite")
        maximum(abs, cv.matrix - cv.matrix') <= 1e-12 * max(1.0, maximum(abs, cv.matrix)) ||
            push!(problems, "$what covariance $j is not symmetric")
        any(i -> cv.matrix[i, i] < 0, 1:6) && push!(problems, "$what covariance $j has a negative variance")
        if !(uppercase(cv.frame) in _COV_RELATIVE_FRAMES)
            ax, _, prob = _frame_axes(cv.frame, s.center)
            ax === nothing && push!(problems, "$what covariance $j frame: $prob")
        end
        j > 1 && _seconds(s.covariances[j - 1].epoch, cv.epoch) <= 0 &&
            push!(problems, "$what covariance epochs do not increase at covariance $j (§5.2.5.7)")
    end
end
