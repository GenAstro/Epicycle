# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# STK ephemeris files (.e), as the Ansys STK help page "Ephemeris (.e) file format" defines them
# (STK 12 and current). Cartesian formats only: EphemerisTimePos, EphemerisTimePosVel and
# EphemerisTimePosVelAcc, with CovarianceTimePos and CovarianceTimePosVel.
#
# Times in the data are seconds from ScenarioEpoch, which is UTC; the seconds are TAI seconds, or
# TDB seconds when `TimeScale TDB` is given. Distances are metres unless DistanceUnit says
# otherwise. Covariance is in metres regardless.

const _STK_WRITE_VERSION = "11.0"
const _STK_MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
const _STK_DATA_KEYS = Dict("EPHEMERISTIMEPOS" => 3, "EPHEMERISTIMEPOSVEL" => 6, "EPHEMERISTIMEPOSVELACC" => 9)
const _STK_COV_KEYS = Dict("COVARIANCETIMEPOS" => 6, "COVARIANCETIMEPOSVEL" => 21)
const _STK_UNSUPPORTED_DATA = ("EPHEMERISLLATIMEPOS", "EPHEMERISLLATIMEPOSVEL", "EPHEMERISMSLLLATIMEPOS",
    "EPHEMERISMSLLLATIMEPOSVEL", "EPHEMERISTERRAINLLATIMEPOS", "EPHEMERISLLRTIMEPOS", "EPHEMERISLLRTIMEPOSVEL",
    "EPHEMERISGEOCENTRICLLATIMEPOS", "EPHEMERISGEOCENTRICLLATIMEPOSVEL", "EPHEMERISGEODETICLLRTIMEPOS",
    "EPHEMERISGEODETICLLRTIMEPOSVEL", "STATEERRORTRANSITION")
const _STK_HEADER_KEYS = ("NUMBEROFEPHEMERISPOINTS", "NUMBEROFCOVARIANCEPOINTS", "SCENARIOEPOCH", "CENTRALBODY",
    "COORDINATESYSTEM", "COORDINATESYSTEMEPOCH", "DISTANCEUNIT", "BLOCKINGFACTOR", "INTERPOLATIONMETHOD",
    "INTERPOLATIONORDER", "INTERPOLATIONSAMPLESM1", "SMOOTHDATA", "COVARIANCEINTERPOLATIONMETHOD",
    "COVARIANCEFORMAT", "TIMEFORMAT", "TIMESCALE", "MESSAGELEVEL", "COMPUTEVELOCITY",
    "STATEERRORTRANSITIONREFERENCEEPOCH", "COVARIANCECOORDINATESYSTEM", "COVARIANCECOORDINATESYSTEMEPOCH",
    "COVARIANCEPROJREFCOORDINATESYSTEM", "COVARIANCEINCLUDEANGLEVELINVELCOVXFM")
# Keywords that change how data are read but that Epicycle does not apply: read, noted, ignored.
const _STK_IGNORED_KEYS = ("BLOCKINGFACTOR", "SMOOTHDATA", "COVARIANCEINTERPOLATIONMETHOD", "MESSAGELEVEL",
    "COMPUTEVELOCITY", "STATEERRORTRANSITIONREFERENCEEPOCH", "COVARIANCEPROJREFCOORDINATESYSTEM",
    "COVARIANCEINCLUDEANGLEVELINVELCOVXFM")
const _STK_DISTANCE = Dict("METERS" => 1e-3, "METER" => 1e-3, "M" => 1e-3,
                           "KILOMETERS" => 1.0, "KILOMETER" => 1.0, "KM" => 1.0)
const _STK_COV_FRAMES = Dict("RIC" => "RIC", "RSW" => "RSW", "RTN" => "RTN", "UVW" => "RIC",
                             "NTC" => "NTC", "TNW" => "TNW", "VNC" => "VNB")

# STK CoordinateSystem ↔ the frame names segments use (the OEM names, see names.jl).
function _stk_frame_to_name(cs::AbstractString, center::AbstractString)
    c = lowercase(cs)
    c == "j2000" && return "EME2000", ""
    c == "icrf" && return "ICRF", ""
    c == "trueofdate" && return center == "EARTH" ? ("TOD", "") :
        ("", "TrueOfDate about $center has no Epicycle axes")
    c == "meanofdate" && return center == "EARTH" ? ("MOD_EARTH", "") :
        ("", "MeanOfDate about $center has no Epicycle axes")
    c == "temeofdate" && return center == "EARTH" ? ("TEME", "") :
        ("", "TEMEOfDate is defined for the Earth only")
    if c == "fixed"
        center == "EARTH" && return "ITRF", ""
        center == "MOON" && return "", "the Moon's Fixed system in STK may be its mean-Earth or principal " *
            "axes; override with ref_frame = \"MOON_ME\" or \"MOON_PA440\" to say which"
        return "FIXED_CB", ""
    end
    c == "inertial" && return "", "STK's Inertial system depends on the central body and STK version; " *
        "override with ref_frame = \"ICRF\" or \"EME2000\" to say which"
    c in ("meanofepoch", "trueofepoch", "temeofepoch", "alignmentatepoch") &&
        return "", "$cs is an epoch-dependent system that Epicycle does not provide"
    (startswith(c, "awb") || c == "custom") &&
        return "", "$cs is a system defined inside an STK scenario; Epicycle cannot know it"
    return "", "$(repr(String(cs))) is not a known STK coordinate system.$(_did_you_mean(cs,
        ("Fixed", "J2000", "ICRF", "Inertial", "TrueOfDate", "MeanOfDate", "TEMEOfDate")))"
end

function _stk_frame_from_name(frame::AbstractString, center::AbstractString)
    f = uppercase(frame)
    f == "EME2000" && return "J2000"
    f == "ICRF" && return "ICRF"
    center == "EARTH" && f == "TOD" && return "TrueOfDate"
    center == "EARTH" && f == "MOD_EARTH" && return "MeanOfDate"
    center == "EARTH" && f == "TEME" && return "TEMEOfDate"
    center == "EARTH" && f == "ITRF" && return "Fixed"
    center != "EARTH" && center != "MOON" && f == "FIXED_CB" && return "Fixed"
    return nothing
end

_stk_body(center::AbstractString) = uppercase(center[1]) * lowercase(center[2:end])

# "1 Jun 2002 12:00:00.000000000" (UTC, Gregorian).
function _stk_parse_date(text::AbstractString)
    m = match(r"^(\d{1,2})\s+([A-Za-z]{3})\s+(\d{4})\s+(\d{1,2}):(\d{2}):(\d{2}(?:\.\d+)?)$", strip(text))
    m === nothing && return nothing, "not an STK date; expected dd mmm yyyy hh:mm:ss.s, e.g. 1 Jun 2026 12:00:00.000"
    mo = findfirst(==(titlecase(lowercase(m[2]))), _STK_MONTHS)
    mo === nothing && return nothing, "month $(repr(String(m[2]))) is not a three-letter English month"
    iso = @sprintf("%s-%02d-%02dT%02d:%s:%s", m[3], mo, parse(Int, m[1]), parse(Int, m[4]), m[5], m[6])
    t, _, problem = _parse_epoch(iso, "UTC")
    return t, problem
end

function _stk_format_date(t::Time, digits::Integer)
    s = _format_epoch(t, "UTC", digits)                    # YYYY-MM-DDThh:mm:ss.f
    y, mo, d = parse(Int, s[1:4]), parse(Int, s[6:7]), parse(Int, s[9:10])
    return @sprintf("%d %s %d %s", d, _STK_MONTHS[mo], y, s[12:end])
end

# ============================================================================================
# Reading
# ============================================================================================

"""
    read_stk(path; strict = false, skip_bad_lines = false,
             ref_frame = nothing, center = nothing) -> EphemerisFileData

Read an STK ephemeris file (`.e`) in a Cartesian format: EphemerisTimePos, EphemerisTimePosVel
or EphemerisTimePosVelAcc, with an optional CovarianceTimePos or CovarianceTimePosVel block.

# Arguments
As for [`read_oem`](@ref). `ref_frame` and `center` take the OEM names (e.g. `"EME2000"`,
`"EARTH"`).

# Notes
- Times are seconds from ScenarioEpoch (UTC), counted in TAI seconds, or TDB seconds with
  `TimeScale TDB`. Positions are converted from DistanceUnit, metres by default, to km.
  Covariance is converted from metres to km.
- SegmentBoundaryTimes split the data into segments; the repeated record at each boundary
  starts the next segment. A repeated time anywhere else is an error unless the two records
  are identical.
- InterpolationSamplesM1 (or the older InterpolationOrder) is the number of points less one:
  Lagrange degree M1, Hermite degree 2·M1 + 1. With neither, STK's default of 5 applies.
- NumberOfEphemerisPoints is the most points STK reads. Fewer rows than it gives is an error, as
  the file may be truncated; extra rows are ignored with a correction recorded, as STK does.
- Formats that are not Cartesian (latitude, longitude and altitude or radius), TimeFormat other
  than seconds, and coordinate systems defined inside an STK scenario are errors that say so.

# Returns
An [`EphemerisFileData`](@ref) with `format = :stk`.

# Examples
```jldoctest
path = tempname() * ".e"
write(path, \"\"\"
stk.v.11.0
BEGIN Ephemeris
ScenarioEpoch 21 Oct 2026 00:00:00.000
CentralBody Earth
CoordinateSystem J2000
DistanceUnit Kilometers
InterpolationMethod Lagrange
InterpolationSamplesM1 2
EphemerisTimePosVel
0.0 7000.0 0.0 0.0 0.0 7.5 0.0
60.0 6998.4 449.9 0.0 -0.05 7.49 0.0
120.0 6993.6 899.3 0.0 -0.1 7.48 0.0
END Ephemeris
\"\"\")
f = EpicycleIO.read_stk(path)
(f.segments[1].frame, length(f.segments[1].times))

# output

("EME2000", 3)
```
"""
function read_stk(path::AbstractString; strict::Bool = false, skip_bad_lines::Bool = false,
                  ref_frame = nothing, center = nothing)
    isfile(path) || throw(ArgumentError("read_stk: no file at $path"))
    text = read(path, String)
    d = _Diagnostics(path; strict)
    if startswith(text, "﻿")
        text = text[nextind(text, 1):end]
        _fixed!(d, 1, "file", "byte-order mark at the start ignored")
    end
    lines = split(text, r"\r\n|\n\r|\r|\n")
    file = EphemerisFileData(format = :stk)
    h = Dict{String, Tuple{String, Int}}()           # header keyword → (value, line)
    boundaries = Tuple{Float64, Int}[]
    rows = Tuple{Float64, Vector{Float64}, Int}[]       # (seconds, values, line)
    cov_rows = Tuple{Float64, Vector{Float64}, Int}[]
    data_kind, cov_kind = "", ""
    state = :version
    pending_cov = Float64[]; pending_cov_line = 0; pending_cov_t = NaN
    for (i, raw) in enumerate(lines)
        s = strip(replace(raw, '\t' => ' '))
        isempty(s) && continue
        startswith(s, "#") && continue                    # comments, as STK files written by GMAT use
        u = uppercase(s)
        if state == :version
            m = match(r"^stk\.v\.(\d+)\.(\d+)$"i, s)
            if m === nothing
                _error!(d, i, "version", "the first line must be the version stamp, e.g. stk.v.11.0; got $(repr(String(first(s, 40))))")
            else
                file.version = "$(m[1]).$(m[2])"
            end
            state = :outside
            m === nothing || continue
        end
        if state == :outside
            if u == "BEGIN EPHEMERIS"
                state = :header
            else
                _error!(d, i, "BEGIN Ephemeris", "expected BEGIN Ephemeris; got $(repr(String(first(s, 40))))")
            end
            continue
        end
        if u == "END EPHEMERIS"
            state = :done
            continue
        end
        state == :done && (_fixed!(d, i, "file", "content after END Ephemeris ignored"); break)
        if u == "BEGIN SEGMENTBOUNDARYTIMES"
            state == :header || _error!(d, i, "SegmentBoundaryTimes", "boundary times after the data; they belong in the header")
            state = :boundaries
            continue
        elseif u == "END SEGMENTBOUNDARYTIMES"
            state = :header
            continue
        elseif u == "BEGIN TRENDINGCONTROL"
            state = :trending
            _note!(d, i, "TrendingControl", "trending control is for STK's smoothing; ignored")
            continue
        elseif u == "END TRENDINGCONTROL"
            state = :header
            continue
        end
        state == :trending && continue
        if state == :boundaries
            x = tryparse(Float64, s)
            x === nothing ? _error!(d, i, "SegmentBoundaryTimes", "$(repr(String(s))) is not a time in seconds") :
                            push!(boundaries, (x, i))
            continue
        end
        if haskey(_STK_DATA_KEYS, u)
            isempty(data_kind) || _error!(d, i, s, "a second ephemeris block; a file holds one")
            data_kind, state = u, :data
            continue
        elseif haskey(_STK_COV_KEYS, u)
            isempty(cov_kind) || _error!(d, i, s, "a second covariance block; a file holds one")
            cov_kind, state = u, :cov
            continue
        elseif u in _STK_UNSUPPORTED_DATA
            _error!(d, i, s, "$s is not a Cartesian format; Epicycle reads EphemerisTimePos, " *
                    "EphemerisTimePosVel and EphemerisTimePosVelAcc")
            state = :skip
            continue
        end
        state == :skip && continue
        if state == :header
            parts = split(s; limit = 2)
            key = uppercase(parts[1])
            value = length(parts) == 2 ? String(strip(parts[2])) : ""
            if !(key in _STK_HEADER_KEYS)
                _fixed!(d, i, parts[1], "not an STK ephemeris keyword; ignored.$(_did_you_mean(key, _STK_HEADER_KEYS))")
            elseif haskey(h, key)
                _error!(d, i, parts[1], "given twice")
            else
                h[key] = (value, i)
            end
            continue
        end
        toks = split(s)
        if state == :data
            n = _STK_DATA_KEYS[data_kind] + 1
            vals = Float64[]
            ok = length(toks) == n
            ok || _bad_stk!(d, i, skip_bad_lines, "data line has $(length(toks)) values; $data_kind has $n")
            if ok
                for tok in toks
                    x = tryparse(Float64, tok)
                    if x === nothing || !isfinite(x)
                        _bad_stk!(d, i, skip_bad_lines, "$(repr(String(tok))) is not a finite number")
                        ok = false
                        break
                    end
                    push!(vals, x)
                end
            end
            ok && push!(rows, (vals[1], vals[2:end], i))
        elseif state == :cov
            # A covariance row may wrap over several lines, as STK's own help shows it.
            n = _STK_COV_KEYS[cov_kind]
            for tok in toks
                x = tryparse(Float64, tok)
                x === nothing && (_error!(d, i, "covariance", "$(repr(String(tok))) is not a number"); break)
                if isnan(pending_cov_t)
                    pending_cov_t, pending_cov_line = x, i
                else
                    push!(pending_cov, x)
                    if length(pending_cov) == n
                        push!(cov_rows, (pending_cov_t, pending_cov, pending_cov_line))
                        pending_cov, pending_cov_t = Float64[], NaN
                    end
                end
            end
        end
    end
    state == :version && _error!(d, 1, "file", "empty file")
    state in (:header, :data, :cov, :boundaries) && _error!(d, length(lines), "END Ephemeris", "END Ephemeris missing; the file may be truncated")
    isnan(pending_cov_t) || _error!(d, pending_cov_line, "covariance", "covariance row cut short")
    isempty(data_kind) && !_has_errors(d) && _error!(d, length(lines), "data", "no ephemeris data block")

    seg = _stk_segment_meta!(d, h, ref_frame, center, file)
    if seg !== nothing && !isempty(rows)
        _stk_build!(d, file, seg, h, data_kind, rows, boundaries, cov_kind, cov_rows)
    end
    hint = any(x -> x.severity == :error && x.field in ("CoordinateSystem", "CentralBody"), d.items) ?
        "Fix the file, or override: read_ephemeris(path; ref_frame = \"…\", center = \"…\")" : ""
    file.diagnostics = _finish(d; hint)
    return file
end

function _bad_stk!(d, i, skip, msg)
    skip ? _fixed!(d, i, "data", "$msg; line skipped and the segment split here") : _error!(d, i, "data", msg)
end

function _stk_segment_meta!(d, h, ref_frame, center, file)
    seg = EphemerisSegmentData()
    ok = true
    # Time system of the data: UTC epoch, TAI or TDB seconds.
    if haskey(h, "TIMEFORMAT") && !(uppercase(h["TIMEFORMAT"][1]) in ("EPSEC", "SEC"))
        _error!(d, h["TIMEFORMAT"][2], "TimeFormat", "TimeFormat $(h["TIMEFORMAT"][1]) is not supported; Epicycle reads times in seconds from ScenarioEpoch")
        ok = false
    end
    seg.time_system = "UTC"
    if haskey(h, "TIMESCALE")
        v, line = h["TIMESCALE"]
        uppercase(v) == "TDB" ? (seg.time_system = "TDB") :
            _error!(d, line, "TimeScale", "TimeScale $(repr(v)) is not TDB, the only value STK defines")
    end
    if !haskey(h, "SCENARIOEPOCH")
        _error!(d, 0, "ScenarioEpoch", "ScenarioEpoch missing; without it the times have no origin (STK would use its scenario's epoch)")
        ok = false
    else
        t, problem = _stk_parse_date(h["SCENARIOEPOCH"][1])
        t === nothing ? (_error!(d, h["SCENARIOEPOCH"][2], "ScenarioEpoch", problem); ok = false) : (seg.start_time = t)
    end
    # Central body.
    cb = haskey(h, "CENTRALBODY") ? h["CENTRALBODY"] : ("Earth", 0)
    haskey(h, "CENTRALBODY") || _note!(d, 0, "CentralBody", "none given; Earth, STK's default, used")
    cname = uppercase(cb[1])
    if center !== nothing
        _note!(d, cb[2], "CentralBody", "file says $(repr(cb[1])), read as $(center) by override")
        cname = uppercase(String(center))
    end
    body, problem = _center_body(cname)
    body === nothing ? (_error!(d, cb[2], "CentralBody", problem); ok = false) : (seg.center = cname)
    # Coordinate system, as an OEM name.
    cs = haskey(h, "COORDINATESYSTEM") ? h["COORDINATESYSTEM"] : ("Fixed", 0)
    haskey(h, "COORDINATESYSTEM") || _fixed!(d, 0, "CoordinateSystem",
        "none given; Fixed, STK's default, used. STK recommends always giving one")
    if ref_frame !== nothing
        _note!(d, cs[2], "CoordinateSystem", "file says $(repr(cs[1])), read as $(ref_frame) by override")
        axes, note, problem = _frame_axes(String(ref_frame), seg.center)
        axes === nothing ? (_error!(d, cs[2], "CoordinateSystem", problem); ok = false) : (seg.frame = uppercase(String(ref_frame)))
    elseif !isempty(seg.center)
        name, problem = _stk_frame_to_name(cs[1], seg.center)
        isempty(name) ? (_error!(d, cs[2], "CoordinateSystem", problem); ok = false) : (seg.frame = name)
    end
    # Interpolation.
    m1line = 0
    m1 = 5
    for key in ("INTERPOLATIONSAMPLESM1", "INTERPOLATIONORDER")
        haskey(h, key) || continue
        v, m1line = h[key]
        x = tryparse(Int, v)
        if x === nothing || x < 1
            _fixed!(d, m1line, key, "$(repr(v)) is not a positive integer; STK's default of 5 used")
        else
            m1 = x
        end
        key == "INTERPOLATIONORDER" && _note!(d, m1line, key, "InterpolationOrder is the older name of InterpolationSamplesM1")
        break
    end
    method = haskey(h, "INTERPOLATIONMETHOD") ? h["INTERPOLATIONMETHOD"] : ("Lagrange", 0)
    mu = uppercase(first(split(method[1])))
    if mu == "HERMITE"
        seg.interpolation = Hermite(2m1 + 1)
    else
        mu == "LAGRANGE" || _note!(d, method[2], "InterpolationMethod",
            "$(method[1]) is not Lagrange or Hermite; the data are read and Lagrange used")
        seg.interpolation = Lagrange(m1)
    end
    for key in _STK_IGNORED_KEYS
        haskey(h, key) && _note!(d, h[key][2], key, "STK-specific setting; ignored")
    end
    seg.object_name = "UNKNOWN"
    return ok ? seg : nothing
end

function _stk_build!(d, file, seg0, h, data_kind, rows, boundaries, cov_kind, cov_rows)
    scale = haskey(h, "DISTANCEUNIT") ? get(_STK_DISTANCE, uppercase(h["DISTANCEUNIT"][1]), nothing) : 1e-3
    if scale === nothing
        _error!(d, h["DISTANCEUNIT"][2], "DistanceUnit", "$(repr(h["DISTANCEUNIT"][1])) is not Meters or Kilometers")
        return
    end
    haskey(h, "DISTANCEUNIT") || _note!(d, 0, "DistanceUnit", "none given; metres, STK's default, used")
    # How many points to read.
    if haskey(h, "NUMBEROFEPHEMERISPOINTS")
        v, line = h["NUMBEROFEPHEMERISPOINTS"]
        n = tryparse(Int, v)
        if n === nothing || n < 1
            _fixed!(d, line, "NumberOfEphemerisPoints", "$(repr(v)) is not a positive integer; all points read")
        elseif n > length(rows)
            _error!(d, line, "NumberOfEphemerisPoints", "says $n points but the file has $(length(rows)); it may be truncated")
            return
        elseif n < length(rows)
            _fixed!(d, line, "NumberOfEphemerisPoints", "says $n points; the $(length(rows) - n) after it are ignored, as STK ignores them")
            resize!(rows, n)
        end
    end
    epoch = seg0.start_time
    base = seg0.time_system == "TDB" ? epoch.tdb : epoch.tai
    at(sec) = Time(base.jd1, base.jd2 + sec / 86400, base.scale, :jd)
    ncols = _STK_DATA_KEYS[data_kind]
    data_kind == "EPHEMERISTIMEPOS" && _note!(d, 0, "EphemerisTimePos",
        "no velocities in the file; they are the derivative of the Lagrange polynomial through the positions")
    seg0.interpolation isa Hermite && ncols == 3 && _error!(d, h["INTERPOLATIONMETHOD"][2], "InterpolationMethod",
        "Hermite needs velocities, which EphemerisTimePos does not have")
    bset = Set(first.(boundaries))
    for k in 2:length(boundaries)
        boundaries[k][1] > boundaries[k - 1][1] || _error!(d, boundaries[k][2], "SegmentBoundaryTimes", "boundary times must ascend")
    end
    # Split rows into segments at the boundaries.
    groups = Vector{Vector{Tuple{Float64, Vector{Float64}, Int}}}([[]])
    for (k, r) in enumerate(rows)
        cur = groups[end]
        if !isempty(cur)
            prev = cur[end]
            if r[1] < prev[1]
                _error!(d, r[3], "data", "time $(r[1]) s is before the previous line's ($(prev[3])); times must ascend")
                continue
            elseif r[1] == prev[1]
                if r[1] in bset
                    push!(groups, [r])                         # the record after the boundary
                    continue
                elseif r[2] == prev[2]
                    _fixed!(d, r[3], "data", "duplicate of line $(prev[3]); one copy kept")
                    continue
                else
                    _error!(d, r[3], "data", "a different state at the same time as line $(prev[3]); " *
                            "STK allows a repeated time only at a SegmentBoundaryTimes entry")
                    continue
                end
            elseif any(b -> prev[1] < b < r[1], bset)
                push!(groups, [])                              # a boundary between points, no repeat
            end
        end
        push!(groups[end], r)
    end
    for (b, line) in boundaries
        any(r -> r[1] == b, rows) || _fixed!(d, line, "SegmentBoundaryTimes", "boundary at $b s matches no data time")
    end
    # Covariance, in metres: km² is m² × 1e-6, for every element.
    covs = CovarianceRecord[]
    if !isempty(cov_kind)
        fmt = haskey(h, "COVARIANCEFORMAT") ? uppercase(h["COVARIANCEFORMAT"][1]) :
              (cov_kind == "COVARIANCETIMEPOS" ? "UT" : "LT")
        fmt = fmt in ("LOWERTRIANGULAR", "LT") ? "LT" : fmt in ("UPPERTRIANGULAR", "UT") ? "UT" : ""
        isempty(fmt) && _error!(d, h["COVARIANCEFORMAT"][2], "CovarianceFormat", "not LowerTriangular or UpperTriangular")
        cframe = seg0.frame
        if haskey(h, "COVARIANCECOORDINATESYSTEM")
            v, line = h["COVARIANCECOORDINATESYSTEM"]
            key = uppercase(first(split(v)))
            if haskey(_STK_COV_FRAMES, key)
                cframe = _STK_COV_FRAMES[key]
            else
                name, problem = _stk_frame_to_name(key, seg0.center)
                isempty(name) ? _error!(d, line, "CovarianceCoordinateSystem", problem) : (cframe = name)
            end
        end
        scale == 1.0 && _note!(d, 0, "DistanceUnit", "covariance read in metres, as STK documents it, though distances are in km")
        dim = cov_kind == "COVARIANCETIMEPOS" ? 3 : 6
        dim == 3 && _note!(d, 0, "CovarianceTimePos", "position covariance only; velocity terms are NaN")
        for (t, vals, line) in cov_rows
            m = fill(NaN, 6, 6)
            k = 0
            for r in 1:dim, c in (fmt == "LT" ? (1:r) : (r:dim))
                k += 1
                m[r, c] = m[c, r] = vals[k] * 1e-6
            end
            any(j -> m[j, j] < 0, 1:dim) && (_error!(d, line, "covariance", "negative variance"); continue)
            push!(covs, CovarianceRecord(at(t), cframe, SMatrix{6, 6}(m)))
        end
    end
    for (k, g) in enumerate(groups)
        isempty(g) && continue
        seg = deepcopy(seg0)
        seg.times = [at(r[1]) for r in g]
        seg.states = [ncols == 3 ? SVector{6}(r[2][1] * scale, r[2][2] * scale, r[2][3] * scale, NaN, NaN, NaN) :
                      SVector{6}((r[2][1:6] .* scale)...) for r in g]
        ncols == 9 && (seg.accelerations = [SVector{3}((r[2][7:9] .* scale)...) for r in g])
        seg.start_time, seg.stop_time = first(seg.times), last(seg.times)
        a, b = first(g)[1], last(g)[1]
        seg.covariances = filter(cv -> a <= _seconds(at(0.0), cv.epoch) <= b, covs)
        n = length(seg.times)
        if n < required_points(seg.interpolation)
            n == 1 && length(groups) > 1 ? _error!(d, g[1][3], "data", "a segment of one point; nothing can be interpolated") :
                _error!(d, g[1][3], "InterpolationSamplesM1", "the segment starting at line $(g[1][3]) has $n points; " *
                        "$(seg.interpolation) needs $(required_points(seg.interpolation))")
        end
        push!(file.segments, seg)
    end
end

# ============================================================================================
# Writing
# ============================================================================================

"""
    write_stk(path, data::EphemerisFileData; digits = 9)

Write `data` as an STK ephemeris file (`.e`), version stamp `stk.v.11.0`.

# Arguments
- `path`: The file to write. It is replaced if it exists.
- `data`: The segments to write, as [`read_stk`](@ref) or [`read_oem`](@ref) returns them.
- `digits`: Decimal places of seconds in ScenarioEpoch, 0 to 9. Default 9.

# Notes
- Everything is checked before the file is opened, and every problem is reported in one
  `ArgumentError`; nothing is written on failure.
- Every keyword STK would otherwise default is written: ScenarioEpoch (the first epoch, UTC),
  CentralBody, CoordinateSystem, DistanceUnit, InterpolationMethod and InterpolationSamplesM1,
  NumberOfEphemerisPoints.
- Distances are in km, unless the data carry covariance: then every distance, and the
  covariance, are in metres, so the file does not depend on whether a reader applies
  DistanceUnit to covariance, which STK's documentation leaves open.
- Consecutive segments must share an endpoint, the time of a maneuver for example: STK
  represents segments by SegmentBoundaryTimes with the time repeated, and has no way to show a gap.
- Times are written as TAI seconds from ScenarioEpoch with 16 significant digits; a time
  read back agrees within 1e-10 s for spans up to a year.
- STK's file has one time system, frame and centre: all segments must agree.

# Examples
```jldoctest
using AstroEpochs, StaticArrays
t0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
times = [Time(t0.jd1, t0.jd2 + k / 1440, :utc, :jd) for k in 0:7]
states = [SVector(7000.0, 450.0k, 0.0, 0.0, 7.5, 0.0) for k in 0:7]
seg = EpicycleIO.EphemerisSegmentData(center = "EARTH", frame = "EME2000", time_system = "UTC",
          times = times, states = states, interpolation = EpicycleIO.Lagrange(5))
path = tempname() * ".e"
EpicycleIO.write_stk(path, EpicycleIO.EphemerisFileData(segments = [seg]))
readline(path)

# output

"stk.v.11.0"
```
"""
function write_stk(path::AbstractString, data::EphemerisFileData; digits::Integer = 9)
    0 <= digits <= 9 || throw(ArgumentError("write_stk: digits must be 0 to 9; got $digits"))
    lines = _stk_lines(data, digits)
    dir = dirname(abspath(path))
    isdir(dir) || throw(ArgumentError("write_stk: folder $dir does not exist"))
    tmp, io = mktemp(dir)
    try
        foreach(l -> write(io, l, '\n'), lines)
        close(io)
        mv(tmp, path; force = true)
    catch
        isopen(io) && close(io)
        rm(tmp; force = true)
        rethrow()
    end
    return path
end

function _stk_lines(data::EphemerisFileData, digits)
    problems = String[]
    segs = data.segments
    isempty(segs) && push!(problems, "there are no segments to write")
    if !isempty(segs)
        s1 = segs[1]
        for (k, s) in enumerate(segs)
            (s.center, s.frame, s.time_system) == (s1.center, s1.frame, s1.time_system) ||
                push!(problems, "segment $k differs from segment 1 in centre, frame or time system; an STK file has one of each")
            (s.interpolation === nothing || typeof(s.interpolation) == typeof(s1.interpolation) &&
             s.interpolation.degree == s1.interpolation.degree) ||
                push!(problems, "segment $k has a different interpolation from segment 1; an STK file has one")
            n = length(s.times)
            n == length(s.states) || push!(problems, "segment $k has $n epochs but $(length(s.states)) states")
            n == 0 && (push!(problems, "segment $k has no states"); continue)
            any(j -> _seconds(s.times[j - 1], s.times[j]) <= 0, 2:n) && push!(problems, "segment $k epochs do not increase strictly")
            any(x -> !all(isfinite, x), s.states) && push!(problems, "segment $k has a state that is not finite")
            (s.accelerations === nothing) == (s1.accelerations === nothing) ||
                push!(problems, "segment $k has accelerations where segment 1 has none, or the reverse")
            s.accelerations === nothing || length(s.accelerations) == n && all(x -> all(isfinite, x), s.accelerations) ||
                push!(problems, "segment $k accelerations are missing or not finite")
            m = something(s.interpolation, Lagrange(5))
            n >= required_points(m) || push!(problems, "segment $k has $n states; $m needs $(required_points(m))")
            for (j, cv) in enumerate(s.covariances)
                all(isfinite, cv.matrix) || push!(problems, "segment $k covariance $j is not finite (STK files carry full position-velocity covariance)")
                maximum(abs, cv.matrix - cv.matrix') <= 1e-12 * max(1.0, maximum(abs, cv.matrix)) ||
                    push!(problems, "segment $k covariance $j is not symmetric")
                cv.frame in (s.frame, "RIC", "RSW", "RTN", "TNW", "NTC") ||
                    push!(problems, "segment $k covariance $j is in $(cv.frame), which an STK file cannot give with the ephemeris in $(s.frame)")
            end
        end
        for k in 2:length(segs)
            (isempty(segs[k - 1].times) || isempty(segs[k].times)) && continue
            abs(_seconds(last(segs[k - 1].times), first(segs[k].times))) <= 1e-9 ||
                push!(problems, "segments $(k - 1) and $k do not share an endpoint; an STK file represents " *
                      "segments by a repeated time at their boundary and cannot show a gap")
        end
        body, prob = _center_body(s1.center)
        body === nothing && push!(problems, "CentralBody: $prob")
        _stk_frame_from_name(s1.frame, s1.center) === nothing &&
            push!(problems, "$(s1.frame) about $(s1.center) has no STK CoordinateSystem that Epicycle writes; " *
                  "supported: EME2000 (J2000), ICRF, and about the Earth ITRF (Fixed), TOD, MOD_EARTH, TEME")
        s1.time_system in ("UTC", "TAI", "TT", "TDB", "TCB", "TCG", "GPS") ||
            push!(problems, "time system $(s1.time_system) is not one Epicycle converts")
        s1.interpolation isa Hermite && any(x -> any(isnan, x), s1.states) &&
            push!(problems, "Hermite needs velocities, and some states have none")
    end
    isempty(problems) || throw(ArgumentError("write_stk: nothing written; $(length(problems)) problem$(length(problems) == 1 ? "" : "s"):\n  " *
                                             join(problems, "\n  ")))
    s1 = segs[1]
    has_cov = any(s -> !isempty(s.covariances), segs)
    scale = has_cov ? 1000.0 : 1.0                             # km → written unit
    epoch_text = _stk_format_date(first(s1.times), digits)
    epoch, _ = _stk_parse_date(epoch_text)                     # times are measured from the epoch as written
    tdb = s1.time_system == "TDB"
    secs(t) = tdb ? (t.tdb - epoch.tdb) * 86400 : (t.tai - epoch.tai) * 86400
    m = something(s1.interpolation, Lagrange(5))
    m1 = m isa Hermite ? (m.degree - 1) ÷ 2 : m.degree
    acc = s1.accelerations !== nothing
    nopos = all(x -> all(isnan, x[4:6]), s1.states)
    npts = sum(length(s.times) for s in segs)
    out = String["stk.v.$_STK_WRITE_VERSION", "BEGIN Ephemeris"]
    kv(k, v) = push!(out, rpad(k, 23) * " " * v)              # at least one space, even for the longest keyword
    kv("NumberOfEphemerisPoints", string(npts))
    kv("ScenarioEpoch", epoch_text)
    tdb && kv("TimeScale", "TDB")
    kv("CentralBody", _stk_body(s1.center))
    kv("CoordinateSystem", _stk_frame_from_name(s1.frame, s1.center))
    kv("DistanceUnit", has_cov ? "Meters" : "Kilometers")
    kv("InterpolationMethod", m isa Hermite ? "Hermite" : "Lagrange")
    kv("InterpolationSamplesM1", string(m1))
    # Covariance keywords are header keywords: they precede every data block.
    covs = [(s, cv) for s in segs for cv in s.covariances]
    if has_cov
        kv("NumberOfCovariancePoints", string(length(covs)))
        kv("CovarianceFormat", "LowerTriangular")
        rel = unique(cv.frame for (_, cv) in covs if cv.frame != s1.frame)
        length(rel) > 1 && throw(ArgumentError("write_stk: nothing written; covariances use frames " *
                                              "$(join(rel, ", ")); an STK file gives one CovarianceCoordinateSystem"))
        if !isempty(rel)
            any(cv -> cv.frame == s1.frame, last.(covs)) && throw(ArgumentError(
                "write_stk: nothing written; some covariances are in $(rel[1]) and some in $(s1.frame); an STK file gives one"))
            kv("CovarianceCoordinateSystem", rel[1] == "RTN" || rel[1] == "RSW" ? "RIC" : rel[1])
        end
    end
    if length(segs) > 1
        push!(out, "BEGIN SegmentBoundaryTimes")
        foreach(k -> push!(out, "    " * _num(secs(first(segs[k].times)))), 2:length(segs))
        push!(out, "END SegmentBoundaryTimes")
    end
    push!(out, "")
    push!(out, nopos ? "EphemerisTimePos" : acc ? "EphemerisTimePosVelAcc" : "EphemerisTimePosVel")
    push!(out, "")
    for (k, s) in enumerate(segs)
        k > 1 && push!(out, "")                                # STK's convention: a blank line between duplicate points
        for (j, t) in enumerate(s.times)
            v = nopos ? s.states[j][1:3] : s.states[j]
            acc && (v = vcat(v, s.accelerations[j]))
            push!(out, _num(secs(t)) * " " * join(_num.(v .* scale), " "))
        end
    end
    if has_cov
        push!(out, "")
        push!(out, "CovarianceTimePosVel")
        for (_, cv) in covs
            vals = [cv.matrix[r, c] * 1e6 for r in 1:6 for c in 1:r]
            push!(out, _num(secs(cv.epoch)) * " " * join(_num.(vals), " "))
        end
    end
    push!(out, "END Ephemeris")
    return out
end
