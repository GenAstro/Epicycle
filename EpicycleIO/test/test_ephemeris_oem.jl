# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# CCSDS OEM reading and writing (ticket #9). Section numbers refer to CCSDS 502.0-B-3.
# Writing is strict: everything is checked before a file is opened. Reading is tolerant but
# conservative: what can be corrected without changing the data is, and recorded; the rest stops
# the read with every problem listed.

using Test
using EpicycleIO
using AstroEpochs
using AstroRoutines: Lagrange, Hermite
import AstroFrames
using StaticArrays

const E = EpicycleIO

# A small valid OEM 3.0 file; `body` replaces the data block, `meta` adds metadata lines.
function oem_text(; version = "3.0", header = "", meta = "", center = "EARTH", frame = "EME2000",
                  system = "UTC", body = nothing, interp = "INTERPOLATION = LAGRANGE\nINTERPOLATION_DEGREE = 3\n")
    data = body === nothing ? join(["2026-10-21T00:0$(k):00 $(7000.0 + k) $(450.0k) 0.0 0.0 7.5 0.$(k)" for k in 0:5], "\n") : body
    """
    CCSDS_OEM_VERS = $version
    $(header)CREATION_DATE = 2026-10-21T00:00:00
    ORIGINATOR = TEST
    META_START
    OBJECT_NAME = SAT
    OBJECT_ID = 2026-001A
    CENTER_NAME = $center
    REF_FRAME = $frame
    TIME_SYSTEM = $system
    START_TIME = 2026-10-21T00:00:00
    STOP_TIME = 2026-10-21T00:05:00
    $(interp)$(meta)META_STOP
    $data
    """
end

function oem_file(text)
    path = tempname() * ".oem"
    write(path, text)
    return path
end

read_quiet(path; kw...) = (@test_logs min_level = Base.CoreLogging.Error E.read_oem(path; kw...))
problems(f) = try f(); String[] catch e e isa E.EphemerisReadError || rethrow(); [sprint(show, p) for p in e.problems] end

utc(s) = Time(s, UTC(), ISOT())

@testset "ephemeris: epochs" begin
    t, notes, prob = E._parse_epoch("2026-10-21T12:34:56.123456789", "UTC")
    @test isempty(prob) && isempty(notes)
    @test E._format_epoch(t, "UTC", 9) == "2026-10-21T12:34:56.123456789"
    # Day of year, the second CCSDS form (§7.5.10).
    t2, _, _ = E._parse_epoch("2026-294T12:34:56.123456789", "UTC")
    @test E._format_epoch(t2, "UTC", 9) == "2026-10-21T12:34:56.123456789"
    @test isempty(E._parse_epoch("2024-366T00:00:00", "UTC")[3])               # leap year
    @test !isempty(E._parse_epoch("2026-366T00:00:00", "UTC")[3])
    # Z is UTC; in another system it is noted and ignored.
    @test isempty(E._parse_epoch("2026-10-21T00:00:00Z", "UTC")[2])
    @test !isempty(E._parse_epoch("2026-10-21T00:00:00Z", "TDB")[2])
    # A leap second reads and writes back as :60, in UTC only.
    ls, _, p = E._parse_epoch("2016-12-31T23:59:60.5", "UTC")
    @test isempty(p) && E._format_epoch(ls, "UTC", 3) == "2016-12-31T23:59:60.500"
    @test !isempty(E._parse_epoch("2016-12-31T23:59:60.5", "TAI")[3])
    # GPS is TAI less 19 s.
    g, _, _ = E._parse_epoch("2026-01-01T00:00:00", "GPS")
    @test E._format_epoch(g, "TAI", 3) == "2026-01-01T00:00:19.000"
    @test E._format_epoch(g, "GPS", 3) == "2026-01-01T00:00:00.000"
    # Rounding carries into the date.
    c, _, _ = E._parse_epoch("2026-12-31T23:59:59.9999999999", "UTC")
    @test E._format_epoch(c, "UTC", 6) == "2027-01-01T00:00:00.000000"
    # Each kind of bad epoch says what is wrong.
    @test occursin("month", E._parse_epoch("2026-13-01T00:00:00", "UTC")[3])
    @test occursin("day", E._parse_epoch("2026-02-30T00:00:00", "UTC")[3])
    @test occursin("hour", E._parse_epoch("2026-01-01T24:00:00", "UTC")[3])
    @test occursin("CCSDS epoch", E._parse_epoch("21/10/2026 00:00", "UTC")[3])
    sp, notes, _ = E._parse_epoch("2026-10-21 00:00:00", "UTC")
    @test sp !== nothing && !isempty(notes)
end

@testset "ephemeris: names" begin
    @test E._frame_axes("EME2000", "EARTH")[1] isa E.MJ2000Eq
    @test E._frame_axes("J2000", "EARTH")[1] isa E.MJ2000Eq
    @test E._frame_axes("tod_earth", "EARTH")[1] isa E.TODEq
    ax, note, _ = E._frame_axes("ITRF2000", "EARTH")
    @test ax isa E.ITRF && occursin("realizations", note)
    _, _, p = E._frame_axes("MOON_PA421", "MOON")
    @test occursin("MOON_PA440", p) && occursin("override", p)
    @test occursin("Did you mean EME2000", E._frame_axes("EME200", "EARTH")[3])
    @test occursin("does not provide", E._frame_axes("B1950", "EARTH")[3])
    @test E._frame_axes("FIXED_CB", "MARS")[1] isa E.CelestialBodyFixed
    @test occursin("barycentre", E._center_body("MARS BARYCENTER")[2])
    @test occursin("Did you mean EARTH", E._center_body("EATH")[2])
    @test E._center_body("moon")[1] === E.moon
    @test E._time_system("et")[1:2] == ("TDB", "ET read as TDB, as SPICE defines it")
    @test occursin("does not convert", E._time_system("UT1")[3])
    for (axes, name) in ((E.ICRF(), "ICRF"), (E.MJ2000Eq(), "EME2000"), (E.ITRF(), "ITRF"),
                         (E.TEME(), "TEME"), (E.TODEq(), "TOD"), (E.MoonPA(), "MOON_PA440"))
        @test E._frame_name(axes, "EARTH") == name
        @test typeof(E._frame_axes(name, "EARTH")[1]) == typeof(axes)       # every name reads back
    end
    @test_throws ArgumentError E._frame_name(AstroFrames.TODEc(), "EARTH")
    @test E._time_system_name(UTC()) == "UTC" && E._time_system_name("gps") == "GPS"
end

@testset "OEM: read a valid file" begin
    f = read_quiet(oem_file(oem_text()))
    @test f.version == "3.0" && f.originator == "TEST" && isempty(f.diagnostics)
    s = only(f.segments)
    @test (s.object_name, s.center, s.frame, s.time_system) == ("SAT", "EARTH", "EME2000", "UTC")
    @test s.interpolation == Lagrange(3)
    @test length(s.times) == 6 && s.states[2] == SVector(7001.0, 450.0, 0.0, 0.0, 7.5, 0.1)
    @test s.accelerations === nothing && isempty(s.covariances)
    # Versions 1.0 and 2.0 read the same way.
    @test read_quiet(oem_file(oem_text(version = "1.0"))).version == "1.0"
    @test read_quiet(oem_file(oem_text(version = "2.0"))).version == "2.0"
    # Windows line endings and blank lines are allowed (§7.3.5, §7.3.7).
    @test length(read_quiet(oem_file(replace(oem_text(), "\n" => "\r\n\r\n"))).segments[1].times) == 6
end

@testset "OEM: accelerations, covariance, segments" begin
    acc = join(["2026-10-21T00:0$(k):00 7000 $(k) 0 0 7.5 0 0.001 0.002 -0.003" for k in 0:5], "\n")
    cov = """
    COVARIANCE_START
    EPOCH = 2026-10-21T00:00:00
    COV_REF_FRAME = RTN
    1.0e-4
    2.0e-6 3.0e-4
    0 0 5.0e-4
    1e-8 0 0 1e-9
    0 1e-8 0 0 1e-9
    0 0 1e-8 0 0 1e-9
    EPOCH = 2026-10-21T00:05:00
    1.0e-4
    0 1.0e-4
    0 0 1.0e-4
    0 0 0 1e-9
    0 0 0 0 1e-9
    0 0 0 0 0 1e-9
    COVARIANCE_STOP
    """
    second = """
    META_START
    OBJECT_NAME = SAT
    OBJECT_ID = 2026-001A
    CENTER_NAME = EARTH
    REF_FRAME = EME2000
    TIME_SYSTEM = UTC
    START_TIME = 2026-10-21T00:05:00
    STOP_TIME = 2026-10-21T00:10:00
    INTERPOLATION = HERMITE
    INTERPOLATION_DEGREE = 5
    META_STOP
    COMMENT after a maneuver
    """ * join(["2026-10-21T00:$(lpad(5 + k, 2, '0')):00 7000 $(k) 0 0 7.6 0" for k in 0:5], "\n") * "\n"
    f = read_quiet(oem_file(oem_text(body = acc) * cov * second))
    a, b = f.segments
    @test a.accelerations[3] == SVector(0.001, 0.002, -0.003)
    @test length(a.covariances) == 2
    c = a.covariances[1]
    @test c.frame == "RTN" && c.matrix[2, 1] == c.matrix[1, 2] == 2.0e-6 && c.matrix[6, 3] == 1e-8
    @test a.covariances[2].frame == "EME2000"                         # REF_FRAME when COV_REF_FRAME absent
    @test b.interpolation == Hermite(5) && b.comments == ["after a maneuver"]
end

@testset "OEM: corrections are made and recorded" begin
    # Each of these reads, with a recorded correction.
    cases = [
        ("lower-case keyword", replace(oem_text(), "ORIGINATOR" => "originator"), "lower case"),
        ("tab", replace(oem_text(), "OBJECT_NAME = SAT" => "OBJECT_NAME =\tSAT"), "tab"),
        ("byte-order mark", "﻿" * oem_text(), "byte-order"),
        ("mixed-case interpolation", replace(oem_text(), "LAGRANGE" => "Lagrange"), "mixed case"),
        ("even Hermite degree", oem_text(interp = "INTERPOLATION = HERMITE\nINTERPOLATION_DEGREE = 4\n"), "degree 5 used"),
        ("Fortran exponent", replace(oem_text(), "7.5 0.0" => "7.5D+00 0.0"), "Fortran"),
        ("start time disagrees", replace(oem_text(), "START_TIME = 2026-10-21T00:00:00" =>
                                                     "START_TIME = 2026-10-20T23:00:00"), "first epoch"),
        ("keywords out of order", replace(oem_text(), "OBJECT_NAME = SAT\nOBJECT_ID = 2026-001A" =>
                                                      "OBJECT_ID = 2026-001A\nOBJECT_NAME = SAT"), "order"),
        ("unknown keyword", oem_text(meta = "OBJECT_COLOR = RED\n"), "not an OEM metadata keyword"),
        ("day-of-year epochs", replace(oem_text(), "2026-10-21T00:0" => "2026-294T00:0"), ""),
    ]
    for (what, text, expect) in cases
        f = E.read_oem(oem_file(text))
        @test length(f.segments[1].times) == 6
        if !isempty(expect)
            @test any(d -> occursin(expect, d.message), f.diagnostics)
            # strict makes the same correction an error.
            @test_throws E.EphemerisReadError E.read_oem(oem_file(text); strict = true)
        end
    end
    # An exact duplicate record is kept once.
    rows = split(oem_text(), "\n")
    i = findfirst(startswith("2026-10-21T00:02"), rows)
    dup = join(vcat(rows[1:i], rows[i], rows[i+1:end]), "\n")
    f = E.read_oem(oem_file(dup))
    @test length(f.segments[1].times) == 6 && any(d -> occursin("duplicate", d.message), f.diagnostics)
    # Overrides are recorded with the file's original value.
    f = E.read_oem(oem_file(oem_text(frame = "EME200")); ref_frame = "EME2000")
    @test f.segments[1].frame == "EME2000"
    @test any(d -> occursin("\"EME200\"", d.message) && occursin("override", d.message), f.diagnostics)
end

@testset "OEM: notes describe a valid file and do not fail strict" begin
    # Valid files that need an assumption to read: no interpolation and few points, a specific
    # ITRF realization, ET, and an override the caller asked for.
    text = oem_text(frame = "ITRF2000", system = "ET", interp = "",
                    body = join(["2026-10-21T00:0$(k):00 7000 $(k) 0 0 7.5 0" for k in 0:2], "
"))
    text = replace(text, "STOP_TIME = 2026-10-21T00:05:00" => "STOP_TIME = 2026-10-21T00:02:00")
    f = @test_logs E.read_oem(oem_file(text); strict = true, center = "EARTH")      # no warning
    @test all(d -> d.severity == :note, f.diagnostics)
    msgs = [d.message for d in f.diagnostics]
    @test any(occursin("realizations", m) for m in msgs) && any(occursin("ET read as TDB", m) for m in msgs)
    @test any(occursin("Lagrange(2) used", m) for m in msgs) && any(occursin("override", m) for m in msgs)
    @test f.segments[1].interpolation == Lagrange(2) && f.segments[1].time_system == "TDB"
end

@testset "OEM: what cannot be corrected stops the read, all at once" begin
    # Two problems in one file are both reported, with line numbers.
    text = replace(oem_text(frame = "EME200"), "2026-10-21T00:03:00 7003.0 1350.0 0.0 0.0 7.5 0.3" =>
                                               "2026-10-21T00:03:00 7003.0 1350.0 0.0 0.0 7.5")
    err = try E.read_oem(oem_file(text)); nothing catch e e end
    @test err isa E.EphemerisReadError
    msg = sprint(showerror, err)
    @test occursin("2 problems", msg) && occursin("EME200", msg) && occursin("Did you mean EME2000", msg)
    # Both are reported, with their lines: the frame, and the short data line below it.
    @test occursin("line 8", msg) && occursin("line 18", msg)
    @test occursin("override", msg)

    bad = Dict(
        "missing time system" => replace(oem_text(), "TIME_SYSTEM = UTC\n" => ""),
        "unsupported time system" => oem_text(system = "UT1"),
        "barycentre" => oem_text(center = "EARTH BARYCENTER"),
        "different states at one epoch" => replace(oem_text(), "2026-10-21T00:02:00 7002.0" =>
                                                               "2026-10-21T00:01:00 7002.0"),
        "too few points for the interpolation" => oem_text(interp = "INTERPOLATION = LAGRANGE\nINTERPOLATION_DEGREE = 9\n"),
        "short covariance" => oem_text() * "COVARIANCE_START\nEPOCH = 2026-10-21T00:00:00\n1.0\n2.0 3.0\nCOVARIANCE_STOP\n",
        "negative variance" => oem_text() * "COVARIANCE_START\nEPOCH = 2026-10-21T00:00:00\n" *
                               "-1.0\n0 1\n0 0 1\n0 0 0 1\n0 0 0 0 1\n0 0 0 0 0 1\nCOVARIANCE_STOP\n",
        "not an OEM" => "hello\nworld\n",
        "empty" => "",
    )
    for (what, text) in bad
        @test !isempty(problems(() -> E.read_oem(oem_file(text))))
    end
    # A bad data line names its line and suggests skip_bad_lines.
    text = replace(oem_text(), "7004.0 1800.0 0.0 0.0 7.5 0.4" => "7004.0 abc 0.0 0.0 7.5 0.4")
    msg = try E.read_oem(oem_file(text)); "" catch e sprint(showerror, e) end
    @test occursin("line 19", msg) && occursin("\"abc\"", msg) && occursin("skip_bad_lines", msg)
end

@testset "OEM: skip_bad_lines splits the segment and rechecks it" begin
    rows = ["2026-10-21T00:$(lpad(k, 2, '0')):00 $(7000.0 + k) 0 0 0 7.5 0" for k in 0:11]
    rows[6] = "2026-10-21T00:05:00 7005.0 0 0 0 7.5"                    # line with a value missing
    text = replace(oem_text(body = join(rows, "\n")), "STOP_TIME = 2026-10-21T00:05:00" => "STOP_TIME = 2026-10-21T00:11:00")
    f = E.read_oem(oem_file(text); skip_bad_lines = true)
    @test length(f.segments) == 2
    @test length(f.segments[1].times) == 5 && length(f.segments[2].times) == 6
    @test any(d -> occursin("segment split", d.message), f.diagnostics)
    # A piece left too short for the declared degree is an error naming the skipped line.
    rows[3] = "2026-10-21T00:02:00 broken"
    text = replace(oem_text(body = join(rows, "\n")), "STOP_TIME = 2026-10-21T00:05:00" => "STOP_TIME = 2026-10-21T00:11:00")
    msg = try E.read_oem(oem_file(text); skip_bad_lines = true); "" catch e sprint(showerror, e) end
    @test occursin("after skipping line", msg) && occursin("needs 4", msg)
end

# A segment as Epicycle would write it: n states a minute apart.
function make_segment(n; acc = false, cov = false, interp = Lagrange(7), t0 = "2026-10-21T00:00:00")
    t = utc(t0)
    times = [Time(t.jd1, t.jd2 + (k * 60.123456789) / 86400, :utc, :jd) for k in 0:n-1]
    w = 2π / 5400
    states = [SVector(7000cos(w * 60k) + 1e-9k, 7000sin(w * 60k), 1234.5678901234567,
                      -7000w * sin(w * 60k), 7000w * cos(w * 60k), 1e-3 / 3) for k in 0:n-1]
    seg = E.EphemerisSegmentData(object_name = "TESTSAT", object_id = "2026-001A", center = "EARTH",
                                 frame = "EME2000", time_system = "UTC", interpolation = interp,
                                 times = times, states = states)
    acc && (seg.accelerations = [SVector(-1e-3, 2e-3, π * 1e-6) for _ in 1:n])
    if cov
        m = SMatrix{6, 6}([i == j ? 1e-3 * i : 1e-7 * (i + j) for i in 1:6, j in 1:6])
        seg.covariances = [E.CovarianceRecord(times[1], "EME2000", m), E.CovarianceRecord(times[end], "RTN", 2m)]
    end
    return seg
end

@testset "OEM: strict write, exact round trip" begin
    data = E.EphemerisFileData(originator = "EPICYCLE", message_id = "TEST-1", classification = "UNCLASSIFIED",
                               comments = ["made by the tests"],
                               segments = [make_segment(20; acc = true, cov = true),
                                           make_segment(10; interp = Hermite(5), t0 = "2026-10-21T00:20:00.234567")])
    data.segments[1].useable_stop = data.segments[1].times[end]
    data.segments[2].useable_start = data.segments[2].times[1]
    path = tempname() * ".oem"
    E.write_oem(path, data)
    lines = readlines(path)
    @test all(l -> length(l) <= 254 && !occursin('\t', l) && all(isascii, l), lines)     # §7.3
    # Every number has a mantissa of at most 16 digits (§7.5.7).
    nums = [m.match for l in lines for m in eachmatch(r"-?\d\.\d+e[+-]\d+", l)]
    @test !isempty(nums) && all(n -> count(isdigit, split(n, 'e')[1]) <= 16, nums)
    back = @test_logs E.read_oem(path; strict = true)                                # no corrections at all
    @test back.message_id == "TEST-1" && back.comments == ["made by the tests"]
    for (a, b) in zip(data.segments, back.segments)
        @test length(a.times) == length(b.times)
        @test maximum(abs(E._seconds(x, y)) for (x, y) in zip(a.times, b.times)) < 1e-9       # ns epochs
        @test all(isapprox.(a.states, b.states; rtol = 1e-15))                                # 16 digits
        @test b.interpolation == a.interpolation
    end
    @test all(isapprox.(data.segments[1].accelerations, back.segments[1].accelerations; rtol = 1e-15))
    c1, c2 = data.segments[1].covariances, back.segments[1].covariances
    @test [c.frame for c in c2] == ["EME2000", "RTN"] && all(isapprox(a.matrix, b.matrix; rtol = 1e-15) for (a, b) in zip(c1, c2))
    # Version 2.0 for older tools: no 3.0-only keywords allowed.
    data.version = "2.0"
    @test_throws ArgumentError E.write_oem(tempname(), data)
    data.message_id = data.classification = ""
    E.write_oem(path, data)
    @test readline(path) == "CCSDS_OEM_VERS       = 2.0" && !any(startswith("MESSAGE_ID"), readlines(path))
    @test E.read_oem(path; strict = true).version == "2.0"
end

@testset "OEM: the writer refuses, lists every problem, and writes nothing" begin
    path = tempname() * ".oem"
    write(path, "existing content")
    function refused(data)
        msg = try E.write_oem(path, data); "" catch e e isa ArgumentError || rethrow(); sprint(showerror, e) end
        @test read(path, String) == "existing content"             # untouched
        return msg
    end
    good() = E.EphemerisFileData(originator = "EPICYCLE", segments = [make_segment(10)])
    d = good(); d.segments[1].times[5] = d.segments[1].times[4]
    @test occursin("increase strictly", refused(d))
    d = good(); d.segments[1].states[3] = SVector(NaN, 0, 0, 0, 0, 0)
    @test occursin("not finite", refused(d))
    d = good(); d.segments[1].interpolation = Lagrange(11)
    @test occursin("needs 12", refused(d))
    d = good(); d.segments[1].frame = "EME200"
    @test occursin("Did you mean EME2000", refused(d))
    d = good(); d.segments[1].center = "MARS BARYCENTER"
    @test occursin("barycentre", refused(d))
    d = good(); push!(d.segments, make_segment(10; t0 = "2026-10-21T00:05:00")); d.segments[2].time_system = "TDB"
    msg = refused(d)
    @test occursin("one per file", msg) && occursin("overlap", msg) && occursin("2 problems", msg)
    d = good(); d.segments[1].covariances = [E.CovarianceRecord(d.segments[1].times[1], "EME2000",
                                             SMatrix{6, 6}([i + 2j for i in 1:6, j in 1:6]))]
    @test occursin("not symmetric", refused(d))
    d = good(); d.originator = ""
    @test occursin("ORIGINATOR is empty", refused(d))
    d = good(); push!(d.comments, "x"^260)
    @test occursin("254", refused(d))
    d = good(); d.version = "1.0"
    @test occursin("\"3.0\" or \"2.0\"", refused(d))
    d = good(); d.segments[1].object_name = "SATé"
    @test occursin("printable ASCII", refused(d))
end

# GMAT's sample OEMs, when GMAT is installed on this machine: every one reads.
const GMAT_OEM = joinpath(homedir(), "Downloads", "gmat-internal-win-R2026a-rc1", "data", "vehicle", "ephem", "ccsds")
if isdir(GMAT_OEM)
    @testset "OEM: GMAT's sample files" begin
        for f in filter(endswith(".oem"), readdir(GMAT_OEM; join = true))
            d = E.read_oem(f)
            @test !isempty(d.segments) && all(s -> length(s.times) > 8, d.segments)
            # A value read equals the text in the file.
            row = split(first(filter(l -> occursin(r"^\d{4}-", l), readlines(f))))
            @test d.segments[1].states[1][1] == parse(Float64, row[2])
        end
    end
end
