# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# STK ephemeris (.e) reading and writing (ticket #9), against the STK help page "Ephemeris (.e)
# file format". Also OEM ↔ STK, since the two share one data model.

using Test
using EpicycleIO
using AstroEpochs
using AstroRoutines: Lagrange, Hermite
using StaticArrays

const ES = EpicycleIO

function stk_text(; header = "", body = nothing, kind = "EphemerisTimePosVel", extra = "",
                  frame = "J2000", unit = "DistanceUnit Kilometers\n")
    rows = body === nothing ?
        join(["$(60.0k) $(7000.0 + k) $(450.0k) 0.0 0.0 7.5 0.$(k)" for k in 0:7], "\n") : body
    """
    stk.v.11.0
    # written by the tests
    BEGIN Ephemeris
    NumberOfEphemerisPoints $(length(split(rows, "\n"; keepempty = false)))
    ScenarioEpoch 21 Oct 2026 00:00:00.000
    CentralBody Earth
    CoordinateSystem $frame
    $(unit)InterpolationMethod Lagrange
    InterpolationSamplesM1 3
    $(header)$kind
    $rows
    $(extra)END Ephemeris
    """
end
stk_file(text) = (p = tempname() * ".e"; write(p, text); p)
stk_utc(s) = Time(s, UTC(), ISOT())

@testset "STK: dates and frames" begin
    t, p = ES._stk_parse_date("1 Jun 2002 12:00:00.000000000")
    @test isempty(p) && ES._format_epoch(t, "UTC", 3) == "2002-06-01T12:00:00.000"
    @test ES._stk_parse_date("01 Jan 2000 11:59:28.000")[1] !== nothing
    @test occursin("month", ES._stk_parse_date("1 Foo 2002 12:00:00")[2])
    @test ES._stk_format_date(stk_utc("2026-10-21T01:02:03.25"), 3) == "21 Oct 2026 01:02:03.250"
    @test ES._stk_frame_to_name("J2000", "EARTH") == ("EME2000", "")
    @test ES._stk_frame_to_name("fixed", "EARTH") == ("ITRF", "")
    @test ES._stk_frame_to_name("Fixed", "MARS") == ("FIXED_CB", "")
    @test occursin("override", ES._stk_frame_to_name("Fixed", "MOON")[2])
    @test occursin("override", ES._stk_frame_to_name("Inertial", "EARTH")[2])
    @test occursin("inside an STK scenario", ES._stk_frame_to_name("AWB TopoCentric Facility/Site", "EARTH")[2])
    @test occursin("Did you mean J2000", ES._stk_frame_to_name("J200", "EARTH")[2])
    for (name, stk) in (("EME2000", "J2000"), ("ICRF", "ICRF"), ("ITRF", "Fixed"), ("TOD", "TrueOfDate"),
                        ("MOD_EARTH", "MeanOfDate"), ("TEME", "TEMEOfDate"))
        @test ES._stk_frame_from_name(name, "EARTH") == stk
        @test ES._stk_frame_to_name(stk, "EARTH")[1] == name                 # and back
    end
end

@testset "STK: read" begin
    f = ES.read_stk(stk_file(stk_text()))
    @test isempty(f.diagnostics) && f.version == "11.0"
    s = only(f.segments)
    @test (s.center, s.frame, s.time_system) == ("EARTH", "EME2000", "UTC")
    @test s.interpolation == Lagrange(3) && length(s.times) == 8
    @test s.states[2] == SVector(7001.0, 450.0, 0.0, 0.0, 7.5, 0.1)
    # Times are TAI seconds from the UTC ScenarioEpoch.
    @test ES._seconds(stk_utc("2026-10-21T00:00:00"), s.times[3]) ≈ 120.0 atol = 1e-9
    # Meters by default.
    m = ES.read_stk(stk_file(stk_text(unit = "")))
    @test m.segments[1].states[2][1] ≈ 7.001 && any(d -> occursin("meters", d.message), m.diagnostics)
    # Hermite: degree 2·M1 + 1.
    h = ES.read_stk(stk_file(replace(stk_text(), "InterpolationMethod Lagrange" => "InterpolationMethod Hermite")))
    @test h.segments[1].interpolation == Hermite(7)
    # The older InterpolationOrder means the same, and keywords are not case sensitive.
    o = ES.read_stk(stk_file(replace(stk_text(), "InterpolationSamplesM1" => "interpolationorder")))
    @test o.segments[1].interpolation == Lagrange(3)
    # STK's default of 5 without either.
    @test ES.read_stk(stk_file(replace(stk_text(), "InterpolationSamplesM1 3\n" => ""))).segments[1].interpolation == Lagrange(5)
    # Accelerations and positions only.
    acc = join(["$(60.0k) 7000 $k 0 0 7.5 0 1e-3 2e-3 3e-3" for k in 0:7], "\n")
    a = ES.read_stk(stk_file(stk_text(body = acc, kind = "EphemerisTimePosVelAcc")))
    @test a.segments[1].accelerations[1] == SVector(1e-3, 2e-3, 3e-3)
    pos = join(["$(60.0k) 7000 $k 0" for k in 0:7], "\n")
    p = ES.read_stk(stk_file(stk_text(body = pos, kind = "EphemerisTimePos")))
    @test isnan(p.segments[1].states[1][4]) && any(d -> occursin("no velocities", d.message), p.diagnostics)
    # A positions-only file writes back as one, positions unchanged.
    path = tempname() * ".e"
    ES.write_stk(path, p)
    @test any(==("EphemerisTimePos"), readlines(path))
    back = ES.read_stk(path)
    @test [x[1:3] for x in back.segments[1].states] == [x[1:3] for x in p.segments[1].states]
    @test all(x -> all(isnan, x[4:6]), back.segments[1].states)
    # TimeScale TDB: TDB seconds from the UTC epoch.
    tdb = ES.read_stk(stk_file(stk_text(header = "TimeScale TDB\n")))
    @test tdb.segments[1].time_system == "TDB"
    @test (tdb.segments[1].times[3] - stk_utc("2026-10-21T00:00:00").tdb) * 86400 ≈ 120.0 atol = 1e-9
end

@testset "STK: segment boundaries" begin
    rows = ["$(60.0k) 7000 $k 0 0 7.5 0" for k in 0:9]
    insert!(rows, 7, "300.0 7000 5 0 0 7.6 0")                    # after the maneuver at 300 s
    insert!(rows, 7, "")                                          # STK's blank line between the two
    text = stk_text(header = "BEGIN SegmentBoundaryTimes\n    3.000000000000000e+02\nEND SegmentBoundaryTimes\n",
                    body = join(rows, "\n"))
    f = ES.read_stk(stk_file(text))
    a, b = f.segments
    @test length(a.times) == 6 && length(b.times) == 5
    @test a.states[end][5] == 7.5 && b.states[1][5] == 7.6               # before and after the maneuver
    @test ES._seconds(a.times[end], b.times[1]) == 0
    # A repeated time that is not a boundary, with different states, is an error.
    msg = try ES.read_stk(stk_file(stk_text(body = join(filter(!isempty, rows), "\n")))); "" catch e sprint(showerror, e) end
    @test occursin("SegmentBoundaryTimes", msg)
    # A boundary on a data time that is not repeated still splits there: the point ends one
    # segment and starts the next, and the correction is recorded.
    single = stk_text(header = "BEGIN SegmentBoundaryTimes\n    3.000000000000000e+02\nEND SegmentBoundaryTimes\n",
                      body = join(["$(60.0k) 7000 $k 0 0 7.5 0" for k in 0:9], "\n"))
    f = ES.read_stk(stk_file(single))
    @test [length(s.times) for s in f.segments] == [6, 5]
    @test f.segments[1].states[end] == f.segments[2].states[1]
    @test any(d -> occursin("one data line", d.message), f.diagnostics)
    @test_throws ES.EphemerisReadError ES.read_stk(stk_file(single); strict = true)
end

# Each problem the reader cannot correct, reported on its line, under its keyword, saying what is
# wrong. In stk_text() the header runs from BEGIN Ephemeris on line 3 to InterpolationSamplesM1 on
# line 10, `header` lines go in from line 11, the data block keyword follows them, then eight
# rows, then `extra`, then END Ephemeris.
@testset "STK: each error names its line, its keyword and the problem" begin
    rec(first = "1.0") = "0.0 " * first * " " * join([i == j ? "1.0" : "0.0" for i in 2:6 for j in 1:i], " ")
    covariance = "CovarianceTimePosVel\n" * rec() * "\n"
    boundaries(times...) = "BEGIN SegmentBoundaryTimes\n" * join(times, "\n") * "\nEND SegmentBoundaryTimes\n"
    pos = join(["$(60.0k) 7000 $k 0" for k in 0:7], "\n")
    cases = [
        ("no BEGIN Ephemeris", replace(stk_text(), "BEGIN Ephemeris\n" => "BEGIN Ephem\n"), 3, "BEGIN Ephemeris", "expected BEGIN Ephemeris"),
        ("boundaries after the data", stk_text(extra = boundaries(60)), 20, "SegmentBoundaryTimes", "after the data"),
        ("a boundary that is not a time", stk_text(header = boundaries("abc")), 12, "SegmentBoundaryTimes", "not a time in seconds"),
        ("boundaries out of order", stk_text(header = boundaries(120, 60)), 13, "SegmentBoundaryTimes", "must ascend"),
        ("a second data block", stk_text(extra = "EphemerisTimePosVel\n"), 20, "EphemerisTimePosVel", "a second ephemeris block"),
        ("a second covariance block", stk_text(extra = covariance * covariance), 22, "CovarianceTimePosVel", "a second covariance block"),
        ("a data block that is not Cartesian", stk_text(kind = "EphemerisLLATimePos", body = "0 10 20 300"),
         11, "EphemerisLLATimePos", "not a Cartesian format"),
        ("a keyword twice", stk_text(header = "CentralBody Earth\n"), 11, "CentralBody", "given twice"),
        ("a covariance value that is not a number", stk_text(extra = "CovarianceTimePosVel\n" * rec("abc") * "\n"),
         21, "covariance", "\"abc\" is not a number"),
        ("a covariance row cut short", stk_text(extra = "CovarianceTimePosVel\n0.0 1 2 3\n"), 21, "covariance", "cut short"),
        ("a negative variance", stk_text(extra = "CovarianceTimePosVel\n" * rec("-1.0") * "\n"), 21, "covariance", "negative variance"),
        ("an empty file", "", 1, "file", "empty file"),
        ("no data block", replace(stk_text(), r"EphemerisTimePosVel\n[\s\S]*END Ephemeris" => "END Ephemeris"),
         11, "data", "no ephemeris data block"),
        ("a TimeScale other than TDB", stk_text(header = "TimeScale UTC\n"), 11, "TimeScale", "not TDB"),
        ("a bad ScenarioEpoch", replace(stk_text(), "21 Oct 2026" => "32 Oct 2026"), 5, "ScenarioEpoch", "day 32"),
        ("an unknown CentralBody", replace(stk_text(), "CentralBody Earth" => "CentralBody Erth"), 6, "CentralBody", "\"ERTH\" is not a center"),
        ("an unknown DistanceUnit", stk_text(unit = "DistanceUnit Feet\n"), 8, "DistanceUnit", "not Meters or Kilometers"),
        ("Hermite with positions only", replace(stk_text(kind = "EphemerisTimePos", body = pos), "Method Lagrange" => "Method Hermite"),
         9, "InterpolationMethod", "Hermite needs velocities"),
        ("times out of order", replace(stk_text(), "120.0 7002.0" => "30.0 7002.0"), 14, "data", "before line 13's 60.0 s"),
        ("an unknown CovarianceFormat", stk_text(header = "CovarianceFormat Diagonal\n", extra = covariance),
         11, "CovarianceFormat", "not LowerTriangular or UpperTriangular"),
        ("an unknown CovarianceCoordinateSystem", stk_text(header = "CovarianceCoordinateSystem XYZ\n", extra = covariance),
         11, "CovarianceCoordinateSystem", "XYZ"),
    ]
    for (what, text, line, field, phrase) in cases
        err = try ES.read_stk(stk_file(text)); nothing catch e e end
        found = err isa ES.EphemerisReadError &&
                any(p -> p.line == line && p.field == field && occursin(phrase, p.message), err.problems)
        found || @info "not as expected: $what" err
        @test found
    end
    # An override that is not a frame is reported on the CoordinateSystem line it replaces.
    err = try ES.read_stk(stk_file(stk_text()); ref_frame = "EME200"); nothing catch e e end
    @test any(p -> p.line == 7 && p.field == "CoordinateSystem" && occursin("EME200", p.message), err.problems)
end

@testset "STK: skip_bad_lines drops the line and splits the segment there" begin
    rows = ["$(60.0k) 7000 $k 0 0 7.5 0" for k in 0:15]
    rows[9] = "480.0 7000 abc 0 0 7.5 0"
    text = stk_text(body = join(rows, "\n"))                   # NumberOfEphemerisPoints counts the bad line
    @test_throws ES.EphemerisReadError ES.read_stk(stk_file(text))
    f = ES.read_stk(stk_file(text); skip_bad_lines = true)
    @test [length(s.times) for s in f.segments] == [8, 7]       # no interpolation across the gap
    @test any(d -> occursin("line skipped and the segment split here", d.message), f.diagnostics)
    @test !any(d -> occursin("truncated", d.message), f.diagnostics)
end

@testset "STK: what cannot be read is reported, all at once" begin
    cases = Dict(
        "no version stamp" => replace(stk_text(), "stk.v.11.0\n" => ""),
        "no ScenarioEpoch" => replace(stk_text(), "ScenarioEpoch 21 Oct 2026 00:00:00.000\n" => ""),
        "unknown frame" => stk_text(frame = "J200"),
        "Moon Fixed" => replace(stk_text(frame = "Fixed"), "CentralBody Earth" => "CentralBody Moon"),
        "latitude-longitude" => stk_text(kind = "EphemerisLLATimePos", body = "0 10 20 300"),
        "TimeFormat" => stk_text(header = "TimeFormat JDate\n"),
        "truncated" => replace(stk_text(), "NumberOfEphemerisPoints 8" => "NumberOfEphemerisPoints 9"),
        "no END" => replace(stk_text(), "END Ephemeris" => ""),
        "bad value" => replace(stk_text(), "7003.0" => "x7003"),
        "too few points" => replace(stk_text(), "InterpolationSamplesM1 3" => "InterpolationSamplesM1 9"),
        "digits of another script in the epoch" => replace(stk_text(), "21 Oct 2026" => "21 Oct ٢٠٢٦"),
    )
    for (what, text) in cases
        err = try ES.read_stk(stk_file(text)); nothing catch e e end
        @test err isa ES.EphemerisReadError
    end
    # A sample count no file uses, large enough to overflow Hermite's degree, is corrected.
    huge = replace(stk_text(), "InterpolationMethod Lagrange\nInterpolationSamplesM1 3" =>
                               "InterpolationMethod Hermite\nInterpolationSamplesM1 9223372036854775807")
    f = ES.read_stk(stk_file(huge))
    @test f.segments[1].interpolation == Hermite(11) && any(d -> occursin("from 0 to 100", d.message), f.diagnostics)
    # Two problems, both listed, with the frame's suggestion and the override hint.
    msg = try ES.read_stk(stk_file(replace(stk_text(frame = "J200"), "7003.0" => "x7003"))); "" catch e sprint(showerror, e) end
    @test occursin("2 problems", msg) && occursin("Did you mean J2000", msg) && occursin("override", msg)
    # The override reads it.
    f = ES.read_stk(stk_file(replace(stk_text(frame = "Fixed"), "CentralBody Earth" => "CentralBody Moon"));
                    ref_frame = "MOON_ME")
    @test f.segments[1].frame == "MOON_ME"
    # Extra rows beyond NumberOfEphemerisPoints are ignored, as STK ignores them, and recorded.
    f = ES.read_stk(stk_file(replace(stk_text(), "NumberOfEphemerisPoints 8" => "NumberOfEphemerisPoints 6")))
    @test length(f.segments[1].times) == 6 && any(d -> occursin("ignored", d.message), f.diagnostics)
end

# A segment as Epicycle writes it.
function stk_segment(n; t0 = "2026-10-21T00:00:00", acc = false, interp = Lagrange(7), dv = 0.0)
    t = stk_utc(t0)
    w = 2π / 5400
    seg = ES.EphemerisSegmentData(object_name = "TESTSAT", object_id = "2026-001A", center = "EARTH",
        frame = "EME2000", time_system = "UTC", interpolation = interp,
        times = [Time(t.jd1, t.jd2 + 30.123456789k / 86400, UTC(), JD()) for k in 0:n-1],
        states = [SVector(7000cos(w * 30k), 7000sin(w * 30k), 1234.5678901234567,
                          -7000w * sin(w * 30k), 7000w * cos(w * 30k) + dv, 1e-3 / 3) for k in 0:n-1])
    acc && (seg.accelerations = [SVector(-1e-3, 2e-3, π * 1e-6) for _ in 1:n])
    return seg
end

@testset "STK: strict write and round trip" begin
    a = stk_segment(12)
    b = stk_segment(12; t0 = ES._format_epoch(last(a.times), "UTC", 9), dv = 0.01)   # shares the endpoint
    b.times[1] = last(a.times)
    data = ES.EphemerisFileData(segments = [a, b])
    path = tempname() * ".e"
    ES.write_stk(path, data)
    text = readlines(path)
    @test text[1] == "stk.v.11.0" && text[end] == "END Ephemeris"
    @test any(==("DistanceUnit            Kilometers"), text) && any(==("InterpolationSamplesM1  7"), text)
    @test any(startswith("BEGIN SegmentBoundaryTimes"), text)
    back = @test_logs ES.read_stk(path; strict = true)
    @test length(back.segments) == 2
    for (x, y) in zip(data.segments, back.segments)
        @test length(x.times) == length(y.times)
        @test maximum(abs(ES._seconds(p, q)) for (p, q) in zip(x.times, y.times)) < 1e-9
        @test all(isapprox.(x.states, y.states; rtol = 1e-15))
        @test y.interpolation == Lagrange(7)
    end
    # Endpoints a fraction of a nanosecond apart, which the writer accepts as shared, still
    # write as a repeated time and read back as two segments.
    b.times[1] = Time(last(a.times).jd1, last(a.times).jd2 + 5e-10 / 86400, UTC(), JD())
    ES.write_stk(path, ES.EphemerisFileData(segments = [a, b]))
    @test length(ES.read_stk(path; strict = true).segments) == 2
    # Hermite 7 is written as 3 samples less one, and reads back as Hermite 7.
    h = stk_segment(12; interp = Hermite(7))
    ES.write_stk(path, ES.EphemerisFileData(segments = [h]))
    @test any(==("InterpolationSamplesM1  3"), readlines(path))
    @test ES.read_stk(path; strict = true).segments[1].interpolation == Hermite(7)
    # With covariance, everything is in meters and comes back in km.
    c = stk_segment(12; acc = true)
    m = SMatrix{6, 6}([i == j ? 1e-3 * i : 1e-7 * (i + j) for i in 1:6, j in 1:6])
    c.covariances = [ES.CovarianceRecord(c.times[1], "EME2000", m), ES.CovarianceRecord(c.times[6], "EME2000", 2m)]
    ES.write_stk(path, ES.EphemerisFileData(segments = [c]))
    @test any(==("DistanceUnit            Meters"), readlines(path))
    # Every keyword is separated from its value, the longest (NumberOfCovariancePoints) included.
    hdr = filter(l -> occursin(r"^[A-Z][A-Za-z0-9]+\s", l) || occursin(r"^[A-Z][A-Za-z]+\d", l), readlines(path))
    @test all(l -> occursin(r"^[A-Za-z0-9]+\s+\S", l) && !occursin(r"^[A-Za-z]+Points\d", l), hdr)
    back = ES.read_stk(path; strict = true).segments[1]
    @test all(isapprox.(c.states, back.states; rtol = 1e-15))
    @test all(isapprox.(c.accelerations, back.accelerations; rtol = 1e-15))
    @test length(back.covariances) == length(c.covariances)
    @test all(abs(ES._seconds(x.epoch, y.epoch)) < 1e-9 && x.frame == y.frame for (x, y) in zip(c.covariances, back.covariances))
    @test all(isapprox(x.matrix, y.matrix; rtol = 1e-14) for (x, y) in zip(c.covariances, back.covariances))
    # A TDB segment is written with TimeScale TDB, its times TDB seconds from the epoch, and reads
    # back at the same instants.
    t = stk_segment(12); t.time_system = "TDB"
    ES.write_stk(path, ES.EphemerisFileData(segments = [t]))
    @test any(l -> occursin(r"^TimeScale\s+TDB$", l), readlines(path))
    back_tdb = only(ES.read_stk(path; strict = true).segments)
    @test back_tdb.time_system == "TDB" && maximum(abs(ES._seconds(a, b)) for (a, b) in zip(t.times, back_tdb.times)) < 1e-9
    # Hermite degree 1 is written as InterpolationSamplesM1 0 and reads back as itself.
    h1 = stk_segment(2; interp = Hermite(1))
    ES.write_stk(path, ES.EphemerisFileData(segments = [h1]))
    @test ES.read_stk(path; strict = true).segments[1].interpolation == Hermite(1)
end

@testset "STK: covariance formats and coordinate systems" begin
    # One symmetric 6×6 covariance in m², written out in each of STK's forms.
    M = [i == j ? 10.0 * i : 0.1 * (i + j) for i in 1:6, j in 1:6]
    lower(m, n) = join([m[r, c] for r in 1:n for c in 1:r], " ")
    upper(m, n) = join([m[r, c] for r in 1:n for c in r:n], " ")
    read_cov(header, kind, values) =
        ES.read_stk(stk_file(stk_text(header = "NumberOfCovariancePoints 1\n" * header,
                                      unit = "DistanceUnit Meters\n", extra = "$kind\n60.0 $values\n")))
    covs(f) = f.segments[1].covariances
    # Lower and upper triangular give the same matrix, in km².
    lt = covs(read_cov("CovarianceFormat LowerTriangular\n", "CovarianceTimePosVel", lower(M, 6)))
    ut = covs(read_cov("CovarianceFormat UpperTriangular\n", "CovarianceTimePosVel", upper(M, 6)))
    @test only(lt).matrix ≈ M * 1e-6 && only(ut).matrix == only(lt).matrix
    # Position only: the velocity terms are NaN, and it says so. Its default format is upper triangular.
    f = read_cov("", "CovarianceTimePos", upper(M, 3))
    p = only(covs(f)).matrix
    @test p[1:3, 1:3] ≈ M[1:3, 1:3] * 1e-6 && all(isnan, p[4:6, :]) && all(isnan, p[:, 4:6])
    @test any(d -> occursin("velocity terms are NaN", d.message), f.diagnostics)
    # Each orbit-relative coordinate system STK names reads as the frame of the same axes.
    for (stk_name, frame) in (("RIC", "RIC"), ("UVW", "RIC"), ("VNC", "VNB"), ("NTC", "NTC"), ("TNW", "TNW"))
        c = covs(read_cov("CovarianceCoordinateSystem $stk_name\n", "CovarianceTimePosVel", lower(M, 6)))
        @test only(c).frame == frame
    end
    # An inertial one names a frame, as CoordinateSystem does.
    @test only(covs(read_cov("CovarianceCoordinateSystem ICRF\n", "CovarianceTimePosVel", lower(M, 6)))).frame == "ICRF"

    # Written in an orbit-relative frame, read back in it: RTN is written as RIC, the same axes,
    # and VNB as STK's VNC.
    path = tempname() * ".e"
    for (frame, keyword, back) in (("RTN", "RIC", "RIC"), ("VNB", "VNC", "VNB"), ("TNW", "TNW", "TNW"))
        s = stk_segment(12)
        s.covariances = [ES.CovarianceRecord(s.times[3], frame, SMatrix{6, 6}(M * 1e-6))]
        ES.write_stk(path, ES.EphemerisFileData(segments = [s]))
        @test any(==("CovarianceCoordinateSystem $keyword"), split.(readlines(path)) .|> l -> join(l, " "))
        c = only(ES.read_stk(path; strict = true).segments[1].covariances)
        @test c.frame == back && c.matrix ≈ M * 1e-6
    end
    # An STK file gives one covariance coordinate system, so mixed frames are refused.
    s = stk_segment(12)
    s.covariances = [ES.CovarianceRecord(s.times[3], "RTN", SMatrix{6, 6}(M * 1e-6)),
                     ES.CovarianceRecord(s.times[5], "TNW", SMatrix{6, 6}(M * 1e-6))]
    @test_throws "one CovarianceCoordinateSystem" ES.write_stk(path, ES.EphemerisFileData(segments = [s]))
    s.covariances[2] = ES.CovarianceRecord(s.times[5], "EME2000", SMatrix{6, 6}(M * 1e-6))
    @test_throws "some covariances are in RTN and some in EME2000" ES.write_stk(path, ES.EphemerisFileData(segments = [s]))
end

@testset "STK: covariance records are counted and none is lost silently" begin
    record(t) = "$t " * join([i == j ? "1.0" : "0.0" for i in 1:6 for j in 1:i], " ")
    with_cov(n, times) = stk_text(header = "NumberOfCovariancePoints $n\n",
                                  extra = "\nCovarianceTimePosVel\n" * join(record.(times), "\n") * "\n")
    # Fewer records than declared is a cut file.
    msg = try ES.read_stk(stk_file(with_cov(2, [0.0])); strict = true); "" catch e sprint(showerror, e) end
    @test occursin("NumberOfCovariancePoints", msg) && occursin("truncated", msg)
    # A record outside the ephemeris is dropped with a correction, which strict refuses.
    f = ES.read_stk(stk_file(with_cov(2, [-10.0, 60.0])))
    @test length(f.segments[1].covariances) == 1
    @test any(d -> occursin("outside the ephemeris", d.message), f.diagnostics)
    @test_throws ES.EphemerisReadError ES.read_stk(stk_file(with_cov(2, [-10.0, 60.0])); strict = true)
    # At a segment boundary a record goes with the segment after it, once.
    rows = ["$(60.0k) 7000 $k 0 0 7.5 0" for k in 0:9]
    insert!(rows, 7, "300.0 7000 5 0 0 7.6 0")
    text = stk_text(header = "BEGIN SegmentBoundaryTimes\n    3.000000000000000e+02\nEND SegmentBoundaryTimes\n" *
                             "NumberOfCovariancePoints 1\n",
                    body = join(rows, "\n"), extra = "\nCovarianceTimePosVel\n" * record(300.0) * "\n")
    a, b = ES.read_stk(stk_file(text); strict = true).segments
    @test isempty(a.covariances) && length(b.covariances) == 1
end

@testset "STK: an empty data block or a missing method is reported, not thrown" begin
    empty_block = replace(stk_text(), r"EphemerisTimePosVel\n[\s\S]*END Ephemeris" => "EphemerisTimePosVel\nEND Ephemeris")
    msg = try ES.read_stk(stk_file(empty_block)); "" catch e sprint(showerror, e) end
    @test occursin("no points", msg)
    f = ES.read_stk(stk_file(replace(stk_text(), "InterpolationMethod Lagrange" => "InterpolationMethod")))
    @test f.segments[1].interpolation == Lagrange(3)
    @test any(d -> occursin("no method given", d.message), f.diagnostics)
end

@testset "STK: the writer refuses, lists every problem, and writes nothing" begin
    path = tempname() * ".e"
    write(path, "existing")
    refused(data) = (msg = try ES.write_stk(path, data); "" catch e sprint(showerror, e) end;
                     @test read(path, String) == "existing"; msg)
    a, b = stk_segment(12), stk_segment(12; t0 = "2026-10-21T01:00:00")
    @test occursin("share an endpoint", refused(ES.EphemerisFileData(segments = [a, b])))
    b = stk_segment(12; t0 = "2026-10-21T01:00:00"); b.frame = "ICRF"
    @test occursin("an STK file has one", refused(ES.EphemerisFileData(segments = [stk_segment(12), b])))
    s = stk_segment(12); s.frame = "MOON_ME"; s.center = "MOON"
    @test occursin("no STK CoordinateSystem", refused(ES.EphemerisFileData(segments = [s])))
    s = stk_segment(12); s.states[2] = SVector(NaN, 0, 0, 0, 0, 0)
    @test occursin("not finite", refused(ES.EphemerisFileData(segments = [s])))
    s = stk_segment(12); s.states[2] = SVector(s.states[2][1:3]..., NaN, NaN, NaN)
    @test occursin("positions only or positions and velocities", refused(ES.EphemerisFileData(segments = [s])))
    s = stk_segment(4)
    @test occursin("needs 8", refused(ES.EphemerisFileData(segments = [s])))
end

@testset "OEM and STK carry the same trajectory" begin
    seg = stk_segment(20; acc = true)
    oem, stk = tempname() * ".oem", tempname() * ".e"
    ES.write_oem(oem, ES.EphemerisFileData(originator = "EPICYCLE", segments = [seg]))
    ES.write_stk(stk, ES.read_oem(oem; strict = true))
    back = ES.read_stk(stk; strict = true)
    oem2 = tempname() * ".oem"
    d = ES.EphemerisFileData(originator = "EPICYCLE", segments = back.segments)
    d.segments[1].object_name, d.segments[1].object_id = seg.object_name, seg.object_id
    ES.write_oem(oem2, d)
    final = ES.read_oem(oem2; strict = true).segments[1]
    @test maximum(abs(ES._seconds(p, q)) for (p, q) in zip(seg.times, final.times)) < 1e-9
    @test all(isapprox.(seg.states, final.states; rtol = 1e-15))
    @test all(isapprox.(seg.accelerations, final.accelerations; rtol = 1e-15))
end

# An stk.v.10.0 file from GMAT R2026a's sample data, trimmed (data/gmat/ASSETS.md): three-digit
# exponents, no interpolation keywords, and its epoch at 11:59:28 UTC.
@testset "STK: GMAT's sample file" begin
    f = ES.read_stk(joinpath(@__DIR__, "data", "gmat", "SampleSTKEphem.e"); strict = true)
    s = only(f.segments)
    @test f.version == "10.0" && length(s.times) == 22 && s.frame == "EME2000" && s.center == "EARTH"
    @test s.interpolation == Lagrange(5)                              # STK's default when the file names none
    @test s.states[1] == SVector(7100.0, 0.0, 1300.0, 0.0, 7.35, 1.0)
    @test ES._format_epoch(s.times[1], "UTC", 3) == "2000-01-01T11:59:28.000"
    @test ES._format_epoch(s.times[end], "UTC", 3) == "2000-01-04T11:59:28.000"
end
