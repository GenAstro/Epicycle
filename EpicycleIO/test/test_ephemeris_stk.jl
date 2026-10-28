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
    # Metres by default.
    m = ES.read_stk(stk_file(stk_text(unit = "")))
    @test m.segments[1].states[2][1] ≈ 7.001 && any(d -> occursin("metres", d.message), m.diagnostics)
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
    )
    for (what, text) in cases
        err = try ES.read_stk(stk_file(text)); nothing catch e e end
        @test err isa ES.EphemerisReadError
    end
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
        times = [Time(t.jd1, t.jd2 + 30.123456789k / 86400, :utc, :jd) for k in 0:n-1],
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
    # Hermite 7 is written as 3 samples less one, and reads back as Hermite 7.
    h = stk_segment(12; interp = Hermite(7))
    ES.write_stk(path, ES.EphemerisFileData(segments = [h]))
    @test any(==("InterpolationSamplesM1  3"), readlines(path))
    @test ES.read_stk(path; strict = true).segments[1].interpolation == Hermite(7)
    # With covariance, everything is in metres and comes back in km.
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
    @test all(isapprox(x.matrix, y.matrix; rtol = 1e-14) for (x, y) in zip(c.covariances, back.covariances))
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

const GMAT_STK = joinpath(homedir(), "Downloads", "gmat-internal-win-R2026a-rc1", "data", "vehicle", "ephem", "stk", "SampleSTKEphem.e")
if isfile(GMAT_STK)
    @testset "STK: GMAT's sample file" begin
        f = ES.read_stk(GMAT_STK)
        s = only(f.segments)
        @test length(s.times) == 2610 && s.frame == "EME2000" && s.center == "EARTH"
        @test s.states[1] == SVector(7100.0, 0.0, 1300.0, 0.0, 7.35, 1.0)
        @test ES._format_epoch(s.times[1], "UTC", 3) == "2000-01-01T11:59:28.000"
    end
end
