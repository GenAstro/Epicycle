# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# JPL SPK reading and writing (ticket #9). CSPICE is the reference: what Epicycle writes, SPICE
# evaluates, and the result must equal Epicycle's own interpolation of the same samples.

using Test
using EpicycleIO
using AstroEpochs, AstroFrames, AstroUniverse, AstroModels, AstroStates
using AstroRoutines: Lagrange, Hermite, interpolate
using StaticArrays
import SPICE

const EK = EpicycleIO

function spk_segment(n; t0 = "2026-10-21T00:00:00", interp = Hermite(7), step = 60.0, frame = "ICRF", name = "TESTSAT")
    t = Time(t0, UTC(), ISOT())
    w = 2π / 5400
    EK.EphemerisSegmentData(object_name = name, object_id = "-900001", center = "EARTH", frame = frame,
        time_system = "UTC", interpolation = interp,
        times = [Time(t.jd1, t.jd2 + step * k / 86400, :utc, :jd) for k in 0:n-1],
        states = [SVector(7000cos(w * step * k), 7000sin(w * step * k), 1000sin(w * step * k / 3),
                          -7000w * sin(w * step * k), 7000w * cos(w * step * k), 1000w / 3 * cos(w * step * k / 3))
                  for k in 0:n-1])
end

# SPICE's own state at `et`, with only `path` loaded for the duration.
function spice_state(path, target, et, center = 399)
    lock(EK._SPICE_LOCK) do
        SPICE.furnsh(path)
        try
            return SPICE.spkgeo(target, et, "J2000", center)[1]
        finally
            SPICE.unload(path)
        end
    end
end

@testset "SPK: ephemeris time converts without loss" begin
    ets = [8.13e8 + 0.123456789k for k in 1:1000]
    @test all(EK._et(EK._from_et(x)) == x for x in ets)
    t = Time("2026-10-21T00:00:00", UTC(), ISOT())
    # An epoch rounds once, to the nearest 64-bit ephemeris time: within half of 0.12 µs here.
    @test abs(EK._seconds(EK._from_et(EK._et(t)), t)) <= eps(EK._et(t)) / 2
end

@testset "SPK: write, read back, and SPICE agrees" begin
    for interp in (Hermite(7), Hermite(5), Lagrange(7), Lagrange(9))
        seg = spk_segment(40; interp)
        path = tempname() * ".bsp"
        EK.write_spk(path, EK.EphemerisFileData(segments = [seg]); naif_id = -900001)
        back = @test_logs EK.read_spk(path)                               # notes only: no warning
        s = only(back.segments)
        @test s.interpolation == interp                                   # type 13 stores the window, not the degree
        @test s.states == seg.states                                      # the stored doubles, exactly
        @test maximum(abs(EK._seconds(a, b)) for (a, b) in zip(seg.times, s.times)) < 1e-6
        @test (s.center, s.frame, s.object_name, s.object_id) == ("EARTH", "ICRF", "TESTSAT", "-900001")
        @test any(d -> occursin("J2000 frame read as ICRF", d.message), back.diagnostics)
        # SPICE's evaluation between samples equals ours, for positions.
        ets = EK._et.(s.times)
        r = [x[1:3] for x in s.states]
        v = [x[4:6] for x in s.states]
        for et in range(ets[1], ets[end]; length = 37)[2:end-1] .+ 0.37
            spice = spice_state(path, -900001, et)
            ours = interp isa Hermite ? interpolate(interp, ets, r, v, et) : interpolate(interp, ets, r, et)
            @test maximum(abs.(spice[1:3] - ours)) < 1e-9                  # km
        end
    end
end

@testset "SPK: coverage, segments and objects" begin
    a = spk_segment(20)
    a.useable_start = a.times[3]                                         # coverage narrower than the samples
    b = spk_segment(20; t0 = "2026-10-21T00:19:00", interp = Lagrange(7), name = "AFTER_BURN")
    path = tempname() * ".bsp"
    EK.write_spk(path, EK.EphemerisFileData(segments = [a, b]); naif_id = -77)
    f = EK.read_spk(path)
    x, y = f.segments
    @test abs(EK._seconds(x.useable_start, a.times[3])) < 1e-6 && x.useable_stop === nothing
    @test y.object_name == "AFTER_BURN" && y.interpolation == Lagrange(7)
    # Two objects in one file: naif_id chooses.
    p2 = tempname() * ".bsp"
    EK.write_spk(p2, EK.EphemerisFileData(segments = [a]); naif_id = -77)
    lock(EK._SPICE_LOCK) do                                              # append a second object
        h = SPICE.spkopa(p2)
        s = spk_segment(10; name = "OTHER")
        SPICE.spkw13(h, -88, 399, "J2000", EK._et(s.times[1]), EK._et(s.times[end]), "OTHER", 7,
                     [collect(v) for v in s.states], EK._et.(s.times))
        SPICE.spkcls(h)
    end
    msg = try EK.read_spk(p2); "" catch e sprint(showerror, e) end
    @test occursin("-77, -88", msg) && occursin("naif_id", msg)
    @test only(EK.read_spk(p2; naif_id = -88).segments).object_name == "OTHER"
    @test occursin("not in the file", try EK.read_spk(p2; naif_id = -5); "" catch e sprint(showerror, e) end)
end

@testset "SPK: a type with no stored samples is sampled through SPICE" begin
    # SPK type 8 (equal-step Lagrange) is not read directly; it goes through SPICE.
    s = spk_segment(60; step = 60.0)
    path = tempname() * ".bsp"
    lock(EK._SPICE_LOCK) do
        h = SPICE.spkopn(path, "TEST", 0)
        SPICE.spkw08(h, -900002, 399, "J2000", EK._et(s.times[1]), EK._et(s.times[end]), "TYPE8", 7,
                     [collect(v) for v in s.states], EK._et(s.times[1]), 60.0)
        SPICE.spkcls(h)
    end
    msg = try EK.read_spk(path); "" catch e sprint(showerror, e) end
    @test occursin("type 8", msg) && occursin("step", msg)
    f = EK.read_spk(path; step = 30.0)
    seg = only(f.segments)
    @test length(seg.times) == 119 && seg.interpolation == Lagrange(7)
    @test EK._seconds(seg.times[end - 1], seg.times[end]) > 29.9         # no near-duplicate at the end
    # The note reports the real interpolation error: recompute it independently against SPICE.
    note = only(filter(d -> occursin("sampled every", d.message), f.diagnostics))
    reported = parse(Float64, match(r"within ([0-9.e+-]+) km", note.message)[1])
    ets = EK._et.(seg.times)
    r = [x[1:3] for x in seg.states]
    actual = maximum(sqrt(sum(abs2, interpolate(Lagrange(7), ets, r, (ets[j] + ets[j + 1]) / 2) -
                                    spice_state(path, -900002, (ets[j] + ets[j + 1]) / 2)[1:3])) for j in 1:length(ets) - 1)
    @test reported ≈ actual rtol = 0.01

    # Reading leaves the caller's kernels as they were: a file the caller loaded stays loaded.
    lock(EK._SPICE_LOCK) do
        SPICE.furnsh(path)
    end
    EK.read_spk(path; step = 30.0)
    loaded = lock(EK._SPICE_LOCK) do
        try SPICE.spkgeo(-900002, EK._et(s.times[10]), "J2000", 399)[1][1] catch; NaN end
    end
    @test loaded ≈ s.states[10][1]
    lock(EK._SPICE_LOCK) do
        SPICE.unload(path)
    end

    # Each segment is sampled from itself, even where a later segment of the file covers the
    # same times and SPICE's own lookup would take that one.
    two = tempname() * ".bsp"
    lock(EK._SPICE_LOCK) do
        h = SPICE.spkopn(two, "TEST", 0)
        for (name, shift) in (("FIRST", 0.0), ("SECOND", 100.0))
            SPICE.spkw08(h, -900002, 399, "J2000", EK._et(s.times[1]), EK._et(s.times[end]), name, 7,
                         [collect(v) .+ [shift, 0, 0, 0, 0, 0] for v in s.states], EK._et(s.times[1]), 60.0)
        end
        SPICE.spkcls(h)
    end
    a, b = EK.read_spk(two; step = 60.0).segments
    @test a.states[5][1] ≈ s.states[5][1] && b.states[5][1] ≈ s.states[5][1] + 100.0
end

@testset "SPK: ECLIPJ2000 and a usable span narrower than the samples" begin
    s = spk_segment(20)
    ets = EK._et.(s.times)
    # An ECLIPJ2000 segment, as another tool would write it, reads as ICRF with SPICE's rotation.
    R = SMatrix{3, 3}(lock(() -> SPICE.pxform("J2000", "ECLIPJ2000", 0.0), EK._SPICE_LOCK))
    eclip = [Vector(vcat(R * x[SVector(1, 2, 3)], R * x[SVector(4, 5, 6)])) for x in s.states]
    path = tempname() * ".bsp"
    lock(EK._SPICE_LOCK) do
        h = SPICE.spkopn(path, "TEST", 0)
        SPICE.spkw13(h, -900003, 399, "ECLIPJ2000", ets[1], ets[end], "ECLIPTIC", 3, eclip, ets)
        SPICE.spkcls(h)
    end
    f = EK.read_spk(path)
    seg = only(f.segments)
    @test seg.frame == "ICRF" && maximum(maximum(abs.(a - b)) for (a, b) in zip(seg.states, s.states)) < 1e-9
    @test any(d -> occursin("ECLIPJ2000", d.message) && occursin("rotated to ICRF", d.message), f.diagnostics)
    @test maximum(abs.(spice_state(path, -900003, ets[7])[1:3] - seg.states[7][1:3])) < 1e-9

    # Writing: an ecliptic segment is refused by write_spk and rotated to ICRF by write_ephemeris,
    # so SPICE reads the position Epicycle's own frames give, not one 23 mas away.
    e = spk_segment(20); e.frame = "J2000_ECLIPTIC"
    @test_throws "ECLIPJ2000" EK.write_spk(tempname() * ".bsp", EK.EphemerisFileData(segments = [e]); naif_id = -900004)
    ec, icrf = CoordinateSystem(earth, MJ2000Ec()), CoordinateSystem(earth, ICRF())
    hist = SpacecraftHistory([HistorySegment(s.times, [CartesianState(Vector(x)) for x in s.states], ec)])
    p = tempname() * ".bsp"
    write_ephemeris(p, hist; naif_id = -900004)
    for k in (1, 10, 20)
        ours = to_vector(CartesianState(Coordinate(CartesianState(Vector(s.states[k])), ec, s.times[k]), icrf))[1:3]
        @test maximum(abs.(spice_state(p, -900004, ets[k])[1:3] - ours)) < 1e-9
    end

    # A usable stop before the last sample: kept on reading, bounding coverage and state_at, while
    # the samples beyond it remain for the interpolation.
    u = spk_segment(20); u.useable_stop = u.times[17]
    path = tempname() * ".bsp"
    EK.write_spk(path, EK.EphemerisFileData(segments = [u]); naif_id = -900005)
    eph = read_ephemeris(path; strict = true)
    m = only(eph.segments).metadata
    @test abs(EK._seconds(m["useable_stop"], u.times[17])) < 1e-6 && length(only(eph.segments).times) == 20
    @test abs(EK._seconds(only(coverage(eph))[2], u.times[17])) < 1e-6
    @test state_at(eph, u.times[16]) isa Coordinate
    @test_throws ArgumentError state_at(eph, u.times[19])
end

@testset "SPK: a damaged file or a bad naif_id is reported, not thrown" begin
    path = tempname() * ".bsp"
    EK.write_spk(path, EK.EphemerisFileData(segments = [spk_segment(20)]); naif_id = -900001)
    @test_throws "integer NAIF ID" EK.read_spk(path; naif_id = 1.5)
    # The type 13 window size is the second-to-last word of the segment; make it impossible.
    ic = only(lock(() -> EK._spk_summaries(EK._Diagnostics(path; strict = false), path), EK._SPICE_LOCK)).ic
    bytes = read(path)
    bytes[(ic[6] - 2) * 8 + 1:(ic[6] - 1) * 8] = reinterpret(UInt8, [-5.0])
    write(path, bytes)
    msg = try EK.read_spk(path); "" catch e sprint(showerror, e) end
    @test occursin("damaged", msg)
end

@testset "SPK: each error names its field and the problem" begin
    s = spk_segment(20; step = 60.0)
    states = [collect(v) for v in s.states]
    et0, et1 = EK._et(s.times[1]), EK._et(s.times[end])
    # A file written straight through SPICE: `segments` is a list of (center, frame, type) to write.
    function spk_file(segments...)
        p = tempname() * ".bsp"
        lock(EK._SPICE_LOCK) do
            h = SPICE.spkopn(p, "TEST", 0)
            for (center, frame, typ) in segments
                typ == 13 ? SPICE.spkw13(h, -900002, center, frame, et0, et1, "T", 7, states, EK._et.(s.times)) :
                            SPICE.spkw08(h, -900002, center, frame, et0, et1, "T", 7, states, et0, 60.0)
            end
            SPICE.spkcls(h)
        end
        return p
    end
    problem(path; kw...) = try EK.read_spk(path; kw...); nothing catch e only(e.problems) end
    text = tempname() * ".bsp"; write(text, "not a DAF file\n")
    cases = [
        ("not an SPK file", problem(text), "file", "could not read it as an SPK file"),
        ("a center Epicycle does not have", problem(spk_file((3, "J2000", 13))), "center", "NAIF body 3"),
        ("a frame Epicycle does not read", problem(spk_file((399, "GALACTIC", 13))), "frame", "SPICE frame 13"),
        ("a sampled type without step", problem(spk_file((399, "J2000", 8))), "step", "type 8, which stores no samples"),
        ("a step that is not positive", problem(spk_file((399, "J2000", 8)); step = 0.0), "step", "must be positive"),
        ("a step too long for the span", problem(spk_file((399, "J2000", 8)); step = 600.0), "step", "at least 8 are needed"),
    ]
    for (what, p, field, phrase) in cases
        found = p !== nothing && p.field == field && occursin(phrase, p.message)
        found || @info "not as expected: $what" p
        @test found
    end
end

@testset "SPK: the writer refuses, lists every problem, and writes nothing" begin
    path = tempname() * ".bsp"
    write(path, "existing")
    refused(data; kw...) = (msg = try EK.write_spk(path, data; kw...); "" catch e sprint(showerror, e) end;
                            @test read(path, String) == "existing"; msg)
    one(s) = EK.EphemerisFileData(segments = [s])
    @test occursin("Rotate the states to ICRF first: EME2000", refused(one(spk_segment(20; frame = "EME2000")); naif_id = -1))
    @test occursin("natural body", refused(one(spk_segment(20)); naif_id = 399))
    @test occursin("integer", refused(one(spk_segment(20)); naif_id = "SAT"))
    s = spk_segment(20); s.accelerations = [SVector(0.0, 0, 0) for _ in 1:20]
    @test occursin("discard_unsupported", refused(one(s); naif_id = -1))
    s = spk_segment(20); s.covariances = [EK.CovarianceRecord(s.times[1], "ICRF", SMatrix{6, 6}(zeros(6, 6)))]
    @test occursin("covariance", refused(one(s); naif_id = -1))
    @test occursin("over 27", refused(one(spk_segment(40; interp = Lagrange(31))); naif_id = -1))
    s = spk_segment(20); s.states[4] = SVector(7000.0, 0, 0, NaN, NaN, NaN)
    @test occursin("positions only", refused(one(s); naif_id = -1))
    a, b = spk_segment(20), spk_segment(20; t0 = "2026-10-21T00:10:00")
    @test occursin("overlap", refused(EK.EphemerisFileData(segments = [a, b]); naif_id = -1))
    # discard_unsupported writes without them.
    s = spk_segment(20); s.accelerations = [SVector(0.0, 0, 0) for _ in 1:20]
    p = tempname() * ".bsp"
    EK.write_spk(p, one(s); naif_id = -1, discard_unsupported = true)
    @test length(EK.read_spk(p).segments[1].times) == 20
end

@testset "SPK: OEM -> SPK -> OEM" begin
    seg = spk_segment(30)
    oem, spk, oem2 = tempname() * ".oem", tempname() * ".bsp", tempname() * ".oem"
    EK.write_oem(oem, EK.EphemerisFileData(originator = "EPICYCLE", segments = [seg]))
    EK.write_spk(spk, EK.read_oem(oem; strict = true); naif_id = -900001)
    back = EK.read_spk(spk)
    back.originator = "EPICYCLE"
    back.segments[1].time_system = "UTC"                                 # the file's system; epochs are absolute
    EK.write_oem(oem2, back)                                             # OEM 3.0, not the SPK file's version
    final = only(EK.read_oem(oem2; strict = true).segments)
    @test all(isapprox.(seg.states, final.states; rtol = 1e-15))
    @test maximum(abs(EK._seconds(a, b)) for (a, b) in zip(seg.times, final.times)) < 1e-6
    @test final.interpolation == Hermite(7)
end
