# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# JPL SPK reading and writing (ticket #9). CSPICE is the reference: what Epicycle writes, SPICE
# evaluates, and the result must equal Epicycle's own interpolation of the same samples.

using Test
using EpicycleIO
using AstroEpochs
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
function spice_state(path, target, et, centre = 399)
    lock(EK._SPICE_LOCK) do
        SPICE.furnsh(path)
        try
            return SPICE.spkgeo(target, et, "J2000", centre)[1]
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
end

@testset "SPK: the writer refuses, lists every problem, and writes nothing" begin
    path = tempname() * ".bsp"
    write(path, "existing")
    refused(data; kw...) = (msg = try EK.write_spk(path, data; kw...); "" catch e sprint(showerror, e) end;
                            @test read(path, String) == "existing"; msg)
    one(s) = EK.EphemerisFileData(segments = [s])
    @test occursin("Rotate EME2000 to ICRF", refused(one(spk_segment(20; frame = "EME2000")); naif_id = -1))
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
