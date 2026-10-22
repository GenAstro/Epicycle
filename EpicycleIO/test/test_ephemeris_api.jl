# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The public ephemeris interface (ticket #9): read_ephemeris, write_ephemeris, state_at,
# coverage, position_difference. Truth is a circular orbit in closed form.

using Test
using EpicycleIO
using AstroEpochs, AstroFrames, AstroUniverse, AstroModels, AstroStates
using AstroRoutines: Lagrange, Hermite
using StaticArrays

const EA = EpicycleIO

const W = 2π / 5400                                                    # rad/s
circ(s) = [7000cos(W * s), 7000sin(W * s), 0.0, -7000W * sin(W * s), 7000W * cos(W * s), 0.0]
const T0 = Time("2026-10-21T00:00:00", UTC(), ISOT())
at(s) = Time(T0.jd1, T0.jd2 + s / 86400, :utc, :jd)

function history(ranges...; cs = CoordinateSystem(earth, ICRF()), interp = nothing)
    segs = HistorySegment[]
    for r in ranges
        md = Dict{String, Any}()
        interp === nothing || (md["interpolation"] = interp)
        push!(segs, HistorySegment([at(s) for s in r], [CartesianState(circ(s)) for s in r], cs; metadata = md))
    end
    return SpacecraftHistory(segs)
end
pos(c) = to_vector(CartesianState(c))[1:3]

@testset "ephemeris API: read and write each format" begin
    h = history(0.0:60.0:3600.0)
    for (ext, kw) in ((".oem", (;)), (".e", (;)), (".bsp", (naif_id = -900001,)))
        p = tempname() * ext
        @test write_ephemeris(p, h; kw...) == p
        eph = read_ephemeris(p)
        @test eph isa SpacecraftHistory && length(epochs(eph)) == 61
        @test maximum(position_difference(eph, h)) < 1e-6                  # km; SPK epochs carry 0.1 µs
        # The format can also be named when the extension does not say.
        q = tempname() * ".dat"
        fmt = ext == ".oem" ? OEM() : ext == ".e" ? STKEphemeris() : SPK()
        write_ephemeris(q, h; format = fmt, kw...)
        @test length(epochs(read_ephemeris(q; format = fmt))) == 61
    end
    @test_throws ArgumentError write_ephemeris(tempname() * ".txt", h)
    @test occursin("NAIF", try write_ephemeris(tempname() * ".bsp", h); "" catch e sprint(showerror, e) end)
    @test_throws ArgumentError write_ephemeris(tempname() * ".oem", SpacecraftHistory())
end

@testset "ephemeris API: what a read history carries" begin
    h = history(0.0:60.0:1200.0; interp = Hermite(7))
    p = tempname() * ".oem"
    write_ephemeris(p, h; object_name = "SAT-1", object_id = "2026-001A", originator = "TESTS")
    eph = read_ephemeris(p)
    s = only(eph.segments)
    @test s.coordinate_system == CoordinateSystem(earth, ICRF())
    @test s.metadata["interpolation"] == Hermite(7)
    @test (s.metadata["object_name"], s.metadata["object_id"]) == ("SAT-1", "2026-001A")
    @test s.metadata["source"].originator == "TESTS" && s.metadata["source"].format == :oem
    @test isempty(diagnostics(eph))
    # A correction made on reading shows in diagnostics.
    write(p, replace(read(p, String), "ORIGINATOR" => "originator"))
    @test !isempty(@test_logs (:warn,) diagnostics(read_ephemeris(p)))
end

@testset "ephemeris API: state_at" begin
    h = history(0.0:60.0:1800.0; interp = Hermite(7))
    for s in (0.0, 61.7, 900.123, 1799.0, 1800.0)
        @test maximum(abs.(pos(state_at(h, at(s))) - circ(s)[1:3])) < 1e-9
    end
    @test state_at(h, at(100.0)).coord_sys == CoordinateSystem(earth, ICRF())
    # Any time scale.
    @test maximum(abs.(pos(state_at(h, at(100.0).tdb)) - circ(100.0)[1:3])) < 1e-8
    # No extrapolation, except within 1 µs of an end.
    @test occursin("no extrapolation", try state_at(h, at(1800.01)); "" catch e sprint(showerror, e) end)
    @test pos(state_at(h, at(1800.0 + 5e-7))) ≈ circ(1800.0)[1:3]
    # A gap: the error names the covered times on both sides.
    g = history(0.0:60.0:600.0, 1200.0:60.0:1800.0)
    msg = try state_at(g, at(900.0)); "" catch e sprint(showerror, e) end
    @test occursin("gap", msg) && occursin("00:10:00", msg) && occursin("00:20:00", msg)
    @test coverage(g) == [(at(0.0), at(600.0)), (at(1200.0), at(1800.0))] && span(g) == (at(0.0), at(1800.0))
    # A shared endpoint, as at a maneuver: the later segment, the state after the event.
    m = history(0.0:60.0:600.0, 600.0:60.0:1200.0)
    m.segments[2].states[1] = CartesianState(circ(600.0) .+ [0, 0, 0, 0, 0, 0.01])
    @test to_vector(CartesianState(state_at(m, at(600.0))))[6] ≈ 0.01
    @test length(coverage(m)) == 1
    # Overlap: an error unless the caller chooses.
    o = history(0.0:60.0:900.0, 600.0:60.0:1500.0)
    o.segments[2].states .= [CartesianState(circ(s) .+ [1.0, 0, 0, 0, 0, 0]) for s in 600.0:60.0:1500.0]
    @test occursin("overlap = :earlier or :later", try state_at(o, at(700.0)); "" catch e sprint(showerror, e) end)
    @test pos(state_at(o, at(700.0); overlap = :earlier))[1] ≈ circ(700.0)[1] atol = 1e-8
    @test pos(state_at(o, at(700.0); overlap = :later))[1] ≈ circ(700.0)[1] + 1 atol = 1e-8
end

@testset "ephemeris API: frames, steps and refusals" begin
    h = history(0.0:30.0:3600.0)
    # Written in another coordinate system, read back, and compared in the original: the same orbit.
    p = tempname() * ".oem"
    write_ephemeris(p, h; coordinate_system = CoordinateSystem(earth, ITRF()))
    eph = read_ephemeris(p)
    @test only(eph.segments).coordinate_system == CoordinateSystem(earth, ITRF())
    @test maximum(position_difference(h, eph)) < 1e-7
    # SPK is written in ICRF axes whatever the source's, without being asked.
    e2000 = history(0.0:30.0:3600.0; cs = CoordinateSystem(earth, MJ2000Eq()))
    q = tempname() * ".bsp"
    write_ephemeris(q, e2000; naif_id = -42)
    @test only(read_ephemeris(q).segments).coordinate_system == CoordinateSystem(earth, ICRF())
    @test maximum(position_difference(e2000, read_ephemeris(q))) < 1e-6
    # step resamples, keeping each segment's ends, and needs no source epoch on the grid.
    r = tempname() * ".oem"
    write_ephemeris(r, history(0.0:7.3:1000.0); step = 120.0)
    t = epochs(read_ephemeris(r))
    last_s = last(0.0:7.3:1000.0)                                     # 992.8 s
    @test length(t) == 10 && abs(EA._seconds(t[end], at(last_s))) < 1e-9   # 0, 120, …, 960, then the end
    @test maximum(position_difference(read_ephemeris(r), history(0.0:7.3:1000.0))) < 1e-6
    # Data a conversion or resampling would leave half-done are refused.
    c = history(0.0:60.0:600.0)
    c.segments[1].metadata["covariances"] = [EA.CovarianceRecord(at(0.0), "ICRF", SMatrix{6, 6}(zeros(6, 6)))]
    @test occursin("rotates states only", try write_ephemeris(tempname() * ".oem", c;
        coordinate_system = CoordinateSystem(earth, MJ2000Eq())); "" catch e sprint(showerror, e) end)
    @test occursin("drop the accelerations and covariance", try write_ephemeris(tempname() * ".oem", c; step = 60.0); "" catch e sprint(showerror, e) end)
    # The same data without conversion keep their covariance through a round trip.
    d = tempname() * ".oem"
    write_ephemeris(d, c)
    @test length(only(read_ephemeris(d).segments).metadata["covariances"]) == 1
end

@testset "ephemeris API: one trajectory through every format" begin
    h = history(0.0:60.0:1800.0, 1800.0:60.0:3600.0; interp = Hermite(7))
    h.segments[2].states[1] = CartesianState(circ(1800.0) .+ [0, 0, 0, 0.001, 0, 0])   # a maneuver
    a, b, c, d = tempname() .* (".oem", ".bsp", ".e", ".oem")
    write_ephemeris(a, h)
    write_ephemeris(b, read_ephemeris(a); naif_id = -900001)
    write_ephemeris(c, read_ephemeris(b))
    write_ephemeris(d, read_ephemeris(c))
    final = read_ephemeris(d)
    @test length(final.segments) == 2
    @test maximum(position_difference(h, final)) < 1e-6
    @test to_vector(final.segments[2].states[1])[4] ≈ to_vector(h.segments[2].states[1])[4]   # the maneuver survives
    @test final.segments[1].metadata["interpolation"] == Hermite(7)
end

@testset "ephemeris API: one-state segments" begin
    # A maneuver records one state, the first of the propagation after it, in the spacecraft's
    # frame: written files leave it out, and state_at evaluates it at its epoch.
    h = history(0.0:60.0:600.0, 600.0:60.0:1200.0)
    burn = HistorySegment([at(600.0)],
        [CartesianState(to_vector(CartesianState(Coordinate(h.segments[2].states[1], h.segments[2].coordinate_system, at(600.0)),
                                                CoordinateSystem(earth, MJ2000Eq()))))],
        CoordinateSystem(earth, MJ2000Eq()); name = "maneuver")
    insert!(h.segments, 2, burn)
    p = tempname() * ".oem"
    write_ephemeris(p, h)
    @test length(read_ephemeris(p).segments) == 2
    @test pos(state_at(SpacecraftHistory([burn]), at(600.0))) ≈ to_vector(burn.states[1])[1:3]
    # A one-state segment that repeats nothing cannot be written.
    lone = SpacecraftHistory([HistorySegment([at(5000.0)], [CartesianState(circ(5000.0))], CoordinateSystem(earth, ICRF()))])
    @test occursin("one state", try write_ephemeris(tempname() * ".oem", lone); "" catch e sprint(showerror, e) end)
end
