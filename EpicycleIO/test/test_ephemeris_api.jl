# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The public ephemeris interface (ticket #9): read_ephemeris, write_ephemeris, state_at,
# coverage, position_difference. Truth is a circular orbit in closed form.

using Test
using EpicycleIO
using AstroEpochs, AstroFrames, AstroUniverse, AstroModels, AstroStates
using AstroRoutines: Lagrange, Hermite
import AstroRoutines
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

@testset "ephemeris API: an OEM in each time system" begin
    # The first epoch, 2026-10-21T00:00:00 UTC, as each system reads it. UTC, TAI, TT and GPS
    # differ by whole offsets (37 leap seconds in 2026); TCG runs ahead of TT at the rate L_G
    # from 1977-01-01T00:00:32.184 TT; TDB is within 2 ms of TT; TCB runs ahead of TDB at L_B
    # from the same origin, less TDB₀.
    function seconds_of_day(text)
        hour, minute, second = parse.(Float64, split(split(text, 'T')[2], ':'))
        return 3600hour + 60minute + second
    end
    jd_tt = 2461334.5 + (37 + 32.184) / 86400
    tcg = 32.184 + 37 + 6.969290134e-10 * (jd_tt - 2443144.5003725) * 86400
    tcb = 32.184 + 37 + 1.550519768e-8 * (jd_tt - 2443144.5003725) * 86400 - 6.55e-5
    expected = Dict("UTC" => (0.0, 1e-9), "TAI" => (37.0, 1e-9), "TT" => (69.184, 1e-9), "GPS" => (18.0, 1e-9),
                    "TCG" => (tcg, 1e-6), "TDB" => (69.184, 2e-3), "TCB" => (tcb, 2e-3))
    h = history(0.0:60.0:600.0)
    for (name, scale) in (("UTC", UTC()), ("TAI", TAI()), ("TT", TT()), ("GPS", "GPS"), ("TCG", TCG()), ("TDB", TDB()), ("TCB", TCB()))
        p = tempname() * ".oem"
        write_ephemeris(p, h; time_system = scale)
        lines = readlines(p)
        @test any(l -> startswith(l, "TIME_SYSTEM") && endswith(l, " $name"), lines)
        first_epoch = first(filter(l -> occursin(r"^2026-", l), lines))
        offset, tolerance = expected[name]
        @test abs(seconds_of_day(split(first_epoch)[1]) - offset) < tolerance
        # Read back, every epoch is the instant it was written for.
        back = read_ephemeris(p; strict = true)
        @test only(back.segments).metadata["source"].format == :oem
        @test maximum(abs(EA._seconds(a, b)) for (a, b) in zip(epochs(back), epochs(h))) < 1e-9
    end
    # ET is TDB by another name; it reads, with a note.
    p = tempname() * ".oem"
    write_ephemeris(p, h; time_system = TDB())
    write(p, replace(read(p, String), "TIME_SYSTEM          = TDB" => "TIME_SYSTEM          = ET"))
    et = read_ephemeris(p)
    @test maximum(abs(EA._seconds(a, b)) for (a, b) in zip(epochs(et), epochs(h))) < 1e-9
end

@testset "ephemeris API: every writing option reads back strictly" begin
    # A file written with any combination of options follows its format exactly, so a strict read
    # accepts it, and it holds the trajectory it was written from: within 1 mm, the resampling
    # and frame rotation included.
    h = history(0.0:60.0:1800.0)
    formats = ((".oem", (;)), (".e", (;)), (".bsp", (naif_id = -900001,)))
    steps = (nothing, 60.0, 37.5)                                          # 37.5 s leaves a short last step
    systems = (nothing, CoordinateSystem(earth, MJ2000Eq()))
    methods = (nothing, Lagrange(5), Hermite(5))
    for (ext, kw) in formats, step in steps, cs in systems, m in methods
        ext == ".bsp" && cs !== nothing && continue                           # SPK is ICRF; refused below
        p = tempname() * ext
        opts = (; kw..., step, coordinate_system = cs, interpolation = m)
        write_ephemeris(p, h; opts...)
        back = try read_ephemeris(p; strict = true) catch e
            @info "strict read refused" ext step cs m e
            nothing
        end
        @test back !== nothing
        back === nothing && continue
        @test maximum(position_difference(back, h)) < 1e-6                       # km
        m === nothing || @test only(back.segments).metadata["interpolation"] == m
        step === nothing || @test length(epochs(back)) == length(0.0:step:1800.0) + (1800.0 % step == 0 ? 0 : 1)
    end
    # SPK holds ICRF or ecliptic axes only, and says how to write the trajectory instead.
    @test_throws "leave coordinate_system out" write_ephemeris(tempname() * ".bsp", h; naif_id = -900001,
                                                               coordinate_system = CoordinateSystem(earth, MJ2000Eq()))
    # OEM alone: the time system and the decimal places of seconds, each read back strictly.
    for sys in (UTC(), TDB(), "GPS"), digits in (0, 3, 9)
        p = tempname() * ".oem"
        # UTC and GPS epochs of this grid fall on whole seconds; TDB's do not, and rounding them
        # would move the trajectory, so the writer refuses.
        if sys isa TDB && digits < 9
            @test_throws "write more decimal places" write_ephemeris(p, h; time_system = sys, digits)
            continue
        end
        write_ephemeris(p, h; time_system = sys, digits)
        back = read_ephemeris(p; strict = true)
        @test maximum(position_difference(back, h)) < 1e-6
        epoch_text = split(first(filter(l -> occursin(r"^2026-", l), readlines(p))))[1]
        @test digits == 0 ? !occursin('.', epoch_text) : length(split(epoch_text, '.')[2]) == digits
    end
end

@testset "ephemeris API: what each format keeps through a round trip" begin
    # Two segments about the Earth and then about the Moon, each with a name, an ID, Hermite 5 and
    # a usable span two minutes inside its data.
    for body in (earth, moon)
        cs = CoordinateSystem(body, ICRF())
        segs = map((0.0:60.0:1800.0, 1800.0:60.0:3600.0)) do r
            md = Dict{String, Any}("object_name" => "PROBE", "object_id" => "2026-099A", "interpolation" => Hermite(5),
                                   "useable_start" => at(first(r) + 120), "useable_stop" => at(last(r) - 120))
            HistorySegment([at(s) for s in r], [CartesianState(circ(s)) for s in r], cs; metadata = md)
        end
        h = SpacecraftHistory(collect(segs))
        near(a, b) = abs(EA._seconds(a, b)) < 1e-6
        # OEM and SPK keep every one of them; SPK names the object by its NAIF ID.
        for (ext, kw, id) in ((".oem", (; time_system = TDB()), "2026-099A"), (".bsp", (; naif_id = -900099), "-900099"))
            p = tempname() * ext
            write_ephemeris(p, h; kw...)
            back = read_ephemeris(p; strict = true)
            @test length(back.segments) == 2
            for (a, b) in zip(h.segments, back.segments)
                m = b.metadata
                @test b.coordinate_system == cs
                @test m["object_name"] == "PROBE" && m["object_id"] == id && m["interpolation"] == Hermite(5)
                @test near(m["useable_start"], a.metadata["useable_start"]) && near(m["useable_stop"], a.metadata["useable_stop"])
                @test near(first(b.times), first(a.times)) && near(last(b.times), last(a.times))
            end
            @test maximum(position_difference(back, h)) < 1e-6
        end
        # STK has neither an object's name nor a usable span: the writer refuses rather than widen
        # the coverage, and with discard_unsupported writes the data, read back as UNKNOWN.
        p = tempname() * ".e"
        @test_throws "an STK file cannot say so" write_ephemeris(p, h)
        write_ephemeris(p, h; discard_unsupported = true)
        back = read_ephemeris(p; strict = true)
        @test length(back.segments) == 2 && all(s -> s.coordinate_system == cs, back.segments)
        @test all(s -> s.metadata["object_name"] == "UNKNOWN" && !haskey(s.metadata, "useable_start"), back.segments)
        @test maximum(position_difference(back, h)) < 1e-6
    end
end

@testset "ephemeris API: resampling keeps the usable span" begin
    h = history(0.0:60.0:3600.0)
    h.segments[1].metadata["useable_start"] = at(600.0)
    h.segments[1].metadata["useable_stop"] = at(3000.0)
    p = tempname() * ".oem"
    write_ephemeris(p, h; step = 120.0)
    s = only(read_ephemeris(p).segments).metadata
    @test abs(EA._seconds(s["useable_start"], at(600.0))) < 1e-9
    @test abs(EA._seconds(s["useable_stop"], at(3000.0))) < 1e-9
end

@testset "ephemeris API: a classification travels with the trajectory" begin
    seg = EA.EphemerisSegmentData(object_name = "SAT", object_id = "2026-001A", center = "EARTH", frame = "EME2000",
                                  time_system = "UTC", times = [at(s) for s in 0.0:60.0:600.0],
                                  states = [SVector{6}(circ(s)) for s in 0.0:60.0:600.0])
    marked = tempname() * ".oem"
    EA.write_oem(marked, EA.EphemerisFileData(originator = "TEST", message_id = "ID-123", classification = "SECRET",
                                              segments = [seg]))
    eph = read_ephemeris(marked; strict = true)
    out = tempname() * ".oem"
    write_ephemeris(out, eph)
    lines = readlines(out)
    @test any(l -> startswith(l, "CLASSIFICATION") && endswith(l, "SECRET"), lines)
    @test !any(startswith("MESSAGE_ID"), lines)                         # a new message, not the source's ID
    write_ephemeris(out, eph; message_id = "ID-124")
    @test any(l -> startswith(l, "MESSAGE_ID") && endswith(l, "ID-124"), readlines(out))
    # STK and SPK cannot carry the marking: refused unless the caller drops it on purpose.
    msg = try write_ephemeris(tempname() * ".e", eph); "" catch e sprint(showerror, e) end
    @test occursin("SECRET", msg) && occursin("classification = \"\"", msg)
    @test isfile(write_ephemeris(tempname() * ".e", eph; classification = ""))
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

@testset "ephemeris API: state_at uses only its window, with the same result" begin
    # The window state_at picks by binary search is the one AstroRoutines picks in the whole
    # table, so the result is identical; uneven spacing exercises the odd-count rule.
    ts = cumsum(vcat(0.0, [40.0 + 25.0 * sin(k) for k in 1:299]))
    full_x = ts
    for m in (Lagrange(7), Lagrange(6), Hermite(7), Hermite(5))
        h = history(ts; interp = m)
        r = [SVector{3}(circ(s)[1:3]) for s in ts]
        v = [SVector{3}(circ(s)[4:6]) for s in ts]
        for s in (0.0, 13.7, ts[2] + 1.0, ts[150] + 3.3, ts[end] - 0.5, ts[end])
            got = pos(state_at(h, at(s)))
            want = m isa Hermite ? EA.interpolate(m, full_x, r, v, s) : EA.interpolate(m, full_x, r, s)
            @test maximum(abs.(got - want)) < 1e-8                      # km; epochs pass through Time
        end
    end
    # Positions only: Lagrange through positions, velocity from its derivative.
    h = history(ts)
    h.segments[1].states .= [CartesianState(vcat(circ(s)[1:3], fill(NaN, 3))) for s in ts]
    c = to_vector(CartesianState(state_at(h, at(ts[100] + 7.0))))
    @test maximum(abs.(c[4:6] - circ(ts[100] + 7.0)[4:6])) < 1e-8
end

@testset "ephemeris API: cost does not grow with the table" begin
    small, large = history(0.0:10.0:20_000.0), history(0.0:10.0:200_000.0)       # 2 001 and 20 001 points
    probe(h) = (t = at(1234.5); state_at(h, t); @elapsed for _ in 1:200; state_at(h, t); end)
    probe(small)
    @test probe(large) < 20 * max(probe(small), 1e-4)                             # O(log N), not O(N)
    @test (@elapsed position_difference(large, large)) < 30.0
end

# The GMAT reference files committed with the examples: every format of one trajectory.
const GMAT_EXAMPLES = joinpath(@__DIR__, "..", "..", "Epicycle", "examples")

# GMAT propagated point-mass gravity with Epicycle's μ, so the truth is the two-body orbit through
# each file's first state. Measured 2026-10-24: the Earth orbit's OEM within 1.2 µm of it at its
# epochs and 46 µm between them (Lagrange 7 at 60 s), its .e and SPK within 6 and 8 mm, the error
# of GMAT's rounded epochs there; the low lunar orbit's three files within 3.4 and 4.1 mm. Each
# tolerance is about twice the measured value.
function two_body(x0, dt, μ)
    a, e, i, Ω, ω, ν = cart_to_kep(x0, μ)
    M = AstroRoutines.true_to_mean_anomaly(ν, e) + sqrt(μ / a^3) * dt
    return kep_to_cart([a, e, i, Ω, ω, AstroRoutines.mean_to_true_anomaly(mod2pi(M), e)], μ)
end
@testset "ephemeris API: GMAT's files are the two-body orbit, at their epochs and between" begin
    gmat_data = joinpath(@__DIR__, "data", "gmat")
    for (path, μ, body, at_epochs, between) in (
            (joinpath(GMAT_EXAMPLES, "Ephemeris_OEM", "gmat_leo.oem"), 398600.4418, earth, 3e-9, 1e-7),     # km
            (joinpath(GMAT_EXAMPLES, "Ephemeris_STK", "gmat_leo.e"), 398600.4418, earth, 1.2e-5, 1.6e-5),
            (joinpath(GMAT_EXAMPLES, "Ephemeris_SPK", "gmat_leo.bsp"), 398600.4418, earth, 1.1e-5, 1.0e-5),
            (joinpath(gmat_data, "gmat_luna.oem"), 4902.8, moon, 7e-6, 8e-6),
            (joinpath(gmat_data, "gmat_luna.e"), 4902.8, moon, 7e-6, 8e-6),
            (joinpath(gmat_data, "gmat_luna.bsp"), 4902.8, moon, 3e-6, 3e-6))
        eph = read_ephemeris(path; strict = true)
        s = only(eph.segments)
        @test s.coordinate_system.origin == body && length(s.times) == 121
        x0, t0 = to_vector(s.states[1]), s.times[1]
        truth(t) = two_body(x0, EA._seconds(t0, t), μ)
        dist(x, y) = sqrt(sum(abs2, x[1:3] - y[1:3]))
        # The stored states, and state_at returning them exactly.
        @test maximum(dist(to_vector(x), truth(t)) for (t, x) in zip(s.times, s.states)) < at_epochs
        @test all(to_vector(CartesianState(state_at(eph, t))) == to_vector(x) for (t, x) in zip(s.times, s.states))
        # Halfway between each pair, interpolated with the file's own method.
        mids = [Time(a.jd1, a.jd2 + EA._seconds(a, b) / 2 / 86400, a.scale, :jd) for (a, b) in zip(s.times[1:end-1], s.times[2:end])]
        @test maximum(dist(to_vector(CartesianState(state_at(eph, t))), truth(t)) for t in mids) < between
    end
end

@testset "ephemeris API: GMAT's OEM, .e and SPK of one orbit agree" begin
    oem = read_ephemeris(joinpath(GMAT_EXAMPLES, "Ephemeris_OEM", "gmat_leo.oem"); strict = true)
    stk = read_ephemeris(joinpath(GMAT_EXAMPLES, "Ephemeris_STK", "gmat_leo.e"); strict = true)
    spk = read_ephemeris(joinpath(GMAT_EXAMPLES, "Ephemeris_SPK", "gmat_leo.bsp"))
    for e in (oem, stk, spk)
        @test length(epochs(e)) == 121 && length(e.segments) == 1
    end
    @test oem.segments[1].metadata["interpolation"] == Lagrange(7)
    @test stk.segments[1].metadata["interpolation"] == Lagrange(5)       # GMAT writes no InterpolationSamplesM1
    @test spk.segments[1].metadata["interpolation"] == Hermite(7)
    @test to_vector(oem.segments[1].states[1]) == [7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]
    # The OEM and .e are both EME2000 and agree to millimetres, GMAT's .e epochs carrying
    # about 1 µs of rounding.
    @test maximum(position_difference(oem, stk)) < 1e-5                # km
    # GMAT's SPK and OEM of the same orbit differ, read by Epicycle under each format's own
    # frame definition (CCSDS EME2000 as FK5 J2000, NAIF J2000 as ICRF), by a fixed rotation of
    # a few tens of mas and GMAT's SPK epochs 25.6 µs from its OEM's: about 1.3 m here. This
    # records the measured difference so a change in either is seen; it is an open question
    # how GMAT's EME2000 relates to ICRF, not a tolerance.
    d = position_difference(oem, spk)
    @test 0.5e-3 < maximum(d) < 2.0e-3
    gap = [EA._seconds(a, b) for (a, b) in zip(epochs(oem), epochs(spk))]
    @test all(g -> -30e-6 < g < -20e-6, gap)
end
