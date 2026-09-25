# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# The 3D path: `orbitview`, the CZML it writes, and the segment extraction
# underneath it.
#
# None of this can be checked by looking at the browser, so it is checked at the
# packet. The CZML took a long time to get right in the prototype — the inertial
# reference frame, the per-segment interval tags that stop Lagrange interpolation
# smoothing across a maneuver, the km-to-m scaling — and none of it is obvious
# enough to survive a refactor unwatched. That is what this file guards.
#
# Most of it needs no spacecraft: `_build_czml` takes segments directly, so the
# writer is tested on data made up here. Only `_extract` needs a real history.
# =============================================================================

using Test
using EpicycleIO
using JSON

using AstroEpochs
using AstroModels
using AstroStates
using AstroFrames
using AstroUniverse

const IO_ = EpicycleIO

IO_.auto_open!(false)

const _OV_EPOCH = Time(2458849.5, 0.0, :tdb, :jd)

"""A spacecraft with two recorded arcs and a 1-sample marker between them.

Built by hand rather than propagated: the point is the shape of the history, and
propagating to get it would make this test depend on the integrator.
"""
function _recorded(; gap_sample = true, marker_name = "TOI")
    cs  = CoordinateSystem(earth, ICRF())
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), _OV_EPOCH;
                     coord_sys = cs, name = "Explorer")

    # The second arc flies faster than the first, so the marker between them
    # stands for a real 0.2 km/s burn rather than a boundary that merely shares
    # a time.
    function arc(name, k0, n; vy = 7.35)
        seg = AstroModels.HistorySegment(cs; name = name)
        for k in k0:(k0 + n - 1)
            push!(seg.times,  Time(2458849.5 + k * 0.01, 0.0, :tdb, :jd))
            push!(seg.states, CartesianState([7000.0 + 10k, 100.0k, 1300.0, 0.0, vy, 1.0]))
        end
        return seg
    end

    push!(sat.history.segments, arc("coast 1", 0, 4))
    if gap_sample
        # A maneuver marker: one sample, which must not become an arc.
        push!(sat.history.segments, arc(marker_name, 4, 1; vy = 7.55))
    end
    push!(sat.history.segments, arc("coast 2", 4, 4; vy = 7.55))
    return sat
end

"""A history with three burns that all carry the same name.

This is the shape `solve_trajectory!` actually produces — it names every marker
it records "maneuver" — and it is what made three burns draw as one point.
"""
function _recorded_three_burns()
    cs  = CoordinateSystem(earth, ICRF())
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), _OV_EPOCH;
                     coord_sys = cs, name = "GeoSat-1")

    function seg!(name, k0, n, vy)
        s = AstroModels.HistorySegment(cs; name = name)
        for k in k0:(k0 + n - 1)
            push!(s.times,  Time(2458849.5 + k * 0.01, 0.0, :tdb, :jd))
            push!(s.states, CartesianState([7000.0 + 10k, 100.0k, 1300.0, 0.0, vy, 1.0]))
        end
        push!(sat.history.segments, s)
    end

    # Three burns of different sizes — 0.2, 0.3 and 0.5 km/s — so a ΔV read from
    # the wrong pair of neighbours shows up as the wrong number rather than
    # passing by coincidence.
    seg!("propagate", 0, 4, 7.35)
    seg!("maneuver",  4, 1, 7.55)
    seg!("propagate", 4, 4, 7.55)
    seg!("maneuver",  8, 1, 7.85)
    seg!("propagate", 8, 4, 7.85)
    seg!("maneuver", 12, 1, 8.35)
    seg!("propagate", 12, 4, 8.35)
    return sat
end

_seg(times, name = "arc", color = (255, 0, 0, 255), width = 1.5) =
    IO_._Segment(times,
                 [(7000.0 + t, 2.0 * t, 3.0 * t) for t in times],
                 name, color, width)

@testset "orbitview" begin

# ── Pure helpers ─────────────────────────────────────────────────────────────
@testset "an ISO string shifts by seconds and keeps exactly one Z" begin
    @test IO_._iso_shift("2020-01-01T00:00:00", 0) == "2020-01-01T00:00:00Z"
    @test IO_._iso_shift("2020-01-01T00:00:00Z", 90) == "2020-01-01T00:01:30Z"
    @test IO_._iso_shift("2020-01-01T00:00:00Z", -30) == "2019-12-31T23:59:30Z"
    @test IO_._iso_shift("2020-01-01T00:00:00Z", 0.5) == "2020-01-01T00:00:00.500Z"

    # Given a Z, do not end up with two.
    @test count(==('Z'), IO_._iso_shift("2020-01-01T00:00:00Z", 1)) == 1
end

@testset "interpolation stays inside the data" begin
    pairs = Tuple{Float64, NTuple{3, Float64}}[
        (0.0, (0.0, 0.0, 0.0)), (10.0, (10.0, 20.0, 30.0))]

    @test IO_._interpolate(pairs, 5.0) == (5.0, 10.0, 15.0)
    @test IO_._interpolate(pairs, 0.0) == (0.0, 0.0, 0.0)
    @test IO_._interpolate(pairs, 10.0) == (10.0, 20.0, 30.0)

    # Outside the span returns nothing rather than extrapolating a position that
    # was never computed — the caller skips the marker instead of inventing one.
    @test IO_._interpolate(pairs, -0.1) === nothing
    @test IO_._interpolate(pairs, 10.1) === nothing
    @test IO_._interpolate(Tuple{Float64, NTuple{3, Float64}}[], 0.0) === nothing
end

@testset "the palette cycles and matches the plots" begin
    @test IO_._pick_color(nothing, 1) == IO_.DEFAULT_PALETTE[1]
    @test IO_._pick_color(nothing, length(IO_.DEFAULT_PALETTE) + 1) == IO_.DEFAULT_PALETTE[1]
    @test IO_._pick_color([(1, 2, 3, 4), (5, 6, 7, 8)], 2) == (5, 6, 7, 8)
    @test IO_._pick_color([(1, 2, 3, 4), (5, 6, 7, 8)], 3) == (1, 2, 3, 4)   # cycles
    @test IO_._pick_color((9, 8, 7, 6), 4) == (9, 8, 7, 6)                   # a lone colour
end

@testset "an entity id survives a URL and a JSON key" begin
    @test IO_._entity_id((; name = "Explorer")) == "Explorer"
    @test IO_._entity_id((; name = "Sat A/B #1")) == "Sat_A_B__1"
    @test IO_._entity_id((; name = "π-craft")) == "_-craft"
end

# ── The CZML writer ──────────────────────────────────────────────────────────
@testset "the document packet carries the clock" begin
    packets = IO_._build_czml([_seg([0.0, 10.0, 20.0])], IO_._Maneuver[],
                              "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 60.0, inertial = true)

    doc = packets[1]
    @test doc["id"] == "document"
    @test doc["version"] == "1.0"

    clock = doc["clock"]
    @test clock["interval"] == "2020-01-01T00:00:00Z/2020-01-01T00:00:20Z"
    @test clock["currentTime"] == "2020-01-01T00:00:00Z"
    @test clock["multiplier"] == 60.0
    @test clock["range"] == "LOOP_STOP"
end

@testset "one packet per arc, in metres" begin
    seg = _seg([0.0, 10.0, 20.0], "coast")
    packets = IO_._build_czml([seg], IO_._Maneuver[], "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 1.0, inertial = true)

    arc = packets[2]
    @test arc["id"] == "sat-seg-1"
    @test arc["name"] == "coast"

    # cartesian is flat: time, x, y, z, time, x, y, z, … and CZML wants metres
    # while the history is in kilometres. Getting this wrong puts the orbit
    # inside the planet, which is why it is asserted rather than eyeballed.
    cart = arc["position"]["cartesian"]
    @test length(cart) == 4 * 3
    @test cart[1] == 0.0
    @test cart[2] == 7000.0 * 1000
    @test cart[3] == 0.0 * 1000
    @test cart[5] == 10.0                    # second sample's time
    @test cart[6] == 7010.0 * 1000

    @test arc["position"]["interpolationAlgorithm"] == "LAGRANGE"
    @test arc["position"]["forwardExtrapolationType"] == "NONE"
    @test arc["path"]["width"] == 1.5
    @test arc["path"]["material"]["solidColor"]["color"]["rgba"] == [255, 0, 0, 255]
end

@testset "the camera sits in the inertial frame, or is told not to" begin
    kw = (; id = "sat", name = "Sat", marker_color = (255, 255, 255, 255),
            marker_pixel_size = 10, multiplier = 1.0)

    inertial = IO_._build_czml([_seg([0.0, 10.0])], IO_._Maneuver[],
                               "2020-01-01T00:00:00Z"; inertial = true, kw...)
    @test inertial[2]["position"]["referenceFrame"] == "INERTIAL"
    @test inertial[end]["position"][1]["referenceFrame"] == "INERTIAL"

    # Without it Cesium defaults to the fixed frame, and the orbit appears to
    # precess because the camera turns with the ground.
    fixed = IO_._build_czml([_seg([0.0, 10.0])], IO_._Maneuver[],
                            "2020-01-01T00:00:00Z"; inertial = false, kw...)
    @test !haskey(fixed[2]["position"], "referenceFrame")
end

@testset "interpolation degree never exceeds the samples available" begin
    kw = (; id = "sat", name = "Sat", marker_color = (255, 255, 255, 255),
            marker_pixel_size = 10, multiplier = 1.0, inertial = true)

    two = IO_._build_czml([_seg([0.0, 10.0])], IO_._Maneuver[],
                          "2020-01-01T00:00:00Z"; kw...)
    @test two[2]["position"]["interpolationDegree"] == 1      # min(5, 2-1)

    many = IO_._build_czml([_seg(collect(0.0:10.0:100.0))], IO_._Maneuver[],
                           "2020-01-01T00:00:00Z"; kw...)
    @test many[2]["position"]["interpolationDegree"] == 5     # capped
end

@testset "the marker is interval-tagged per arc" begin
    # This is the subtle one. The marker's position is an ARRAY of interval-tagged
    # properties, one per arc, so Lagrange interpolation never runs across the
    # boundary between two arcs. A single property would smooth the marker
    # straight through a velocity discontinuity.
    segs = [_seg([0.0, 10.0, 20.0], "coast 1"), _seg([20.0, 30.0, 40.0], "coast 2")]
    packets = IO_._build_czml(segs, IO_._Maneuver[], "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 1.0, inertial = true)

    @test length(packets) == 1 + 2 + 1            # document, two arcs, one marker

    marker = packets[end]
    @test marker["id"] == "sat"
    @test marker["label"]["text"] == "Sat"
    @test marker["position"] isa AbstractVector
    @test length(marker["position"]) == 2
    @test marker["position"][1]["interval"] ==
          "2020-01-01T00:00:00Z/2020-01-01T00:00:20Z"
    @test marker["position"][2]["interval"] ==
          "2020-01-01T00:00:20Z/2020-01-01T00:00:40Z"
end

@testset "a maneuver point lands on the trajectory" begin
    segs = [_seg([0.0, 10.0, 20.0])]
    mnv  = [IO_._Maneuver("TOI", 10.0, (0.0, 0.1, 0.0))]
    packets = IO_._build_czml(segs, mnv, "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 1.0, inertial = true)

    pt = packets[end]
    @test pt["id"] == "sat-mnv-1-TOI"          # indexed, so repeated names cannot collide
    @test pt["position"]["cartesian"] == [7010.0 * 1000, 20.0 * 1000, 30.0 * 1000]
    @test occursin("TOI", pt["label"]["text"])
    @test occursin("0.1", pt["label"]["text"])              # the ΔV magnitude

    # A maneuver outside the recorded span is skipped rather than placed wrongly.
    outside = IO_._build_czml(segs, [IO_._Maneuver("Late", 999.0, (0.0, 0.1, 0.0))],
                              "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 1.0, inertial = true)
    @test length(outside) == 1 + 1 + 1                      # no maneuver packet
end

@testset "the writer refuses what it cannot draw" begin
    kw = (; id = "sat", name = "Sat", marker_color = (255, 255, 255, 255),
            marker_pixel_size = 10, multiplier = 1.0, inertial = true)

    @test_throws ArgumentError IO_._build_czml(IO_._Segment[], IO_._Maneuver[],
                                                "2020-01-01T00:00:00Z"; kw...)
    # One sample is a point, not an arc — Lagrange has nothing to work with.
    @test_throws ArgumentError IO_._build_czml([_seg([0.0])], IO_._Maneuver[],
                                                "2020-01-01T00:00:00Z"; kw...)
end

@testset "the packets are JSON, which is how they reach the page" begin
    packets = IO_._build_czml([_seg([0.0, 10.0])], IO_._Maneuver[],
                              "2020-01-01T00:00:00Z";
                              id = "sat", name = "Sat",
                              marker_color = (255, 255, 255, 255), marker_pixel_size = 10,
                              multiplier = 1.0, inertial = true)
    round_trip = JSON.parse(JSON.json(packets))
    @test round_trip[1]["id"] == "document"
    @test round_trip[2]["position"]["cartesian"][2] == 7000.0 * 1000
end

# ── Extraction from a real history ───────────────────────────────────────────
@testset "arcs come off the spacecraft, markers do not become arcs" begin
    sat = _recorded()
    epoch_iso, segments, maneuvers = IO_._extract(sat)

    # Three history segments, but the 1-sample maneuver marker is not an arc.
    @test length(sat.history.segments) == 3
    @test length(segments) == 2
    @test [s.name for s in segments] == ["coast 1", "coast 2"]

    # Times are seconds from the earliest sample, not absolute.
    #
    # The tolerance is not decoration. These come from subtracting two Julian
    # dates near 2.46e6, where a Float64 ulp is about 40 microseconds, so the
    # difference carries that much cancellation error however exact the inputs
    # look. It is far below anything a plotted trajectory can show, but it means
    # an equality test here would be a flake waiting to happen.
    @test segments[1].times_s[1] == 0.0
    @test segments[1].times_s[2] ≈ 0.01 * 86400 atol = 1e-3
    @test segments[2].times_s[1] > segments[1].times_s[end]

    # Positions are kilometres, as 3-tuples.
    @test segments[1].positions_km[1] == (7000.0, 0.0, 1300.0)
    @test length(segments[1].positions_km) == 4

    # Each arc gets the next colour, so two arcs are visibly different.
    @test segments[1].color != segments[2].color
    @test segments[1].color == IO_.DEFAULT_PALETTE[1]
    @test segments[2].color == IO_.DEFAULT_PALETTE[2]

    # The epoch handed to Cesium is UTC, whatever scale the history was recorded
    # in. This history is TDB, so the string is 69.184 s earlier than the epoch
    # that built it — TT-TAI plus the 37 leap seconds standing in 2020. Cesium
    # reads CZML as UTC and has no idea what TDB is, so skipping this conversion
    # would slide the whole trajectory by over a minute against the rotating
    # Earth, which is about 30 km of ground track.
    @test endswith(epoch_iso, "Z")
    @test startswith(epoch_iso, "2019-12-31T23:58:50")

    # The 1-sample marker is a burn: 0.2 km/s in y, read from the velocity jump
    # between the arcs either side rather than from the marker's own state.
    @test length(maneuvers) == 1
    @test maneuvers[1].name == "TOI"
    @test maneuvers[1].dv_kms[2] ≈ 0.2 rtol = 1e-9
    # The burn is stamped at the marker's own time, which is where the second arc
    # begins rather than where the first one ends.
    @test maneuvers[1].time_s ≈ segments[2].times_s[1] atol = 1e-3
end

@testset "a boundary without a velocity change is not a burn" begin
    # A marker sitting between two arcs that fly at the same speed is a segment
    # boundary, not a maneuver. Drawing a 0 km/s burn there would be noise.
    cs  = CoordinateSystem(earth, ICRF())
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), _OV_EPOCH;
                     coord_sys = cs, name = "Steady")
    for (nm, k0, n) in (("coast 1", 0, 4), ("boundary", 4, 1), ("coast 2", 4, 4))
        s = AstroModels.HistorySegment(cs; name = nm)
        for k in k0:(k0 + n - 1)
            push!(s.times,  Time(2458849.5 + k * 0.01, 0.0, :tdb, :jd))
            push!(s.states, CartesianState([7000.0 + 10k, 100.0k, 1300.0, 0.0, 7.35, 1.0]))
        end
        push!(sat.history.segments, s)
    end

    _, _, maneuvers = IO_._extract(sat)
    @test isempty(maneuvers)
end

@testset "three burns with one name still draw as three points" begin
    # This is the bug. solve_trajectory! names every marker it records
    # "maneuver", so the packet id — which was built from the name alone —
    # was identical for all three. Cesium merges packets by id, so three burns
    # became one point, and maneuver rendering was switched off rather than
    # diagnosed.
    sat = _recorded_three_burns()
    _, segments, maneuvers = IO_._extract(sat)

    @test length(segments) == 4
    @test length(maneuvers) == 3
    @test all(m -> m.name == "maneuver", maneuvers)      # the names really do repeat

    # Each burn is differenced against its own neighbours, so the three sizes
    # come back distinct and in order.
    @test [round(m.dv_kms[2], digits = 6) for m in maneuvers] == [0.2, 0.3, 0.5]

    packets = IO_.build_orbitview(sat)
    ids = [p["id"] for p in packets]
    mnv = filter(x -> occursin("-mnv-", x), ids)

    @test length(mnv) == 3
    @test length(unique(mnv)) == 3                       # the fix: ids are indexed
    @test length(unique(ids)) == length(ids)             # and nothing else collides either
end

@testset "segments recorded in different time scales still draw" begin
    # A spacecraft built in TAI and put through a targeting solve records its first segments in
    # TAI and its propagated ones in TT. `Time - Time` refuses a difference across scales, and
    # orbitview raised on exactly that history. The first arc here is TAI in isot format and the
    # second is TT in jd format, 0.01 day later, so both a scale and a format differ.
    cs  = CoordinateSystem(earth, ICRF())
    t0  = Time("2020-01-01T00:00:00.000", TAI(), ISOT())
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), t0;
                     coord_sys = cs, name = "Mixed")

    tai = AstroModels.HistorySegment(cs; name = "coast 1")
    for k in 0:3
        push!(tai.times,  AstroEpochs._time_jd(t0.jd1, t0.jd2 + k * 0.001, :tai, :isot))
        push!(tai.states, CartesianState([7000.0 + 10k, 100.0k, 1300.0, 0.0, 7.35, 1.0]))
    end
    tt0 = t0.tt
    tt  = AstroModels.HistorySegment(cs; name = "coast 2")
    for k in 0:3
        push!(tt.times,  Time(tt0.jd1, tt0.jd2 + 0.01 + k * 0.001, :tt, :jd))
        push!(tt.states, CartesianState([7100.0 + 10k, 100.0k, 1300.0, 0.0, 7.35, 1.0]))
    end
    push!(sat.history.segments, tai, tt)

    epoch_iso, segments, _ = IO_._extract(sat)
    @test length(segments) == 2
    @test startswith(epoch_iso, "2019-12-31T23:59:23")          # the TAI arc, shown in UTC
    @test segments[2].times_s[1] ≈ 0.01 * 86400 atol = 1e-3      # offset measured in one scale
end

@testset "a spacecraft with nothing recorded says so" begin
    sat = Spacecraft(CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]), _OV_EPOCH;
                     coord_sys = CoordinateSystem(earth, ICRF()), name = "Fresh")
    err = try
        IO_._extract(sat); nothing
    catch e
        e
    end
    # An ArgumentError, because a spacecraft with nothing recorded is a bad argument rather than
    # something that went wrong partway through. The message names the spacecraft and the call
    # that fixes it.
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("history", msg)
    @test occursin("Fresh", msg)
    @test occursin("propagate!", msg)
end

@testset "playback fits the trajectory into about two minutes" begin
    sat = _recorded()

    # A quarter day is ~21600 s, so :auto asks for roughly 21600/120 = 180x.
    auto = IO_.build_orbitview(sat)
    span = 0.07 * 86400                       # earliest to latest sample, seconds
    @test auto[1]["clock"]["multiplier"] ≈ clamp(span / 120.0, 1.0, 86_400.0) rtol = 1e-6

    # An explicit rate is passed through untouched.
    @test IO_.build_orbitview(sat; speed = 3600)[1]["clock"]["multiplier"] == 3600.0

    # Colour and width are the caller's if they want them.
    custom = IO_.build_orbitview(sat; colors = [(1, 2, 3, 4)], linewidth = 4)
    @test custom[2]["path"]["width"] == 4.0
    @test custom[2]["path"]["material"]["solidColor"]["color"]["rgba"] == [1, 2, 3, 4]
end

# ── The public call ──────────────────────────────────────────────────────────
@testset "orbitview puts a scene on the dashboard" begin
    IO_.reset_panels!()
    sat = _recorded()

    orbitview("Trajectory", sat)
    p = IO_.get_panel("Trajectory")
    @test p.family === :scene
    @test p.czml !== nothing
    @test p.czml[1]["id"] == "document"
    @test p.traces == []                       # a scene carries CZML, not traces

    # The payload a scene publishes is its CZML, not a Plotly figure.
    @test JSON.parse(IO_.payload(p))[1]["id"] == "document"

    # And the manifest tells the page it is a scene, so it loads Cesium rather
    # than handing the packets to Plotly.
    @test any(e -> e[:title] == "Trajectory" && e[:kind] == "scene", IO_.manifest()[:panels])

    # Without a name it still lands somewhere.
    orbitview(sat)
    @test any(p -> p.family === :scene, values(IO_.panels()))
end

@testset "orbitview says what it cannot do yet" begin
    IO_.reset_panels!()
    sat = _recorded()

    # Two spacecraft in one view needs the writer to merge entity sets. It does
    # not, so it says so rather than drawing one and dropping the other.
    err = try
        orbitview!("Trajectory", sat); nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("own panel", sprint(showerror, err))

    # The wrong number of arguments names the two shapes that work.
    err2 = try
        orbitview("Trajectory", sat, sat); nothing
    catch e
        e
    end
    @test err2 isa ArgumentError
    @test occursin("orbitview(sat)", sprint(showerror, err2))
end

end # testset
