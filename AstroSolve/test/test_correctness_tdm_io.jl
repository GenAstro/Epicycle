# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# CCSDS 503.0-B-2 tracking data message reading and writing.
#
# Truth is a round trip: records written and read back must be the records that went in. That is
# an exact property rather than an approximate one, which makes it the strongest truth available
# for a file format short of an external reference file — and unlike the estimator tests it needs
# no physics, no propagation and no optimizer.
#
# The round trip is not quite the identity, and the two places it is not are contract, so they are
# asserted rather than tolerated. Epochs are written to millisecond precision, so they come back
# on the millisecond grid. `computed` is not a TDM data field, so it comes back `NaN` whatever went
# in. Everything else must survive exactly, including which station each record came from, which is
# what lets an estimator pick the matching measurement specification.

using AstroSolve
using AstroEpochs
using Test

const _TDM = AstroSolve.TrackingDataIO

using AstroSolve.TrackingDataIO: ObservationRecord, TrackingDataFile, TDMHeader,
                                 TDMSegmentMeta, CCSDS_KVN, read_records, write_records,
                                 SUPPORTED_OBSERVABLES, ALLOWED_MODES, ALLOWED_TIME_SYSTEMS

const _TDM_DIR = mktempdir()
_tdm_path(name) = joinpath(_TDM_DIR, name)

const _TDM_HDR = TDMHeader("2020-03-01T12:00:00", "GEN ASTRO")

_tdm_meta(p1; scale = "TT", p2 = "SAT-1") =
    TDMSegmentMeta(time_system = scale, participant_1 = p1, participant_2 = p2)

# Epochs survive to the millisecond, which is what the writer formats to. One millisecond in days.
const _TDM_MS = 1 / 86_400_000

_tdm_epoch(iso; scale = TT()) = Time(iso, scale, ISOT())

# Compare epochs by absolute Julian date, which sidesteps the format tag entirely.
_tdm_close(a::Time, b::Time; tol = _TDM_MS) = abs(a.jd - b.jd) <= tol

"""Write `text` to a scratch file and return the `TrackingDataFile` for it."""
function _tdm_file(name, text)
    path = _tdm_path(name)
    open(io -> print(io, text), path, "w")
    return TrackingDataFile(path)
end

@testset "TDM — a single segment round-trips" begin
    epochs = [_tdm_epoch("2020-03-01T00:0$(i):00.000") for i in 0:4]
    recs = [ObservationRecord(:RANGE, epochs[i], 42164.0 + i, "DSS-14") for i in 1:5]
    file = TrackingDataFile(_tdm_path("single.tdm"))

    # The writer hands back the file it was given, so a write composes with a read.
    @test write_records(file, recs, _TDM_HDR, _tdm_meta("DSS-14")) === file

    back, header, metas = read_records(file)

    @test length(back) == 5
    @test length(metas) == 1

    # Header keywords survive verbatim.
    @test header.creation_date == _TDM_HDR.creation_date
    @test header.originator    == _TDM_HDR.originator

    # Metadata survives, including the fields the caller did not set and took as defaults.
    @test metas[1].time_system   == "TT"
    @test metas[1].participant_1 == "DSS-14"
    @test metas[1].participant_2 == "SAT-1"
    @test metas[1].mode          == "SEQUENTIAL"
    @test metas[1].path          == "1,2,1"
    @test metas[1].range_units   == "km"

    for (r, orig) in zip(back, recs)
        @test r.measurement_type == :RANGE
        @test r.observed == orig.observed          # exact: Float64 printing round-trips
        @test r.participant_1 == "DSS-14"          # tagged from the enclosing segment
        @test _tdm_close(r.t_receive, orig.t_receive)
        @test isnan(r.computed)                    # not a TDM data field
    end
end

@testset "TDM — several stations round-trip into their own segments" begin
    # One file, two passes. This is the case the participant_1 tagging exists for: an estimator
    # has to know which station each measurement came from to pick its measurement spec.
    e(i) = _tdm_epoch("2020-03-01T00:0$(i):00.000")
    recs = ObservationRecord[]
    for i in 0:3
        push!(recs, ObservationRecord(:RANGE,   e(i), 1000.0 + i, "DSS-14"))
        push!(recs, ObservationRecord(:DOPPLER, e(i), -2.0 - i,   "DSS-43"))
    end

    metas = [_tdm_meta("DSS-14"), _tdm_meta("DSS-43")]
    file  = TrackingDataFile(_tdm_path("two_station.tdm"))
    write_records(file, recs, _TDM_HDR, metas)

    back, _, back_metas = read_records(file)

    @test length(back) == 8
    @test length(back_metas) == 2
    @test [m.participant_1 for m in back_metas] == ["DSS-14", "DSS-43"]

    # Records arrive grouped by segment rather than interleaved as they went in, which is what
    # the reader's docstring says and why build_od_closures sorts by epoch.
    @test [r.participant_1 for r in back] ==
          ["DSS-14", "DSS-14", "DSS-14", "DSS-14", "DSS-43", "DSS-43", "DSS-43", "DSS-43"]

    # Each station keeps its own observable and its own values.
    goldstone = filter(r -> r.participant_1 == "DSS-14", back)
    canberra  = filter(r -> r.participant_1 == "DSS-43", back)
    @test all(r -> r.measurement_type === :RANGE,   goldstone)
    @test all(r -> r.measurement_type === :DOPPLER, canberra)
    @test [r.observed for r in goldstone] == [1000.0, 1001.0, 1002.0, 1003.0]
    @test [r.observed for r in canberra]  == [-2.0, -3.0, -4.0, -5.0]
end

@testset "TDM — the writer orders records and derives the segment span" begin
    # Records are sorted by epoch within a segment, and START_TIME and STOP_TIME come from the
    # records rather than from the meta, so they always describe what was written.
    shuffled = [_tdm_epoch("2020-03-01T00:0$(i):00.000") for i in (3, 0, 4, 1, 2)]
    recs = [ObservationRecord(:RANGE, shuffled[i], Float64(i), "DSS-14") for i in 1:5]

    file = TrackingDataFile(_tdm_path("ordering.tdm"))
    write_records(file, recs, _TDM_HDR, _tdm_meta("DSS-14"))
    text = read(file.path, String)

    @test occursin("START_TIME    = 2020-03-01T00:00:00.000", text)
    @test occursin("STOP_TIME     = 2020-03-01T00:04:00.000", text)

    back, _, _ = read_records(file)
    @test issorted([r.t_receive.jd for r in back])

    # Sorting moves the values with the epochs rather than only the epochs.
    @test [r.observed for r in back] == [2.0, 4.0, 5.0, 1.0, 3.0]
end

@testset "TDM — an untagged record belongs to the only segment" begin
    # Callers that predate participant_1 leave it empty. With one meta there is no ambiguity, so
    # the record is routed to it and comes back tagged.
    t = _tdm_epoch("2020-03-01T00:00:00.000")
    recs = [ObservationRecord(:RANGE, t, 7000.0, "")]
    file = TrackingDataFile(_tdm_path("untagged.tdm"))
    write_records(file, recs, _TDM_HDR, _tdm_meta("DSS-14"))

    back, _, _ = read_records(file)
    @test length(back) == 1
    @test back[1].participant_1 == "DSS-14"

    # With two metas the same record has no segment to go to, and says so rather than picking one.
    file2 = TrackingDataFile(_tdm_path("untagged_two.tdm"))
    err = try
        write_records(file2, recs, _TDM_HDR, [_tdm_meta("DSS-14"), _tdm_meta("DSS-43")])
    catch e
        e
    end
    @test err isa ArgumentError
    msg = sprint(showerror, err)
    @test occursin("participant_1", msg) && occursin("2 metas", msg)
end

@testset "TDM — every supported time scale round-trips in its own scale" begin
    # The epoch is parsed in the scale its segment declares. A file whose segments declare
    # different scales must read each in its own, which a single global scale would break.
    for sys in ALLOWED_TIME_SYSTEMS
        tag  = _TDM._scale_tag(sys)
        t    = Time("2020-03-01T00:00:00.000", tag, ISOT())
        recs = [ObservationRecord(:RANGE, t, 1234.5, "DSS-14")]
        file = TrackingDataFile(_tdm_path("scale_$(sys).tdm"))
        write_records(file, recs, _TDM_HDR, _tdm_meta("DSS-14"; scale = sys))

        back, _, metas = read_records(file)
        @test metas[1].time_system == sys
        @test back[1].t_receive.scale == t.scale
        @test _tdm_close(back[1].t_receive, t)
    end

    # Two segments, two scales, one file. Both epochs are the same instant written two ways, so
    # after reading they must still be the same instant.
    t_tai = Time("2020-03-01T00:00:00.000", TAI(), ISOT())
    t_tt  = t_tai.tt
    recs = [ObservationRecord(:RANGE, t_tai, 1.0, "DSS-14"),
            ObservationRecord(:RANGE, t_tt,  2.0, "DSS-43")]
    file = TrackingDataFile(_tdm_path("mixed_scales.tdm"))
    write_records(file, recs, _TDM_HDR,
                  [_tdm_meta("DSS-14"; scale = "TAI"), _tdm_meta("DSS-43"; scale = "TT")])

    back, _, metas = read_records(file)
    @test [m.time_system for m in metas] == ["TAI", "TT"]
    @test back[1].t_receive.scale == :tai
    @test back[2].t_receive.scale == :tt
    @test abs(back[1].t_receive.jd - back[2].t_receive.tai.jd) <= _TDM_MS
end

@testset "TDM — epochs format on the millisecond grid, including the carry" begin
    # _format_isot rounds to the millisecond and carries 60.000 s into the next minute, which
    # AstroEpochs' own _to_isot does not. An epoch a hair under a minute boundary is the case
    # that distinguishes them, and an emitted `:60.000` is not a legal ISO 8601 second field.
    base = Time("2020-03-01T00:00:59.9999", TT(), ISOT())
    @test _TDM._format_isot(base) == "2020-03-01T00:01:00.000"

    # Midnight carries the civil day as well as the minute.
    eod = Time("2020-03-01T23:59:59.9999", TT(), ISOT())
    @test _TDM._format_isot(eod) == "2020-03-02T00:00:00.000"

    # Nothing else moves: an epoch already on the grid formats to itself.
    on_grid = Time("2020-03-01T12:34:56.789", TT(), ISOT())
    @test _TDM._format_isot(on_grid) == "2020-03-01T12:34:56.789"

    # And the carry survives a round trip rather than only a format call.
    recs = [ObservationRecord(:RANGE, base, 500.0, "DSS-14")]
    file = TrackingDataFile(_tdm_path("carry.tdm"))
    write_records(file, recs, _TDM_HDR, _tdm_meta("DSS-14"))
    back, _, _ = read_records(file)
    @test _tdm_close(back[1].t_receive, base)
end

@testset "TDM — the writer refuses what it cannot write" begin
    # An ArgumentError naming the field, the constraint and the
    # value, with the valid options enumerated where the set is small and closed.
    t    = _tdm_epoch("2020-03-01T00:00:00.000")
    recs = [ObservationRecord(:RANGE, t, 1.0, "DSS-14")]
    meta = _tdm_meta("DSS-14")
    file = TrackingDataFile(_tdm_path("reject.tdm"))

    msg_of(f) = try; f(); ""; catch e; sprint(showerror, e); end

    @test_throws ArgumentError write_records(file, ObservationRecord[], _TDM_HDR, [meta])
    @test occursin("at least one observation",
                   msg_of(() -> write_records(file, ObservationRecord[], _TDM_HDR, [meta])))

    @test_throws ArgumentError write_records(file, recs, _TDM_HDR, TDMSegmentMeta[])

    # An unsupported time system enumerates the ones that work.
    bad_scale = TDMSegmentMeta(time_system = "GPS", participant_1 = "DSS-14",
                               participant_2 = "SAT-1")
    m = msg_of(() -> write_records(file, recs, _TDM_HDR, [bad_scale]))
    @test occursin("TIME_SYSTEM", m) && occursin("GPS", m) && occursin("UTC", m)

    bad_mode = TDMSegmentMeta(time_system = "TT", participant_1 = "DSS-14",
                              participant_2 = "SAT-1", mode = "CONTINUOUS")
    m = msg_of(() -> write_records(file, recs, _TDM_HDR, [bad_mode]))
    @test occursin("MODE", m) && occursin("SEQUENTIAL", m) && occursin("CONTINUOUS", m)

    # Two segments claiming the same station would make routing ambiguous.
    m = msg_of(() -> write_records(file, recs, _TDM_HDR, [meta, _tdm_meta("DSS-14")]))
    @test occursin("distinct participant_1", m) && occursin("DSS-14", m)

    # An observable the format subset does not carry.
    odd = [ObservationRecord(:ANGLE_1, t, 1.0, "DSS-14")]
    m = msg_of(() -> write_records(file, odd, _TDM_HDR, [meta]))
    @test occursin("measurement_type", m) && occursin("ANGLE_1", m) && occursin("RANGE", m)

    # A record whose station no segment claims. The message names the stations that exist, which
    # is what turns a typo into a one-look fix.
    stray = [ObservationRecord(:RANGE, t, 1.0, "DSS-99")]
    m = msg_of(() -> write_records(file, stray, _TDM_HDR, [meta]))
    @test occursin("DSS-99", m) && occursin("DSS-14", m)

    # And a meta set that no record matches at all.
    m = msg_of(() -> write_records(file, recs, _TDM_HDR, [_tdm_meta("DSS-43")]))
    @test occursin("DSS-43", m)
end

@testset "TDM — the reader refuses a malformed file" begin
    # Each branch gets a file that is well-formed except for the one thing under test, so a
    # failure names the defect rather than whichever check happens to fire first.
    good_meta = """
    META_START
    TIME_SYSTEM   = TT
    START_TIME    = 2020-03-01T00:00:00.000
    STOP_TIME     = 2020-03-01T00:00:00.000
    PARTICIPANT_1 = DSS-14
    PARTICIPANT_2 = SAT-1
    MODE          = SEQUENTIAL
    PATH          = 1,2,1
    RANGE_UNITS   = km
    META_STOP
    """
    hdr(v = "2.0") = "CCSDS_TDM_VERS = $v\nCREATION_DATE = 2020-03-01T12:00:00\nORIGINATOR    = GEN ASTRO\n\n"
    data = "DATA_START\nRANGE = 2020-03-01T00:00:00.000 42164.0\nDATA_STOP\n"
    whole = hdr() * good_meta * "\n" * data

    msg_of(f) = try; f(); ""; catch e; sprint(showerror, e); end
    reads(name, text) = read_records(_tdm_file(name, text))

    # The baseline parses, which is what makes each mutation below a controlled change.
    recs, _, _ = reads("ok.tdm", whole)
    @test length(recs) == 1

    # A file that is not there.
    m = msg_of(() -> read_records(TrackingDataFile(_tdm_path("absent.tdm"))))
    @test occursin("readable file", m) && occursin("absent.tdm", m)

    # Wrong or missing version.
    m = msg_of(() -> reads("badvers.tdm", hdr("1.0") * good_meta * "\n" * data))
    @test occursin("CCSDS_TDM_VERS", m) && occursin("2.0", m)
    m = msg_of(() -> reads("novers.tdm",
                           "CREATION_DATE = x\nORIGINATOR = y\n\n" * good_meta * "\n" * data))
    @test occursin("CCSDS_TDM_VERS", m)

    # A META segment missing a required key enumerates what is required.
    stripped = replace(good_meta, "PATH          = 1,2,1\n" => "")
    m = msg_of(() -> reads("nopath.tdm", hdr() * stripped * "\n" * data))
    @test occursin("PATH", m) && occursin("RANGE_UNITS", m)

    # Unsupported MODE and TIME_SYSTEM, each enumerating the supported set.
    m = msg_of(() -> reads("badmode.tdm",
                           hdr() * replace(good_meta, "SEQUENTIAL" => "CONTINUOUS") *
                           "\n" * data))
    @test occursin("MODE", m) && occursin("SEQUENTIAL", m)
    m = msg_of(() -> reads("badsys.tdm",
                           hdr() * replace(good_meta, "TIME_SYSTEM   = TT" =>
                                                      "TIME_SYSTEM   = GPS") * "\n" * data))
    @test occursin("TIME_SYSTEM", m) && occursin("GPS", m)

    # An observable outside the subset.
    m = msg_of(() -> reads("badobs.tdm",
                           hdr() * good_meta * "\nDATA_START\n" *
                           "ANGLE_1 = 2020-03-01T00:00:00.000 1.0\nDATA_STOP\n"))
    @test occursin("observable", m) && occursin("ANGLE_1", m)

    # A data line missing its value.
    m = msg_of(() -> reads("shortline.tdm",
                           hdr() * good_meta * "\nDATA_START\n" *
                           "RANGE = 2020-03-01T00:00:00.000\nDATA_STOP\n"))
    @test occursin("<epoch> <value>", m)

    # A metadata line with no '=' at all.
    m = msg_of(() -> reads("nokv.tdm", hdr() * "META_START\nTIME_SYSTEM TT\nMETA_STOP\n"))
    @test occursin("KEY = VALUE", m)

    # Segment nesting and termination, in both directions.
    @test occursin("still open",
                   msg_of(() -> reads("nested_meta.tdm",
                                      hdr() * "META_START\nMETA_START\n")))
    @test occursin("META_START",
                   msg_of(() -> reads("stray_metastop.tdm", hdr() * "META_STOP\n")))
    @test occursin("DATA_START",
                   msg_of(() -> reads("stray_datastop.tdm",
                                      hdr() * good_meta * "\nDATA_STOP\n")))
    @test occursin("META segment",
                   msg_of(() -> reads("data_first.tdm", hdr() * "DATA_START\n")))
    @test occursin("ends inside a META segment",
                   msg_of(() -> reads("open_meta.tdm", hdr() * "META_START\nTIME_SYSTEM = TT\n")))
    @test occursin("ends inside a DATA segment",
                   msg_of(() -> reads("open_data.tdm", hdr() * good_meta * "\nDATA_START\n" *
                                      "RANGE = 2020-03-01T00:00:00.000 1.0\n")))

    # A header key the format does not define, outside any segment. COMMENT is not that case:
    # CCSDS allows it anywhere and the reader skips it, which the next testset covers.
    m = msg_of(() -> reads("strayhdr.tdm", hdr() * "ORIGIN = hello\n" * good_meta * "\n" * data))
    @test occursin("CCSDS_TDM_VERS", m) && occursin("ORIGIN", m)

    # A file with a header and nothing else.
    @test occursin("META segment", msg_of(() -> reads("hdronly.tdm", hdr())))

    # A file with metadata but no observations.
    @test occursin("at least one observation",
                   msg_of(() -> reads("nodata.tdm", hdr() * good_meta * "\n")))
end

@testset "TDM — the constructors normalize what they are given" begin
    # ObservationRecord has a four-argument form that leaves `computed` unset and a five-argument
    # form that takes it. Both accept any Real and store Float64, so an Int observation from a
    # caller does not change the record's type.
    t = _tdm_epoch("2020-03-01T00:00:00.000")

    four = ObservationRecord(:RANGE, t, 42164, "DSS-14")
    @test four.observed === 42164.0
    @test isnan(four.computed)

    five = ObservationRecord(:RANGE, t, 42164, 42163, "DSS-14")
    @test five.observed === 42164.0
    @test five.computed === 42163.0
    @test five.participant_1 == "DSS-14"

    # TrackingDataFile defaults to the KVN format and accepts one explicitly. Both store the path
    # as a String, so a SubString from a path split does not leak into the type.
    implied = TrackingDataFile(_tdm_path("x.tdm"))
    @test implied.format isa CCSDS_KVN
    @test implied.path isa String

    # A String plus a format hits the struct's own constructor; any other AbstractString goes
    # through the converting method, which is the case a path from `split` or `chop` lands in.
    explicit = TrackingDataFile(_tdm_path("x.tdm"), CCSDS_KVN())
    @test explicit.format isa CCSDS_KVN
    @test explicit.path == implied.path

    sub = TrackingDataFile(SubString(_tdm_path("x.tdm"), 1), CCSDS_KVN())
    @test sub.path isa String
    @test sub.path == implied.path
end

@testset "TDM — COMMENT lines are ignored wherever they appear" begin
    # CCSDS allows COMMENT anywhere in a TDM. The reader skips them, so a commented file and the
    # same file without comments must parse to the same records; otherwise a perfectly legal file
    # written by another tool fails to load.
    head = "CCSDS_TDM_VERS = 2.0\nCREATION_DATE = 2020-03-01T12:00:00\nORIGINATOR    = GEN ASTRO\n\n"
    meta = """
    META_START
    TIME_SYSTEM   = TT
    START_TIME    = 2020-03-01T00:00:00.000
    STOP_TIME     = 2020-03-01T00:00:00.000
    PARTICIPANT_1 = DSS-14
    PARTICIPANT_2 = SAT-1
    MODE          = SEQUENTIAL
    PATH          = 1,2,1
    RANGE_UNITS   = km
    META_STOP
    """
    data = "DATA_START\nRANGE = 2020-03-01T00:00:00.000 42164.0\nDATA_STOP\n"

    plain, _, _ = read_records(_tdm_file("nocomment.tdm", head * meta * "\n" * data))

    commented = head * "COMMENT top-level note\n" *
                replace(meta, "META_START\n" => "META_START\n    COMMENT inside metadata\n") *
                "\n" * replace(data, "DATA_START\n" => "DATA_START\nCOMMENT inside data\n")
    withcom, _, metas = read_records(_tdm_file("comment.tdm", commented))

    @test length(withcom) == length(plain) == 1
    @test withcom[1].observed == plain[1].observed
    @test withcom[1].t_receive.jd == plain[1].t_receive.jd
    @test metas[1].participant_1 == "DSS-14"
end

@testset "TDM — _scale_tag maps every allowed system and refuses the rest" begin
    # The mapping is a small closed set, so the message enumerates it rather than sending the
    # reader to the docstring.
    @test _TDM._scale_tag("UTC") == UTC()
    @test _TDM._scale_tag("TAI") == TAI()
    @test _TDM._scale_tag("TT")  == TT()
    @test _TDM._scale_tag("TDB") == TDB()
    @test _TDM._scale_tag("TCG") == TCG()
    @test _TDM._scale_tag("TCB") == TCB()

    # Case and surrounding space are not the caller's problem.
    @test _TDM._scale_tag("  tt  ") == TT()

    @test_throws ArgumentError _TDM._scale_tag("GPS")
    m = try; _TDM._scale_tag("GPS"); catch e; sprint(showerror, e); end
    @test occursin("TIME_SYSTEM", m) && occursin("GPS", m) && occursin("TDB", m)
end
