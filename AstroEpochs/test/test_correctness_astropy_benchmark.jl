# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# AstroEpochs against Astropy, across the scales and formats, at 1668 epochs.
#
# Truth: Astropy 8.0.1 (pyerfa 2.0.1.5, ERFA 2.0.1), in test/astropy/reference.csv, written by
# test/astropy/make_reference.py; this file never runs Python. For each case the file has an input
# epoch in one scale and format, and Astropy's two-part JD and ISOT string for it in every scale.
#
# AstroEpochs uses Astropy's algorithms (ERFA's scale transforms, TDB − TT series, calendar and
# formatting, and Astropy's normalisation of the two-part date), so it agrees to the last bit of
# the date: the tolerance, 2e-11 s, is two units in the last place of a date's fraction of a day,
# and was measured at zero. ISOT strings agree exactly.
#
# Case kinds (see make_reference.py): random epochs 1972-2100 in every scale; either side of every
# leap second in UTC and TAI; instants inside a leap second (UTC 23:59:60.x); ISOT strings that
# round across a second, minute, hour or day; J2000, the 1977 TCG/TCB epoch, MJD 0 and 1600-2500.
#
# Known gap, by decision: UTC before 1972, where Astropy applies the drifting pre-1972 offsets and
# AstroEpochs uses TAI − UTC = 0 (with a warning). Those cases are `@test_broken`, so they show in
# the summary and fail loudly if the gap is ever closed without this file being updated.
#
# Round trips: every case through every pair of scales and back. A route of k hops rounds up to 2k
# times, and the longest (TCB to UTC) has four each way, so a round trip is held to 1e-10 s, about
# ten units in the last place; the worst measured was 2.9e-11 s.
#
# Epochs before 1972 warn on each UTC conversion (see leap_seconds.jl). That warning is tested in
# the known-gap testset; the two large testsets run with it silenced.

using Test
using AstroEpochs
using Logging: with_logger, NullLogger

const _TOL_S = 2e-11
const _TOL_ROUND_TRIP_S = 1e-10
const _SCALES = (:tai, :tt, :utc, :tdb, :tcb, :tcg)
const _TAGS = Dict(:tai => TAI(), :tt => TT(), :utc => UTC(), :tdb => TDB(), :tcb => TCB(),
                   :tcg => TCG())

function _read_reference()
    lines = filter(l -> !startswith(l, "#"), readlines(joinpath(@__DIR__, "astropy", "reference.csv")))
    hdr = split(lines[1], ",")
    return [Dict(zip(hdr, split(l, ","))) for l in lines[2:end]]
end

_input_time(r) = r["in_format"] == "jd" ?
    Time(parse(Float64, r["in_value"]), parse(Float64, r["in_value2"]), _TAGS[Symbol(r["in_scale"])], JD()) :
    Time(String(r["in_value"]), _TAGS[Symbol(r["in_scale"])], ISOT())

_err_s(t, jd1, jd2) = ((t.jd1 - jd1) + (t.jd2 - jd2)) * 86400

const _REFERENCE = _read_reference()

with_logger(NullLogger()) do
@testset "Astropy benchmark — every scale, every case kind" begin
    for kind in ("random", "leap", "leap_inst", "rounding", "special")
        @testset "$kind" begin
            for r in filter(r -> r["kind"] == kind, _REFERENCE), out in _SCALES
                t = getproperty(_input_time(r), out)
                err = _err_s(t, parse(Float64, r["$(out)_jd1"]), parse(Float64, r["$(out)_jd2"]))
                @test abs(err) <= _TOL_S
                @test t.isot == r["$(out)_isot"]
            end
        end
    end
end
end

@testset "Astropy benchmark — UTC before 1972 (known gap)" begin
    # AstroEpochs uses TAI − UTC = 0 before 1972 and warns; Astropy uses the pre-1972 drift. The
    # difference measured here is 9.9 s. Marked broken rather than skipped, so it stays visible.
    for r in filter(r -> r["kind"] == "pre1972", _REFERENCE)
        t = (@test_logs (:warn,) match_mode = :any getproperty(_input_time(r), :tai))
        err = _err_s(t, parse(Float64, r["tai_jd1"]), parse(Float64, r["tai_jd2"]))
        @test abs(err) < 11.0                   # the gap is the pre-1972 offset, and no more
        @test_broken abs(err) <= _TOL_S
    end
end

with_logger(NullLogger()) do
@testset "Round trips — every case, every pair of scales" begin
    for r in filter(r -> r["kind"] != "pre1972", _REFERENCE)
        t0 = _input_time(r)
        for a in _SCALES, b in _SCALES
            a == b && continue
            ta = getproperty(t0, a)
            back = getproperty(getproperty(ta, b), a)
            @test abs(_err_s(back, ta.jd1, ta.jd2)) <= _TOL_ROUND_TRIP_S
        end
    end
end
end
