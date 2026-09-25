# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Leap seconds: TAI − UTC from the IERS list, kept current.
#
# Truth, by section:
#   the list in use        the IERS list as IANA publishes it: the last change is 2017-01-01 to
#                          37 s, and the list has not expired. Needs the network the first time.
#   agreement with Tempo   Tempo.jl's built-in table, which is correct through 2017-01-01; every
#                          date either side of every change must agree.
#   the change itself      the definition: TAI − UTC steps at 0h UTC on the date in the list.
#   a future leap second   a list with one more entry, parsed from text: dates after it read the
#                          new value, which is the reason for downloading the list at all.

using Test
using AstroEpochs
using AstroEpochs: Tempo, leap_second_table, tai_minus_utc, parse_leap_seconds,
                   offset_utc2tai, offset_tai2utc, LeapSecondTable, _LEAP

const _J2000 = 2451545.0

@testset "leap seconds — the list in use is current" begin
    table = leap_second_table()
    @test table.delta[end] == 37.0
    @test table.starts[end] ≈ 2457754.5 - _J2000              # 2017-01-01T00:00:00 UTC
    @test table.expires > time() / 86_400 + 2440587.5 - _J2000   # has not expired
    @test issorted(table.starts)
end

@testset "leap seconds — agree with Tempo either side of every change" begin
    starts = Tempo.LEAPSECONDS.jd2000
    # The table lookup agrees exactly. Tempo's offset functions pass the epoch through one
    # floating-point Julian date and back, which costs up to about 1e-7 s; ours return the
    # integer offset, so the offsets agree to that rounding and no further.
    for d in starts, side in (-0.5, 0.5)
        utc_days = d + side
        @test tai_minus_utc(utc_days) == Tempo.leapseconds(utc_days)
        @test offset_utc2tai(utc_days * 86_400) ≈ Tempo.offset_utc2tai(utc_days * 86_400) atol = 1e-6
        tai_seconds = utc_days * 86_400 + tai_minus_utc(utc_days)
        @test offset_tai2utc(tai_seconds) ≈ Tempo.offset_tai2utc(tai_seconds) atol = 1e-6
    end
end

@testset "leap seconds — TAI − UTC steps at 0h UTC on the date" begin
    new_year_2017 = 2457754.5 - _J2000
    @test tai_minus_utc(new_year_2017 - 1 / 86_400) == 36.0      # 2016-12-31T23:59:59 UTC
    @test tai_minus_utc(new_year_2017) == 37.0                   # 2017-01-01T00:00:00 UTC

    # Through Time: UTC to TAI adds 37 s today, and the round trip is exact
    utc = Time("2024-06-01T12:00:00", UTC(), ISOT())
    @test (utc.tai - AstroEpochs._time_jd(utc.jd1, utc.jd2, :tai, :jd)) * 86_400 ≈ 37.0 atol = 1e-6
    @test utc.tai.utc.isot == utc.isot
end

@testset "leap seconds — a leap second added to the list takes effect" begin
    # The published list's format, with one invented entry: 2030-01-01 (NTP 4102444800) to 38 s
    text = """
    # a comment line
    #@	4133980800
    2272060800	10	# 1 Jan 1972
    3692217600	37	# 1 Jan 2017
    4102444800	38	# 1 Jan 2030, invented for this test
    """
    table = parse_leap_seconds(text)
    @test length(table.delta) == 3 && table.delta[end] == 38.0

    saved = _LEAP[]
    try
        _LEAP[] = table
        new_year_2030 = 2462502.5 - _J2000                          # 2030-01-01T00:00:00 UTC
        @test tai_minus_utc(new_year_2030 - 1 / 86_400) == 37.0
        @test tai_minus_utc(new_year_2030) == 38.0
        utc = Time("2030-06-01T00:00:00", UTC(), ISOT())
        @test (utc.tai - AstroEpochs._time_jd(utc.jd1, utc.jd2, :tai, :jd)) * 86_400 ≈ 38.0 atol = 1e-6
    finally
        _LEAP[] = saved
    end

    @test_throws ArgumentError parse_leap_seconds("# no entries\n#@\t4133980800\n")
    @test_throws ArgumentError parse_leap_seconds("2272060800\t10\n")          # no expiry line
end

nothing
