# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Leap seconds: TAI − UTC from the IERS list, kept current.
#
# Truth, by section:
#   the list in use        the IERS list as IANA publishes it: the last change is 2017-01-01 to
#                          37 s, and the list has not expired. Needs the network the first time.
#   the built-in table     the IERS list again: the table used offline, adapted from Tempo.jl, has
#                          every change the list has through 2017-01-01, and every date either side
#                          of every change reads the same value from both.
#   the change itself      the definition: TAI − UTC steps at 0h UTC on the date in the list.
#   a future leap second   a list with one more entry, parsed from text: dates after it read the
#                          new value, which is the reason for downloading the list at all.

using Test
using AstroEpochs
using AstroEpochs: leap_second_table, tai_minus_utc, parse_leap_seconds,
                   utctai, taiutc, LeapSecondTable, _LEAP, _builtin_leap_table

const _J2000 = 2451545.0

@testset "leap seconds — the list in use is current" begin
    table = leap_second_table()
    @test table.delta[end] == 37.0
    @test table.starts[end] ≈ 2457754.5 - _J2000              # 2017-01-01T00:00:00 UTC
    @test table.expires > time() / 86_400 + 2440587.5 - _J2000   # has not expired
    @test issorted(table.starts)
end

@testset "leap seconds — the built-in table agrees with the IERS list" begin
    builtin = _builtin_leap_table()
    list    = leap_second_table()
    n       = length(builtin.starts)
    @test builtin.starts == list.starts[1:n]
    @test builtin.delta  == list.delta[1:n]

    # Either side of every change, through the transforms a conversion uses (ERFA's utctai and
    # taiutc): 0h on the day that ends with the leap second, and a quarter-day into the next. At 0h
    # the UTC and TAI dates differ by the table's value exactly. Later on a leap-second day they
    # differ by more, because that day's fraction runs over 86401 s; that is ERFA's quasi-JD UTC,
    # held to Astropy in test_correctness_astropy_benchmark.jl.
    for (i, d) in enumerate(builtin.starts), (days, idx) in ((-1.0, i - 1), (0.25, i))
        idx < 1 && continue                     # before 1972 there is no value, and it warns
        expected = builtin.delta[idx]
        utc = (_J2000 + d + days, 0.0)            # the date in the large part, as ERFA expects
        @test tai_minus_utc(d + days) == expected
        tai = utctai(utc...)
        @test ((tai[1] - utc[1]) + (tai[2] - utc[2])) * 86_400 ≈ expected atol = 1e-9
        back = taiutc(tai...)
        @test ((back[1] - utc[1]) + (back[2] - utc[2])) * 86_400 ≈ 0.0 atol = 1e-9
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
