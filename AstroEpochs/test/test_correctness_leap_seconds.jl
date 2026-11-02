# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Leap seconds: TAI − UTC from the IERS list, kept current.
#
# Truth, by section:
#   the list in use        the IERS changes through 2017-01-01, which every list holds, whether
#                          it came from the network, a stored copy or the built-in table. A later
#                          leap second or an expired copy does not fail this; no test needs the
#                          network.
#   the built-in table     the IERS list again: the table used offline, adapted from Tempo.jl, has
#                          every change the list has through 2017-01-01, and every date either side
#                          of every change reads the same value from both.
#   the change itself      the definition: TAI − UTC steps at 0h UTC on the date in the list.
#   a future leap second   a list with one more entry, parsed from text: dates after it read the
#                          new value, which is the reason for downloading the list at all.
#   a damaged list         the published format: the hash on its last line, and the changes every
#                          list holds. A real list, as IANA published it, is the one that must pass.

using Test
using AstroEpochs
using AstroEpochs: leap_second_table, tai_minus_utc, parse_leap_seconds,
                   utctai, taiutc, LeapSecondTable, _LEAP, _builtin_leap_table,
                   _BUILTIN_LEAP_SECONDS, cal2jd

include("leap_list_helpers.jl")

const _J2000 = 2451545.0

@testset "leap seconds — the list in use holds every change through 2017" begin
    table = leap_second_table()
    i2017 = findfirst(==(2457754.5 - _J2000), table.starts)     # 2017-01-01T00:00:00 UTC
    @test i2017 !== nothing
    @test table.delta[i2017] == 37.0
    @test issorted(table.starts)
    @test all(==(1.0), diff(table.delta[2:end]))                # one second at a time, since 1972
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
    text = _valid_leap_list(4133980800; extra = [(4102444800, 38)])
    table = parse_leap_seconds(text)
    @test length(table.delta) == 29 && table.delta[end] == 38.0

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

end

@testset "leap seconds — the hash is computed as the published file computes it" begin
    # The head and tail of leap-seconds.list as IANA published it, whose hash the IERS computed.
    published = """
    #\$\t3992312697
    #@\t4023129600
    2272060800\t10\t# 1 Jan 1972
    2287785600\t11\t# 1 Jul 1972
    2303683200\t12\t# 1 Jan 1973
    2335219200\t13\t# 1 Jan 1974
    2366755200\t14\t# 1 Jan 1975
    2398291200\t15\t# 1 Jan 1976
    2429913600\t16\t# 1 Jan 1977
    2461449600\t17\t# 1 Jan 1978
    2492985600\t18\t# 1 Jan 1979
    2524521600\t19\t# 1 Jan 1980
    2571782400\t20\t# 1 Jul 1981
    2603318400\t21\t# 1 Jul 1982
    2634854400\t22\t# 1 Jul 1983
    2698012800\t23\t# 1 Jul 1985
    2776982400\t24\t# 1 Jan 1988
    2840140800\t25\t# 1 Jan 1990
    2871676800\t26\t# 1 Jan 1991
    2918937600\t27\t# 1 Jul 1992
    2950473600\t28\t# 1 Jul 1993
    2982009600\t29\t# 1 Jul 1994
    3029443200\t30\t# 1 Jan 1996
    3076704000\t31\t# 1 Jul 1997
    3124137600\t32\t# 1 Jan 1999
    3345062400\t33\t# 1 Jan 2006
    3439756800\t34\t# 1 Jan 2009
    3550089600\t35\t# 1 Jul 2012
    3644697600\t36\t# 1 Jul 2015
    3692217600\t37\t# 1 Jan 2017
    #h\ta9bad145 84c31c70 758402aa b37bfd54 5923836a
    """
    table = parse_leap_seconds(published)
    @test length(table.delta) == 28 && table.delta[end] == 37.0
    @test table.starts == _builtin_leap_table().starts
    @test parse_leap_seconds(_valid_leap_list(4023129600)).delta == table.delta
end

@testset "leap seconds — a damaged list is refused" begin
    good = _valid_leap_list(4133980800)
    @test parse_leap_seconds(good).delta[end] == 37.0

    entries = _list_builtin_entries()
    bad = Dict(
        # Damage the hash catches. The hash is the last line, so any cut removes it.
        "cut off after 2006"  => join(split(good, '\n')[1:27], '\n'),
        "cut off mid-line"    => good[1:div(length(good), 2)],
        "an offset altered"   => replace(good, "\t37\t" => "\t99\t"),
        "a wrong hash"        => _valid_leap_list(4133980800; file_hash = "0"^40),
        "no entries"          => "#\$\t3992312697\n#@\t4133980800\n#h\t0\n",
        "not a list"          => "<html>not a leap-second list</html>\n",
        # Lists that hash correctly and are still not the list.
        "changes missing"     => _valid_leap_list(4133980800; entries = entries[[1, 28]]),
        "the last changes gone" => _valid_leap_list(4133980800; entries = entries[1:24]),
        "an offset of 5000 s" => _valid_leap_list(4133980800; extra = [(4102444800, 5000)]),
        "a jump of 2 s"       => _valid_leap_list(4133980800; extra = [(4102444800, 39)]),
        "a fractional step"   => _valid_leap_list(4133980800; extra = [(4102444800, 37.5)]),
        "a date repeated"     => _valid_leap_list(4133980800; extra = [(3692217600, 38)]),
        "an entry from 1900"  => _valid_leap_list(4133980800; entries = [(0, 5); entries]),
        "a date out of order" => _valid_leap_list(4133980800; entries = entries[[2, 1, 3:28...]]),
        "an offset not a number" => _valid_leap_list(4133980800; extra = [(4102444800, "NaN")]),
        "an expiry of Inf"    => _valid_leap_list("Inf"),
    )
    for (name, text) in bad
        @test try parse_leap_seconds(text); "accepted: $name" catch e; e isa ArgumentError ? "refused" : "$name threw $(typeof(e))" end == "refused"
    end

    # Each line the file must have
    for prefix in ("#\$", "#@", "#h")
        text = join(filter(l -> !startswith(l, prefix), split(good, '\n')), '\n')
        @test_throws ArgumentError parse_leap_seconds(text)
    end

    # A line with no offset says so, rather than failing on an index
    text = _valid_leap_list(4133980800; extra = [(4102444800, "")])
    @test_throws ArgumentError parse_leap_seconds(replace(text, "4102444800\t\t# entry" => "4102444800"))

    # A later leap second going either way is real, and passes.
    @test parse_leap_seconds(_valid_leap_list(4133980800; extra = [(4102444800, 36)])).delta[end] == 36.0
end

nothing
