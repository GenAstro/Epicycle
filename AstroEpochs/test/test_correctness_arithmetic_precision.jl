# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Below a nanosecond: arithmetic, equality and the normalisation of the two-part date. Holding these
# separately from the conversions matters, because the reason to split a date is that a difference
# far smaller than a single float at the size of a Julian date survives.
#
# Truth, by section:
#   differences and sums    the definition: the two parts are added by exact sums, so the offset
#                           that went in comes back, to the spacing of a float near 0.5 day
#                           (about 5 ps).
#   equality                the definition: a different instant is a different time.
#   normalisation           Astropy 8.0.2, astropy.time.utils.day_frac, at the exact ties.
#   observer terms          pyerfa 2.0.1.5 dtdb, written by astropy/make_dtdb_observer.py.

using Test
using AstroEpochs
using AstroEpochs: _rebalance, dtdb

const _PS = 1e-12 / 86400.0             # one picosecond, in days

@testset "precision — an offset below a nanosecond comes back from a difference" begin
    for jd in (2378497.0, 2451545.0, 2460000.0, 2524593.0), frac in (-0.4, -0.1, 0.0, 0.25, 0.49)
        a = Time(jd, frac, TT(), JD())
        for offset_ps in (5_000, 86, 50, 10)           # 5 ns down to 10 ps
            offset = offset_ps * _PS
            b = Time(jd, frac + offset, TT(), JD())
            @test abs((b - a) - offset) < 5 * _PS
            @test abs((a - b) + offset) < 5 * _PS
        end
    end
end

@testset "precision — adding an interval keeps it to picoseconds" begin
    for jd in (2378497.0, 2451545.0, 2460000.0), frac in (-0.4, 0.0, 0.25)
        a = Time(jd, frac, TT(), JD())
        for offset_ps in (5_000, 1_000, 86, 50)
            offset = offset_ps * _PS
            @test abs(((a + offset) - a) - offset) < 5 * _PS
            @test abs(((a - offset) - a) + offset) < 5 * _PS
            @test abs(((offset + a) - a) - offset) < 5 * _PS
        end
        # Whole and part of a day together: the part is not lost to the size of the whole.
        @test abs(((a + (10.0 + 100 * _PS)) - a) - (10.0 + 100 * _PS)) < 1e-11
    end
end

@testset "precision — many small steps do not accumulate more than the rounding allows" begin
    t0 = Time(2451545.0, 0.0, TT(), JD())
    t  = t0
    for _ in 1:86_400
        t += 1 / 86_400
    end
    @test abs((t - t0) - 1.0) * 86_400 < 1e-6
end

@testset "precision — equality sees a difference below a nanosecond" begin
    a = Time(2451545.0, 0.25, TT(), JD())
    @test a == Time(2451545.0, 0.25, TT(), JD())
    @test a != Time(2451545.0, 0.25 + 1e-15, TT(), JD())      # 86 ps
    @test a != Time(2451545.0, 0.25 - 1e-15, TT(), JD())
    @test a != a + 100 * _PS
    @test Time(2451545.25, 0.0, TT(), JD()) == Time(2451545.0, 0.25, TT(), JD())   # one instant, split differently
end

@testset "precision — the date is split as Astropy splits it, including at exact ties" begin
    cases = [((2451545.0, 0.5),            (2451546.0, -0.5)),
             ((2451544.0, 0.5),            (2451544.0, 0.5)),
             ((2451545.0, -0.5),           (2451544.0, 0.5)),
             ((2451544.0, -0.5),           (2451544.0, -0.5)),
             ((2451545.0, 0.5000000000000002), (2451546.0, -0.4999999999999998)),
             ((2451545.0, 0.4999999999999999), (2451545.0, 0.4999999999999999)),
             ((2451545.5, 0.0),            (2451546.0, -0.5)),
             ((2451544.5, 0.0),            (2451544.0, 0.5)),
             ((2451545.25, 0.25),          (2451546.0, -0.5)),
             ((1e-20, 2451544.5),          (2451544.0, 0.5)),
             # A half-day where the rounding error decides which way the day goes
             ((2451544.5, 5e-17),          (2451545.0, -0.49999999999999994)),
             ((2451544.5, -5e-17),         (2451544.0, 0.49999999999999994)),
             ((2451545.5, 5e-17),          (2451546.0, -0.49999999999999994)),
             ((2451545.5, -5e-17),         (2451545.0, 0.49999999999999994)),
             ((2451544.5, 1e-20),          (2451544.0, 0.5)),
             ((2451545.5, -1e-20),         (2451546.0, -0.5))]
    for ((a, b), (day, frac)) in cases
        @test _rebalance(a, b) == (day, frac)
    end
end

@testset "precision — TDB − TT with observer terms agrees with pyerfa" begin
    lines = readlines(joinpath(@__DIR__, "astropy", "dtdb_observer.csv"))[2:end]
    @test length(lines) == 100
    for line in lines
        date1, date2, ut, elong, u, v, truth = parse.(Float64, split(line, ','))
        @test dtdb(date1, date2, ut, elong, u, v) ≈ truth atol = 1e-17 rtol = 0
    end
end

nothing
