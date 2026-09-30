# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Edge cases found in review of the move to Astropy's algorithms, each once a wrong answer or a
# nuisance: a two-part MJD, subtraction under ForwardDiff, equality, the warnings, and the UTC day
# before 1972. Truth is Astropy where it has a value, and the definition otherwise.

using Test
using AstroEpochs
using AstroEpochs: leap_second_table
using ForwardDiff

@testset "a two-part MJD is exact, whichever part holds the fraction" begin
    # Astropy 8.0.1: Time(0.1, 0.0, format="mjd", scale="tt") has jd2 = -0.4 exactly. The zero
    # point used to be added before normalising, which put it 8 µs away.
    @test Time(0.1, 0.0, TT(), MJD()).jd2 == -0.4
    # The same decimal date, split two ways, agrees to the spacing of Float64 at 58000, 0.63 µs,
    # within which 58000.123456789 itself is stored; it used to differ by 15.6 µs. The exact
    # comparison with Astropy is the benchmark's `mjd` rows.
    a = Time(58000.0, 0.123456789, TT(), MJD())
    b = Time(58000.123456789, 0.0, TT(), MJD())
    @test abs((a - b) * 86400) < eps(58000.0) * 86400
end

@testset "subtracting days works under ForwardDiff, as adding does" begin
    t = Time(2451545.0, TT(), JD())
    @test ForwardDiff.derivative(x -> (t - x).jd, 0.5) == -1.0
    @test ForwardDiff.derivative(x -> (t + x).jd, 0.5) == 1.0
    @test (t - 0.25) == (t + (-0.25))
end

@testset "equality is the instant, to the last bit" begin
    # Summing the parts first rounded at the magnitude of a Julian date, about 40 µs.
    @test Time(2451545.0, 1e-12, TT(), JD()) != Time(2451545.0, 0.0, TT(), JD())
    # The same instant, split differently, is normalised to the same parts.
    @test Time(2451545.0, 0.25, TT(), JD()) == Time(2451544.0, 1.25, TT(), JD())
    @test Time(2451545.0, 0.25, TT(), JD()) != Time(2451545.0, 0.25, TAI(), JD())
end

@testset "UTC before 1972 warns once, however many conversions" begin
    t = Time("1965-06-15T12:00:00.000", TT(), ISOT())
    @test_logs (:warn, r"before 1972-01-01") match_mode = :all begin
        t.utc; t.utc.isot; t.utc.tai
    end
end

@testset "a Float32 time keeps its type until a scale conversion" begin
    # As the Time docstring says: the scale transforms compute in Float64.
    t = Time(2451545.0f0, 0.25f0, TT(), JD())
    @test t.jd1 isa Float32 && t.jd2 isa Float32
    @test (t + 1.0f0).jd1 isa Float32
    @test t.tdb isa Time{Float64}
    @test abs((t.tdb - Time(2451545.0, 0.25, TT(), JD()).tdb) * 86400) < 1e-6
end

@testset "the day before 1972 is an ordinary UTC day" begin
    # TAI − UTC is 0 before the table and 10 s from 1972-01-01. That step is a jump at 0h, not a
    # leap second at the end of 1971-12-31, which made that day 86410 s long and accepted 23:59:65.
    @test_throws ArgumentError Time("1971-12-31T23:59:65.000", UTC(), ISOT())
    @test_throws ArgumentError Time("1971-12-31T23:59:60.000", UTC(), ISOT())
    @test_logs (:warn,) match_mode = :any begin
        noon = Time(2441316.5, 0.5, UTC(), JD())                       # 1971-12-31T12:00 UTC
        @test noon.isot == "1971-12-31T12:00:00.000"
        @test noon.tai.isot == "1971-12-31T12:00:00.000"
        @test Time("1971-12-31T23:59:59.500", UTC(), ISOT()).tai.isot == "1971-12-31T23:59:59.500"
    end
    @test Time("1972-01-01T00:00:00.000", UTC(), ISOT()).tai.isot == "1972-01-01T00:00:10.000"
end

@testset "UTC well past the leap-second list's expiry warns once" begin
    # Only a downloaded list has an expiry; the built-in table warned when it was loaded. The
    # warning starts five years past it, as ERFA's does past its release, so a year ahead is quiet.
    if isfinite(leap_second_table().expires)
        soon = Time(2451545.0 + leap_second_table().expires + 365.0, 0.0, UTC(), JD())
        @test_logs soon.tai
        t = Time("2127-01-01T00:00:00.000", UTC(), ISOT())
        @test_logs (:warn, r"expires") match_mode = :all begin
            t.tai; t.tai.utc; t.tt
            @test (t.tai - Time("2127-01-01T00:00:00.000", TAI(), ISOT())) * 86400 ≈ 37.0
        end
    end
end
