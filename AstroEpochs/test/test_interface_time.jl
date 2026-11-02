# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The public interface other packages use: tags for scale and format, changing a held time's scale,
# format and numeric type, calendar fields, and the offsets between scales.
#
# Truth: the definitions, and the leap-second list for the leap-second cases.

using Test
using AstroEpochs
using ForwardDiff

@testset "interface — scale and format are tags" begin
    t = Time(2451545.0, TT(), JD())
    @test t.scale === TT()
    @test t.format === JD()
    @test t.scale == TT() && t.scale != TDB()
    @test t.tdb.scale === TDB()
    @test Time(t, ISOT()).format === ISOT()
    for (s, name) in ((TT(), :tt), (TAI(), :tai), (TDB(), :tdb), (UTC(), :utc), (TCB(), :tcb), (TCG(), :tcg))
        @test Time(2451545.0, s, MJD()).scale === s
        @test getproperty(t, name).scale === s
    end
end

@testset "interface — symbols are not a way in" begin
    @test_throws ArgumentError Time(2451545.0, :tt, :jd)
    @test_throws ArgumentError Time(2451545.0, 0.0, :tt, :jd)
    @test_throws ArgumentError Time("2024-01-01T00:00:00", :tt, :isot)
end

@testset "interface — integers are promoted, so a whole-number date works" begin
    @test Time(60000, UTC(), MJD()) isa Time{Float64}
    @test Time(60000, UTC(), MJD()) == Time(60000.0, UTC(), MJD())
    @test Time(60000, UTC(), MJD()).isot == "2023-02-25T00:00:00.000"
    @test Time(2451545, TT(), JD()) isa Time{Float64}
    @test Time(2451545, TT(), JD()) == Time(2451545.0, TT(), JD())
    @test Time(60000, 0, UTC(), MJD()) == Time(60000.0, UTC(), MJD())
    @test Time(2451545, 0, TT(), JD()) isa Time{Float64}
    @test (Time(2451545, 0, TT(), JD()) + 0.5).jd == 2451545.5
    @test Time(Int32(60000), TT(), MJD()) isa Time{Float64}
end

@testset "interface — Time(t, scale) converts, for a scale held in a variable" begin
    a = Time(Time(2451545.0, 0.1, TT(), JD()), ISOT())
    @test Time(a, TDB()) == a.tdb
    @test Time(a, TDB()).format === ISOT()
    @test Time(a, TT()) == a

    b = Time(2451545.5, 0.1, TDB(), JD())
    @test Time(b, a.scale) - a == b.tt - a
end

@testset "interface — Time(t, format) keeps the date to the last bit" begin
    t = Time(2460000.0, 0.1234567890123, TT(), JD())
    for f in (JD(), MJD(), ISOT())
        r = Time(t, f)
        @test r.format === f
        @test r.scale === TT()
        @test (r.jd1, r.jd2) === (t.jd1, t.jd2)
        @test r - t == 0.0
    end
end

@testset "interface — Time{T}(t) changes the numeric type and nothing else" begin
    t = Time(2451545.0, 0.25, TDB(), MJD())
    b = Time{BigFloat}(t)
    @test b isa Time{BigFloat}
    @test (b.scale, b.format) === (TDB(), MJD())
    @test b.jd1 == t.jd1 && b.jd2 == t.jd2

    @test convert(Time{BigFloat}, t) isa Time{BigFloat}
    @test convert(Time{Float64}, t) === t

    D = ForwardDiff.Dual{Nothing, Float64, 1}
    d = Time{D}(t)
    @test d isa Time{D}
    @test ForwardDiff.value(d.jd2) == t.jd2
end

@testset "interface — isot(t; digits) writes as many decimals as asked" begin
    t = Time(2451545.0, 0.25, TT(), JD())
    @test isot(t) == t.isot == "2000-01-01T18:00:00.000"
    @test isot(t; digits = 0) == "2000-01-01T18:00:00"
    @test isot(t; digits = 6) == "2000-01-01T18:00:00.000000"

    # Rounding carries, and a leap second reads 60, at any number of decimals.
    @test isot(Time("2020-03-01T23:59:59.9999", TT(), ISOT()); digits = 3) == "2020-03-02T00:00:00.000"
    @test isot(Time("2020-03-01T23:59:59.9999", TT(), ISOT()); digits = 6) == "2020-03-01T23:59:59.999900"
    @test isot(Time("2016-12-31T23:59:60.25", UTC(), ISOT()); digits = 2) == "2016-12-31T23:59:60.25"
end

@testset "interface — a time from calendar fields is the time from its string" begin
    a = Time(2024, 3, 1, 12, 30, 15.25, TAI(), JD())
    b = Time("2024-03-01T12:30:15.25", TAI(), ISOT())
    @test (a.jd1, a.jd2) == (b.jd1, b.jd2)
    @test a.format === JD()
    @test Time(2016, 12, 31, 23, 59, 60.5, UTC(), ISOT()).isot == "2016-12-31T23:59:60.500"

    # A second of 60 is real only on a day that ends with a leap second, and only in UTC.
    @test_throws ArgumentError Time(2024, 6, 1, 23, 59, 60.0, UTC(), ISOT())
    @test_throws ArgumentError Time(2016, 12, 31, 23, 59, 60.5, TT(), ISOT())
end

@testset "interface — times have an order, in one scale" begin
    a = Time(2451545.0, 0.25, TT(), JD())
    b = a + 1e-15                                        # 86 ps later
    c = Time(a, MJD()) + 1.0
    @test a < b && !(b < a) && !(a < a)
    @test a <= a && b > a && b >= a
    @test a < c                                          # across formats
    @test sort([c, b, a]) == [a, b, c]
    @test minimum([c, b, a]) == a && maximum([c, b, a]) == c
    @test min(a, c) == a && max(a, c) == c
    @test issorted([a, b, c])

    err = try a < a.tdb; nothing catch e; e end
    @test err isa ArgumentError
    @test occursin("cannot compare times in different scales", sprint(showerror, err))
    @test occursin("convert one first", sprint(showerror, err))
    @test_throws ArgumentError sort([a, a.tdb])
    @test_throws ArgumentError min(a, a.tdb)

    # A half-day either way of a split still orders by the instant
    @test Time(2451544.0, 0.5, TT(), JD()) < Time(2451545.0, -0.25, TT(), JD())
    @test Time(2451546.0, -0.5, TT(), JD()) > Time(2451545.0, 0.25, TT(), JD())
end

@testset "interface — isapprox compares in seconds, in one scale" begin
    a = Time(2460000.123, TT(), JD())
    b = Time(59999.623, TT(), MJD())
    @test a != b
    @test !(a ≈ b)                                       # 12 µs apart, the default is 1 ns
    @test isapprox(a, b; atol = 1e-3)
    @test !isapprox(a, b; atol = 1e-6)
    @test a ≈ a
    @test a ≈ a + 1e-10 / 86400                          # 0.1 ns
    @test !(a ≈ a + 1e-8 / 86400)                        # 10 ns
    @test_throws ArgumentError a ≈ a.tdb
end

@testset "interface — a time is one thing in a broadcast" begin
    t = Time(2451545.0, 0.25, TT(), JD())
    @test t .+ [0.0, 1.0] == [t, t + 1.0]
    @test [0.0, 1.0] .+ t == [t, t + 1.0]
    @test t .+ (0:0.5:1) == [t, t + 0.5, t + 1.0]
    @test t .- [0.5] == [t - 0.5]
    @test (t + 1.0) .- [t, t + 1.0] == [1.0, 0.0]
    @test t .+ 1.0 == t + 1.0
end

@testset "interface — any string type makes a time" begin
    text = split("2024-01-01T00:00:00 trailing")[1]
    @test text isa SubString
    @test Time(text, TT(), ISOT()) == Time("2024-01-01T00:00:00", TT(), ISOT())
    @test_throws ArgumentError Time(text, TT(), JD())
end

@testset "interface — propertynames lists what works" begin
    t = Time(2451545.0, TT(), JD())
    names = propertynames(t)
    for name in names
        @test getproperty(t, name) !== nothing
        @test hasproperty(t, name)
    end
    @test all(n -> n in names, (:jd, :mjd, :isot, :tt, :tai, :tdb, :utc, :tcg, :tcb, :scale, :format))
    @test !hasproperty(t, :_jd1) && !hasproperty(t, :_scale)
    @test :_jd1 in propertynames(t, true)
end

@testset "interface — no property loses its type" begin
    t = Time(2451545.0, 0.25, TT(), JD())
    @test (@inferred (x -> x.jd)(t)) isa Float64
    @test (@inferred (x -> x.mjd)(t)) isa Float64
    @test (@inferred (x -> x.jd1)(t)) isa Float64
    @test (@inferred (x -> x.isot)(t)) isa String
    @test (@inferred (x -> x.tdb)(t)) isa Time{Float64}
    @test (@inferred (x -> x.utc)(t)) isa Time{Float64}
    D = ForwardDiff.Dual{Nothing, Float64, 1}
    @test (@inferred (x -> x.jd)(Time{D}(t))) isa D
end

@testset "interface — a time with no finite date passes through and prints" begin
    for bad in (Time(NaN, TT(), JD()), Time(NaN, UTC(), JD()), Time(Inf, TT(), JD()),
                Time(2451545.0, TT(), JD()) + NaN, Time(Time(2451545.0, UTC(), JD()), ISOT()) + Inf)
        for name in (:tt, :tai, :tdb, :utc, :tcg, :tcb)
            @test isnan(getproperty(bad, name).jd)
        end
        @test isnan(bad.jd)
        @test isnan((bad + 1.0).jd)
        @test occursin("value", sprint(show, bad))
        err = try bad.isot; nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("no calendar date", sprint(showerror, err))
    end
    @test occursin("no calendar date", sprint(show, Time(Time(NaN, TT(), JD()), ISOT())))

    # So does a date too far off for the calendar, and printing still works.
    huge = Time(Time(2451545.0, TT(), JD()) + 1e10, ISOT())
    @test occursin("no calendar date", sprint(show, huge))

    # A derivative survives to the end of the computation, as NaN.
    s = ForwardDiff.Dual{Nothing}(NaN, 1.0)
    @test isnan(ForwardDiff.value((Time(2451545.0, UTC(), JD()) + s).tdb.jd))
end

@testset "interface — TAI − UTC at a time, whatever its scale" begin
    @test tai_minus_utc(Time("2024-01-01T00:00:00", UTC(), ISOT())) == 37.0
    @test tai_minus_utc(Time("2024-01-01T00:00:00", TT(), ISOT())) == 37.0
    @test tai_minus_utc(Time("2016-12-31T23:59:59", UTC(), ISOT())) == 36.0
    @test tai_minus_utc(Time("2017-01-01T00:00:00", UTC(), ISOT())) == 37.0
end

@testset "interface — TDB − TT at a time, at the geocentre and at an observer" begin
    t = Time("2024-01-01T00:00:00", TT(), ISOT())

    # The geocentric value is what the scale conversion applies.
    @test (t.tdb - Time(t.jd1, t.jd2, TDB(), JD())) * 86400 ≈ tdb_minus_tt(t) atol = 1e-10

    # An observer changes it by microseconds at most, and the geocentre is where the terms vanish.
    r = [-2353.6, -4641.3, 3677.1]
    @test 0 < abs(tdb_minus_tt(t, r) - tdb_minus_tt(t)) < 3e-6
    @test tdb_minus_tt(t, [0.0, 0.0, 0.0]) ≈ tdb_minus_tt(t) atol = 1e-15

    # The same instant in another scale gives the same answer.
    @test tdb_minus_tt(t.tai, r) ≈ tdb_minus_tt(t, r) atol = 1e-11
end

nothing
