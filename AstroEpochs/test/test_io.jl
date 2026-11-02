# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

#=
Note, testing of valid properties is implicitly tested in test_formatconversions
and test_scaleconversions.  These only test other areas not covered in those tests. 
=#

@testset "getproperty() fallback behavior" begin
    t = Time("2020-01-01T00:00:00", TT(), ISOT())

    # Accessing valid field works
    @test typeof(t.jd1) == Float64

    # Accessing invalid field throws expected error
    @test_throws ErrorException t.foobar
end

@testset "Test Time show method" begin
    t1 = Time("2024-02-29T12:34:56.123", TAI(), ISOT())
    t2 = Time(60369.52426068287, TAI(), MJD())
    t3 = Time(2452090.024260688, TAI(), JD())

    output1 = sprint(show, t1)
    output2 = sprint(show, t2)
    output3 = sprint(show, t3)

    @test occursin("AstroEpochs.Time", output1)
    @test occursin("scale  = TAI()", output1)
    @test occursin("format = ISOT()", output1)
    @test occursin("2024-02-29T12:34:56.123", output1)

    @test occursin("60369.52426068287", output2)
    @test occursin("format = MJD()", output2)

    @test occursin("2.452090024260688e6", output3)
    @test occursin("format = JD()", output3)
end

@testset "User Input Validation Tests" begin

    #Invalid constructors
    @test_throws ArgumentError Time(2451545.0, 0.0, :junk, JD())
    @test_throws ArgumentError Time(2451545.0, 0.0, TT(), :banana)
    @test_throws ArgumentError Time(24515.0, TT(), ISOT())
    @test_throws ArgumentError Time("1980-01-06T00:00:00", TT(), MJD())
    @test_throws ArgumentError Time(2456666.0, "b", TT(), MJD())
    @test_throws ArgumentError Time("apple", 1.0 , TT(), MJD())
    @test_throws ArgumentError Time("apple", "pear" , TT(), MJD())
    @test_throws ArgumentError Time("20253-05-31T12:34:56.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-13-31T12:34:56.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-06-51T12:34:56.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-06-23T26:34:56.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-11-30T12:62:56.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-07-11T12:34:61.123", TAI(), ISOT())
    @test_throws ArgumentError Time("2025-07-31T12:34:056.123", TAI(), ISOT())

end
