# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Two-number construction and calendar validation.
#
# Truth: definitions. MJD 60000 is JD 2460000.5, 2023-02-25T00:00:00, as AstroPy reads
# Time(60000.0, 0.0, format = "mjd"). A date the calendar does not have is refused with a message
# that names it.

using Test
using AstroEpochs

@testset "two numbers are parts of a date in the format named" begin
    mjd = Time(60000.0, 0.0, TT(), MJD())
    @test mjd.jd ≈ 2460000.5 atol = 1e-12
    @test mjd.mjd ≈ 60000.0 atol = 1e-12
    @test mjd.isot == "2023-02-25T00:00:00.000"
    @test mjd == Time(60000.0, TT(), MJD())                 # the single-value form agrees

    # A fraction carried in the second part keeps its precision
    split = Time(60000.0, 0.25, TT(), MJD())
    @test split.jd ≈ 2460000.75 atol = 1e-12

    jd = Time(2460000.0, 0.5, TT(), JD())
    @test jd.jd ≈ 2460000.5 atol = 1e-12
    @test jd - mjd ≈ 0.0 atol = 1e-12

    # ISOT is a string format; two numbers are refused, with the string form in the message
    err = try Time(2460000.0, 0.5, TT(), ISOT()) catch e e end
    @test err isa ArgumentError
    @test occursin("ISOT time is a string", sprint(showerror, err))
end

@testset "a date the calendar does not have is refused by name" begin
    for (text, month) in (("2024-02-30T00:00:00", "February 2024 has 29 days"),
                          ("2023-02-29T00:00:00", "February 2023 has 28 days"),
                          ("2025-04-31T00:00:00", "April 2025 has 30 days"))
        err = try Time(text, UTC(), ISOT()) catch e e end
        @test err isa ArgumentError
        @test occursin(text, sprint(showerror, err))
        @test occursin(month, sprint(showerror, err))
    end
    # The last day of each month, and 29 February in a leap year, are dates
    @test Time("2024-02-29T00:00:00", UTC(), ISOT()).isot == "2024-02-29T00:00:00.000"
    @test Time("2000-02-29T00:00:00", UTC(), ISOT()).isot == "2000-02-29T00:00:00.000"   # divisible by 400
    @test_throws ArgumentError Time("1900-02-29T00:00:00", UTC(), ISOT())                 # divisible by 100
    @test Time("2025-12-31T00:00:00", UTC(), ISOT()).isot == "2025-12-31T00:00:00.000"
end

nothing
