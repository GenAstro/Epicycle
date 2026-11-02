# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test

@testset "Basic Time Construction" begin

    # Construction from parts
    t1 = Time(2451545.0, 0.2, TT(), JD())
    @test t1.jd == 2451545.2
    @test t1.scale == TT()
    @test t1.format == JD()
    @test t1.jd1 == 2451545.0
    @test t1.jd2 == 0.2

    # Type-based construction
    t1 = Time(2451545.0, 0.2, TT(), JD())
    @test t1.jd == 2451545.2
    @test t1.scale == TT()
    @test t1.format == JD()
    @test t1.jd1 == 2451545.0
    @test t1.jd2 == 0.2

    # Construction from parts
    t2 = Time(58000.0, TAI(), MJD())
    @test isapprox(t2.mjd, 58000.0; atol=1e-6)
    @test t2.scale == TAI()

    # Type-based construction
    t2 = Time(58000.0, TAI(), MJD())
    @test isapprox(t2.mjd, 58000.0; atol=1e-6)
    @test t2.scale == TAI()

    # Symbol based construction
    t3 = Time("2023-05-02T12:00:00", TDB(), ISOT())
    @test t3.format == ISOT()
    @test t3.scale == TDB()

    # Type based construction
    t3 = Time("2023-05-02T12:00:00", TDB(), ISOT())
    @test t3.format == ISOT()
    @test t3.scale == TDB()

    # Mixed-type constructor (Real, Real) promotes and builds Time{T}
    t_mix = Time(2451545, 0.25f0, TT(), JD())
    @test t_mix.scale === TT()
    @test t_mix.format === JD()
    @test t_mix.jd ≈ 2451545.25

end
