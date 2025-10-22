using Test

# TODO. When refactoring is done, modify the symbol tests below as needed.
@testset "Basic Time Construction" begin

    # Symbol-based construction
    t1 = Time(2451545.0, 0.2, :tt, :jd)
    @test t1.jd == 2451545.2
    @test t1.scale == :tt
    @test t1.format == :jd
    @test t1.jd1 == 2451545.0
    @test t1.jd2 == 0.2

    # Type-based construction
    t1 = Time(2451545.0, 0.2, TT(), JD())
    @test t1.jd == 2451545.2
    @test t1.scale == :tt
    @test t1.format == :jd
    @test t1.jd1 == 2451545.0
    @test t1.jd2 == 0.2

    # Symbol-based construction
    t2 = Time(58000.0, :tai, :mjd)
    @test isapprox(t2.mjd, 58000.0; atol=1e-6)
    @test t2.scale == :tai

    # Type-based construction
    t2 = Time(58000.0, TAI(), MJD())
    @test isapprox(t2.mjd, 58000.0; atol=1e-6)
    @test t2.scale == :tai

    # Symbol based construction
    t3 = Time("2023-05-02T12:00:00", :tdb, :isot)
    @test t3.format == :isot
    @test t3.scale == :tdb

    # Type based construction
    t3 = Time("2023-05-02T12:00:00", TDB(), ISOT())
    @test t3.format == :isot
    @test t3.scale == :tdb

    # Mixed-type constructor (Real, Real) promotes and builds Time{T}
    t_mix = Time(2451545, 0.25f0, :tt, :jd)
    @test t_mix.scale === :tt
    @test t_mix.format === :jd
    @test t_mix.jd ≈ 2451545.25

end
