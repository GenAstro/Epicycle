using Test
using AstroEpochs

@testset "Time Arithmetic Tests" begin
    # Construct a base time
    t0 = Time("2025-07-31T12:00:00", :tt, :isot)
    
    # Scalar addition
    t1 = t0 + 0.25  # add 6 hours
    @test isapprox(t1.jd - t0.jd, 0.25; atol=1e-12)

    # Scalar subtraction
    t2 = t1 - 0.5  # subtract 12 hours
    @test isapprox(t2.jd - t0.jd, -0.25; atol=1e-12)

    # Commutative addition
    t3 = 1.0 + t0
    @test isapprox(t3.jd - t0.jd, 1.0; atol=1e-12)

    # Difference between times
    dt = t1 - t0
    @test isapprox(dt, 0.25; atol=1e-12)

    # Rebalancing test near day boundaries
    t4 = Time(t0.jd1, 0.499999999, :tt, :jd)
    t5 = t4 + 1e-6
    @test abs(t5.jd2) < 0.5  # Should rebalance so jd2 is still within [-0.5, 0.5]

    # Time scale mismatch should error on subtraction
    t6 = Time(t0.jd1, t0.jd2, :tdb, :jd)
    @test_throws ErrorException t6 - t0

    # Time format mismatch should error on subtraction
    t7 = Time(t0.jd1, t0.jd2, :tt, :mjd)
    @test_throws ErrorException t7 - t0
end

# Unit test for _rebalance function.  while jd2 <= half
@test begin
    jd1, jd2 = AstroEpochs._rebalance(0.0, -1.2)
    jd1 == -1.0 && -0.5 <= jd2 < 0.5
end

# Unit test for _rebalance function.  while jd2 >= half
@test begin
    jd1, jd2 = AstroEpochs._rebalance(0.0, 1.2)
    jd1 == 1.0 && 0.5 >= jd2 < 0.5
end

@test begin
    # Heterogeneous _rebalance fallback (promote then rebalance)
    jd1b, jd2b = AstroEpochs._rebalance(1.0, big(0.75))
    (jd1b == big(2.0)) && (jd2b == big(-0.25))
end

@testset "Time == and isequal tests" begin
    # Equality across numeric types and split variants
    a = Time(2451545.0, 0.25, :tt, :jd)
    b = Time(big(2451545.0), big(0.25), :tt, :jd)  # different T, same value
    @test a == b

    # Same value but different split (still equal by jd sum)
    c = Time(a.jd1 + 0.75, a.jd2 - 0.75, :tt, :jd)
    @test a == c

    # Inequality due to different scale and different format (early returns)
    d = Time(a.jd1, a.jd2, :tai, :jd)
    e = Time(a.mjd, :tt, :mjd)
    @test a != d
    @test a != e
end