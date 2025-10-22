"""
    isapproxrel(a::Float64, b::Float64; rtol=1e-12, atol=1e-12) -> Bool

Returns `true` if `a` and `b` are approximately equal using relative tolerance,
or absolute tolerance if both are small (< 1.0).
"""
function isapproxrel(a::Float64, b::Float64; rtol::Float64=1e-14, atol::Float64=1e-14)::Bool
    if abs(a) < 1.0 && abs(b) < 1.0
        return abs(a - b) ≤ atol
    else
        return abs(a - b) ≤ rtol * max(abs(a), abs(b))
    end
end

@testset "Leap Day Format Conversion - 2024-02-29T12:34:56.123 TAI" begin
    # Reference values from AstroPy
    ref_jd1 = 2460370.0
    ref_jd2 = 0.024259259259259314
    ref_jd  = 2460370.0242592595
    ref_mjd = 60369.52425925926
    ref_iso = "2024-02-29T12:34:56.000"

    # Symbol based construction
    t = Time("2024-02-29T12:34:56", :tai, :isot)
    @test isapprox(t.jd1, ref_jd1; atol=1e-12)
    @test isapprox(t.jd2, ref_jd2; atol=1e-15)
    @test isapprox(t.jd,  ref_jd;  atol=1e-12)
    @test isapprox(t.mjd, ref_mjd; atol=1e-12)
    @test t.isot == ref_iso

    # Type based construction
    t = Time("2024-02-29T12:34:56", TAI(), ISOT())
    @test isapprox(t.jd1, ref_jd1; atol=1e-12)
    @test isapprox(t.jd2, ref_jd2; atol=1e-15)
    @test isapprox(t.jd,  ref_jd;  atol=1e-12)
    @test isapprox(t.mjd, ref_mjd; atol=1e-12)
    @test t.isot == ref_iso
end

@testset "ISO Format Conversion — nanosecond Precision" begin

    # Reference values from AstroPy
    ref_jd1 = 2452090.0
    ref_jd2 = 0.024260688157280108
    ref_jd  = 2452090.024260688
    ref_mjd = 52089.52426068816
    ref_isot = "2001-06-29T12:34:56.123"

    # Symbol based construction from high-precision ISO string
    t = Time("2001-06-29T12:34:56.123456789", :tai, :isot)
    @test isapprox(t.jd1, ref_jd1; atol=1e-12)
    @test isapprox(t.jd2, ref_jd2; atol=1e-15)
    @test isapprox(t.jd,  ref_jd;  atol=1e-12)
    @test isapprox(t.mjd, ref_mjd; atol=1e-12)
    @test t.isot == ref_isot  # rounded to ms

    # Type based construction from high-precision ISO string
    t = Time("2001-06-29T12:34:56.123456789", TAI(), ISOT())
    @test isapprox(t.jd1, ref_jd1; atol=1e-12)
    @test isapprox(t.jd2, ref_jd2; atol=1e-15)
    @test isapprox(t.jd,  ref_jd;  atol=1e-12)
    @test isapprox(t.mjd, ref_mjd; atol=1e-12)
    @test t.isot == ref_isot  # rounded to ms
end

@testset "Time Format I/O Permutations" begin
    # Truth values from Astropy
    jd_ref  = 2450628.5712382756
    jd1_ref = 2450629.0
    jd2_ref = -0.42876172453703704
    mjd_ref = 50628.071238275465
    isot_ref = "1997-06-29T01:42:34.987"

    # Symbol based ISOT construction
    t1 = Time("1997-06-29T01:42:34.987", :tai, :isot)
    @test isapproxrel(t1.jd, jd_ref)
    @test isapproxrel(t1.jd1, jd1_ref)
    @test isapproxrel(t1.jd2, jd2_ref)
    @test isapproxrel(t1.mjd, mjd_ref)
    @test t1.isot == isot_ref

    # Type based ISOT construction
    t1 = Time("1997-06-29T01:42:34.987", TAI(), ISOT())
    @test isapproxrel(t1.jd, jd_ref)
    @test isapproxrel(t1.jd1, jd1_ref)
    @test isapproxrel(t1.jd2, jd2_ref)
    @test isapproxrel(t1.mjd, mjd_ref)
    @test t1.isot == isot_ref

    # Symbol based construction from jd
    t2 = Time(jd_ref, :tai, :jd)
    @test isapproxrel(t2.jd, jd_ref)
    @test isapproxrel(t2.jd1, jd1_ref)
    @test isapproxrel(t2.jd2, jd2_ref; atol = 1e-9, rtol = 1e-9) # Setting time using a JD limits precision in jd2
    @test isapproxrel(t2.mjd, mjd_ref)
    @test t2.isot == isot_ref

    # Type based construction from jd
    t2 = Time(jd_ref, TAI(), JD())
    @test isapproxrel(t2.jd, jd_ref)
    @test isapproxrel(t2.jd1, jd1_ref)
    @test isapproxrel(t2.jd2, jd2_ref; atol = 1e-9, rtol = 1e-9) # Setting time using a JD limits precision in jd2
    @test isapproxrel(t2.mjd, mjd_ref)
    @test t2.isot == isot_ref

    # Symbol based construction from mjd
    t3 = Time(mjd_ref, :tai, :mjd)
    @test isapproxrel(t3.jd, jd_ref)
    @test isapproxrel(t3.jd1, jd1_ref - 0.5)
    @test isapproxrel(t3.jd2, jd2_ref + 0.5; atol = 1e-11, rtol = 1e-11) # Setting a time in MJD limits precision in jd2
    @test isapproxrel(t3.mjd, mjd_ref)
    @test t3.isot == isot_ref

    # Type based construction from mjd
    t3 = Time(mjd_ref, TAI(), MJD())
    @test isapproxrel(t3.jd, jd_ref)
    @test isapproxrel(t3.jd1, jd1_ref - 0.5)
    @test isapproxrel(t3.jd2, jd2_ref + 0.5; atol = 1e-11, rtol = 1e-11) # Setting a time in MJD limits precision in jd2
    @test isapproxrel(t3.mjd, mjd_ref)
    @test t3.isot == isot_ref
end
