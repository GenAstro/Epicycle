# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using ForwardDiff

# Zygote is checked when the environment has it, as `Pkg.test` does, and not in the shared
# development environment, which does not carry it.
const HAS_ZYGOTE = Base.find_package("Zygote") !== nothing
HAS_ZYGOTE && @eval using Zygote

@testset "Dual Time construction" begin
    @testset "TAI + MJD from Dual value" begin
        x = ForwardDiff.Dual{Nothing}(60369.52426068287, 1.0)
        t = Time(x, TAI(), MJD())
        @test t.scale == TAI()
        @test t.format == MJD()
        @test t.jd1 isa typeof(x)
        @test t.jd2 isa typeof(x)
        @test t.mjd == x
        # Arithmetic preserves Dual
        t2 = t + 1.0
        @test t2.jd1 isa typeof(x)
        @test t2.mjd == x + one(x)
    end

    @testset "UTC + JD from Dual value" begin
        y = ForwardDiff.Dual{Nothing}(2.452090024260688e6, 1.0)
        t = Time(y, UTC(), JD())
        @test t.scale == UTC()
        @test t.format == JD()
        @test t.jd1 isa typeof(y)
        @test t.jd2 isa typeof(y)
        @test t.jd == y
        # Subtraction preserves Dual
        t2 = t - 2.0
        @test t2.jd1 isa typeof(y)
        @test t2.jd == y - ForwardDiff.Dual{Nothing}(2.0, 0.0)
    end

    @testset "TDB + precision JD from Dual jd1/jd2" begin
        jd1d = ForwardDiff.Dual{Nothing}(2.45209e6, 1.0)
        jd2d = ForwardDiff.Dual{Nothing}(0.024260688107460737, 0.0)
        t = Time(jd1d, jd2d, TDB(), JD())
        @test t.scale == TDB()
        @test t.format == JD()
        # The date is normalised as Astropy's is (day_frac): a whole day and a fraction in
        # [-0.5, 0.5]. The whole day is a rounded number and carries no derivative, so the
        # derivative moves to the fraction, and the date as a whole keeps value and derivative.
        @test t.jd1 isa typeof(jd1d) && t.jd2 isa typeof(jd2d)
        @test ForwardDiff.value(t.jd1) == round(ForwardDiff.value(jd1d + jd2d))
        @test ForwardDiff.partials(t.jd1)[1] == 0.0
        @test ForwardDiff.partials(t.jd2)[1] == 1.0
        @test t.jd === jd1d + jd2d
        # Round-trip through no-op add
        t2 = t + 0.0
        @test t2.jd1 isa typeof(jd1d)
        @test t2.jd2 isa typeof(jd2d)
        @test t2.jd === t.jd
    end
end

@testset "Test ForwardDiff differentiation" begin
    # Base time in TAI
    jd1 = 2.45209e6
    jd2 = 0.024260688107460737
    t0 = Time(jd1, jd2, TAI(), JD())

    # Helper: add seconds in a given scale without touching production code
    # - Convert t0 to `scale`, add sec as days, convert back to t0.scale, return jd
    jd_after_add_in_scale(t::Time, sec, scale::Symbol) = begin
        t_in = getproperty(t, scale)
        dt_days = sec / (one(sec) * 86400.0)  # keep AD type by using one(sec)
        t_back = Time(t_in + dt_days, t.scale)
        t_back.jd
    end

    # 1) Type-stability with Dual seconds input (only dt is Dual)
    s = ForwardDiff.Dual{Nothing}(10.0, 1.0)
    let
        t_in = getproperty(t0, :tdb)
        dt_days = s / (one(s) * 86400.0)
        t_back = Time(t_in + dt_days, t0.scale)
        @test t_back.scale == TAI()
        @test t_back.format == JD()
        @test t_back.jd1 isa typeof(s)
        @test t_back.jd2 isa typeof(s)
    end

    # 2) Differentiate jd (numeric JD) w.r.t. seconds added in TDB
    g(sec) = jd_after_add_in_scale(t0, sec, :tdb)
    d_fwd = ForwardDiff.derivative(g, 0.0)
    @test d_fwd > 0
    @test isapprox(d_fwd, 1/86400; rtol=1e-6)

    # Central difference in seconds, on the two-part date so it is not limited by one float at the
    # size of a Julian date. TDB and TAI differ in rate by about 1e-8, so this tells them apart.
    d_fd = (Time(t0.tdb + 1 / 86400, TAI()) - Time(t0.tdb - 1 / 86400, TAI())) / 2
    @test isapprox(d_fwd, d_fd; rtol=1e-9)

    # Zygote gradient w.r.t. seconds
    if HAS_ZYGOTE
        d_zyg, = Zygote.gradient(g, 0.0)
        @test d_zyg > 0
        @test isapprox(d_zyg, 1/86400; rtol=1e-6)
        @test isapprox(d_zyg, d_fwd; rtol=1e-9)
    end

    # Every scale pair: seconds on one scale against seconds on another, against a central
    # difference on the two-part date. The rates differ by parts in 1e8, so this tells them apart.
    base = Time("2024-06-01T06:00:00", TT(), ISOT())
    scales = (:tt, :tai, :tdb, :utc, :tcg, :tcb)
    for a in scales, b in scales
        a === b && continue
        elapsed(s) = begin
            from = getproperty(base, a)
            (getproperty(from + s / 86400, b) - getproperty(from, b)) * 86400
        end
        @test ForwardDiff.derivative(elapsed, 0.0) ≈ (elapsed(1.0) - elapsed(-1.0)) / 2 atol = 1e-10
    end

    # 3) Round-trip sanity: adding 0 seconds leaves epoch unchanged (value identity)
    t_same = Time(t0.tdb + (0.0 / 86400.0), t0.scale)
    @test t_same.scale == t0.scale
    @test t_same.format == t0.format
    @test t_same.jd == t0.jd
end
nothing
