using Test
using ForwardDiff
using FiniteDiff
using Zygote

@testset "Dual Time construction" begin
    @testset "TAI + MJD from Dual value" begin
        x = ForwardDiff.Dual{Nothing}(60369.52426068287, 1.0)
        t = Time(x, TAI(), MJD())
        @test t.scale == :tai
        @test t.format == :mjd
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
        @test t.scale == :utc
        @test t.format == :jd
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
        @test t.scale == :tdb
        @test t.format == :jd
        @test t.jd1 === jd1d
        @test t.jd2 === jd2d
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
        t_back = getproperty(t_in + dt_days, t.scale)
        t_back.jd
    end

    # 1) Type-stability with Dual seconds input (only dt is Dual)
    s = ForwardDiff.Dual{Nothing}(10.0, 1.0)
    let
        t_in = getproperty(t0, :tdb)
        dt_days = s / (one(s) * 86400.0)
        t_back = getproperty(t_in + dt_days, t0.scale)
        @test t_back.scale == :tai
        @test t_back.format == :jd
        @test t_back.jd1 isa typeof(s)
        @test t_back.jd2 isa typeof(s)
    end

    # 2) Differentiate jd (numeric JD) w.r.t. seconds added in TDB
    g(sec) = jd_after_add_in_scale(t0, sec, :tdb)
    d_fwd = ForwardDiff.derivative(g, 0.0)
    @test d_fwd > 0
    @test isapprox(d_fwd, 1/86400; rtol=1e-6)

    # FiniteDiff cross-check (central difference in seconds)
    #d_fd = FiniteDiff.finite_difference_derivative(g, 0.0, Val(:central), Float64, nothing; absstep=1e-3)
    #@test isapprox(d_fwd, d_fd; rtol=1e-9)

    # Zygote gradient w.r.t. seconds
    d_zyg, = Zygote.gradient(g, 0.0)
    @test d_zyg > 0
    @test isapprox(d_zyg, 1/86400; rtol=1e-6)
    @test isapprox(d_zyg, d_fwd; rtol=1e-9)

    # 3) Round-trip sanity: adding 0 seconds leaves epoch unchanged (value identity)
    t_same = getproperty(getproperty(t0, :tdb) + (0.0 / 86400.0), t0.scale)
    @test t_same.scale == t0.scale
    @test t_same.format == t0.format
    @test t_same.jd == t0.jd
end
#=
@testset "Test ForwardDiff differentiation" begin
    # Base time in TAI
    jd1 = 2.45209e6
    jd2 = 0.024260688107460737
    t0 = Time(jd1, jd2, TAI(), JD())

    # 1) Type-stability with Dual seconds input
    s = ForwardDiff.Dual{Nothing}(10.0, 1.0)
    t1 = add_seconds_in_scale(t0, s, :tdb)
    @test t1.scale == :tai
    @test t1.format == :jd
    @test t1.jd1 isa typeof(s)
    @test t1.jd2 isa typeof(s)

    # 2) Differentiate jd (numeric JD) w.r.t. seconds added in TDB
    #    Expect derivative ~ 1/86400 with small deviation due to scale offsets.
    dt = 86400.0  # 1 day in seconds
    g(sec) = add_seconds_in_scale(t0, sec, :tdb).jd
    d_fwd   = ForwardDiff.derivative(g, dt)
    #println("ForwardDiff derivative (jd w.r.t. tdb seconds): $d_fwd")
    @test d_fwd > 0
    @test isapprox(d_fwd, 1/86400; rtol=1e-6)

    # Zygote gradient w.r.t. seconds
    d_zyg, = Zygote.gradient(g, dt)
    #println("Zygote gradient (jd w.r.t. tdb seconds): $d_zyg")
    @test d_zyg > 0
    @test isapprox(d_zyg, 1/86400; rtol=1e-6)

    # Compare forward and reverse diff results
    @test isapprox(d_fwd, d_zyg; rtol=1e-9)

    # Round-trip sanity: adding 0 seconds leaves epoch unchanged
    t_same = add_seconds_in_scale(t0, 0.0, :tdb)
    @test t_same.scale == t0.scale
    @test t_same.format == t0.format
    @test t_same.jd == t0.jd
end
=#
nothing
