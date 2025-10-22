using Test
using SPICE

using AstroUniverse
using AstroEpochs

# Helper: are SPICE kernels loaded?
has_kernels() = try
    SPICE.ktotal("ALL") > 0
catch
    false
end

@testset "CelestialBody constructors" begin
    # Positional constructor
    b = CelestialBody("X", 1.0, 2.0, 0.1, 123)
    @test b.name == "X"
    @test b.mu == 1.0
    @test b.equatorial_radius == 2.0
    @test b.flattening == 0.1
    @test b.naifid == 123

    # Keyword constructor defaults to EARTH_DEFAULTS
    d = CelestialBody()
    @test d.name == "unnamed"
    @test d.mu == AstroUniverse.EARTH_DEFAULTS.mu
    @test d.equatorial_radius == AstroUniverse.EARTH_DEFAULTS.equatorial_radius
    @test d.flattening == AstroUniverse.EARTH_DEFAULTS.flattening
    @test d.naifid == AstroUniverse.EARTH_DEFAULTS.naifid

    # Type promotion (BigFloat)
    bb = CelestialBody(name="B", mu=big(1.0), equatorial_radius=2.0, flattening=0.5, naifid=1)
    @test bb.mu isa BigFloat
    @test bb.equatorial_radius isa BigFloat
    @test bb.flattening isa BigFloat
    @test bb.naifid isa Int
end

@testset "CelestialBody show methods" begin
    b = CelestialBody("ShowMe", 3.0, 4.0, 0.2, 42)

    # MIME text/plain
    s_plain = sprint(show, MIME"text/plain"(), b)
    @test occursin("CelestialBody:", s_plain)
    @test occursin("name               = ShowMe", s_plain)
    @test occursin("NAIF ID            = 42", s_plain)

    # Generic delegator
    s_generic = sprint(show, b)
    @test s_generic == s_plain
end

@testset "Canonical body instances" begin
    @test earth.name == "Earth"
    @test earth.naifid == 399
    @test moon.name == "Moon"
    @test moon.naifid == 301
    @test sun.name == "Sun"
    @test sun.naifid == 10
    @test mars.name == "Mars"
    @test mars.naifid == 499
    @test venus.name == "Venus"
    @test venus.naifid == 299
    @test mercury.name == "Mercury"
    @test mercury.naifid == 199
    @test jupiter.name == "Jupiter"
    @test jupiter.naifid == 599
    @test saturn.name == "Saturn"
    @test saturn.naifid == 699
    @test uranus.name == "Uranus"
    @test uranus.naifid == 799
    @test neptune.name == "Neptune"
    @test neptune.naifid == 899
    @test pluto.name == "Pluto"
    @test pluto.naifid == 999
end

@testset "translate (calls work - no truth data)" begin
    if has_kernels()
        jd = 2458018.0

        # Returns a 3-element vector
        r_em = translate(earth, moon, jd)
        @test length(r_em) == 3

        # Same body → zero vector
        r_ee = translate(earth, earth, jd)
        @test r_ee == zeros(length(r_ee))
    else
        @info "SPICE kernels not loaded; skipping translate tests."
        @test true  # keep a passing test to not mark set as empty
    end
end

@testset "translate test values against GMAT" begin

    # Earth w/r/t Sun at 2015-09-21T12:23:12.000 TDB
    t = Time("2015-09-21T12:23:12.000", TDB(), ISOT())
    # pos,vel from GMAT R2022a
    pos = [150111759.0438753, -4823271.352340085, -2092193.447543865];
    #vel = [0.5463431391132124, 27.21155496102735, 11.7965546582813];
    r_em = translate(sun, earth, t.jd)

    tol = 1e-10
    @test isapprox(r_em[1], pos[1]; rtol=tol)
    @test isapprox(r_em[2], pos[2]; rtol=tol)
    @test isapprox(r_em[3], pos[3]; rtol=tol)

    # pos,vel from GMAT R2022a
    t = Time("1995-05-10T00:00:00.000", TDB(), ISOT())
    pos = [365758.2395306571, -107619.1828216814, -16569.18869984438];
    vel = [0.2181927510507531, 0.9448998938751025, 0.3403413852066187];
    r_em = translate(moon, earth, t.jd)

    tol = 1e-10
    @test isapprox(r_em[1], pos[1]; rtol=tol)
    @test isapprox(r_em[2], pos[2]; rtol=tol)
    @test isapprox(r_em[3], pos[3]; rtol=tol)

end

@testset "CelestialBody constructor validation" begin
    # Valid flattening values
    @test CelestialBody("B1", 398600.0, 6378.0, 0.0, 1) isa CelestialBody
    @test CelestialBody("B2", 398600.0, 6378.0, 0.5, 2) isa CelestialBody

    # Helper to capture error messages
    _msg(f) = try
        f()
        nothing
    catch e
        sprint(showerror, e)
    end

    # μ must be finite and > 0
    for bad in (0.0, -1.0, Inf, NaN)
        m = _msg(() -> CelestialBody("X", bad, 6378.0, 0.1, 1))
        @test occursin("μ must be finite and > 0", m)
    end

    # equatorial_radius must be finite and > 0
    for bad in (0.0, -1.0, Inf, NaN)
        m = _msg(() -> CelestialBody("X", 398600.0, bad, 0.1, 1))
        @test occursin("equatorial_radius must be finite and > 0", m)
    end

    # flattening must be finite and in [0, 1)
    for bad in (-1.0, 1.0, Inf, NaN)
        m = _msg(() -> CelestialBody("X", 398600.0, 6378.0, bad, 1))
        @test occursin("flattening must be finite and in [0, 1)", m)
    end

    # Keyword constructor mirrors checks
    @test occursin("μ must be finite and > 0",
                   _msg(() -> CelestialBody(name="X", mu=0.0, equatorial_radius=6378.0, flattening=0.1, naifid=1)))
    @test occursin("equatorial_radius must be finite and > 0",
                   _msg(() -> CelestialBody(name="X", mu=398600.0, equatorial_radius=0.0, flattening=0.1, naifid=1)))
    @test occursin("flattening must be finite and in [0, 1)",
                   _msg(() -> CelestialBody(name="X", mu=398600.0, equatorial_radius=6378.0, flattening=1.0, naifid=1)))
end


@testset "gravparam API: get_gravparam/set_gravparam!" begin
    # Construct a simple body
    body = CelestialBody("TestBody", 100.0, 1000.0, 0.0, 42)

    # get_gravparam returns μ
    @test get_gravparam(body) == 100.0

    # set_gravparam! updates μ and preserves the field's element type
    Tμ = typeof(body.mu)
    rv = set_gravparam!(body, BigFloat(200))  # pass different numeric type
    @test rv === body
    @test body.mu == 200.0
    @test typeof(body.mu) == Tμ

    # Subsequent update works
    @test set_gravparam!(body, 150.0) === body
    @test body.mu == 150.0

    # Invalid μ values throw
    @test_throws ArgumentError set_gravparam!(body, 0.0)
    @test_throws ArgumentError set_gravparam!(body, -1.0)
    @test_throws ArgumentError set_gravparam!(body, Inf)
    @test_throws ArgumentError set_gravparam!(body, NaN)
end

nothing