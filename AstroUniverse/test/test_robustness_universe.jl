# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The production audit's findings, each pinned where it was found: wrong answers that used to come
# back without a word, and failures that used to give no way forward.
#
# Truth: the raw IERS and NAIF files where a number is asserted, named at each test.

using Test
using AstroUniverse
using AstroEpochs
using ForwardDiff
using LinearAlgebra
using SHA

const _AU = AstroUniverse

_msg(f) = try
    f()
    ""
catch e
    sprint(showerror, e)
end

const _T_R = Time(2458849.5, TDB(), JD())

# --- The rows of both series for 2024-01-01, as IERS and USNO publish them -------------------------

const _CSV_HEAD = "MJD;Year;Month;Day;Type;x_pole;sigma_x_pole;y_pole;sigma_y_pole;x_rate;" *
    "sigma_x_rate;y_rate;sigma_y_rate;Type;UT1-UTC;sigma_UT1-UTC;LOD;sigma_LOD;Type;dPsi;" *
    "sigma_dPsi;dEpsilon;sigma_dEpsilon;dX;sigma_dX;dY;sigma_dY;Type;bulB/x_pole;bulB/y_pole;" *
    "Type;bulB/UT-UTC;Type;bulB/dPsi;bulB/dEpsilon;bulB/dX;bulB/dY"
const _CSV_1980 = "60310;2024;01;01;final;0.136912;0.000012;0.202190;0.000019;;;;;final;0.0087837;" *
    "0.0000084;0.2375;0.0072;final;-110.440;0.800;-5.644;0.140;;;;;final;0.136894;0.202185;final;" *
    "0.0087572;final;-110.452;-5.772;;"
const _CSV_2000A = "60310;2024;01;01;final;0.136912;0.000012;0.202190;0.000019;;;;;final;0.0087837;" *
    "0.0000084;0.2375;0.0072;final;;;;;0.295;0.318;-0.095;0.140;final;0.136894;0.202185;final;" *
    "0.0087572;final;;;0.283;-0.183"

_second_day(row) = replace(row, "60310;2024;01;01" => "60311;2024;01;02")

function _csv_file(row)
    path = tempname() * ".csv"
    write(path, join((_CSV_HEAD, row, _second_day(row)), "\n") * "\n")
    return path
end

function _finals_line(mjd, day, dpsi_or_dx, deps_or_dy)
    line = fill(' ', 134)
    put!(r, s) = (line[r] = collect(lpad(s, length(r))); nothing)
    put!(1:2, "24"); put!(3:4, "1"); put!(5:6, string(day)); put!(8:15, string(mjd, ".00"))
    line[17] = 'I'; put!(19:27, "0.136912"); put!(28:36, "0.000012")
    put!(38:46, "0.202190"); put!(47:55, "0.000019"); line[58] = 'I'
    put!(59:68, "0.0087837"); put!(69:78, "0.0000084"); put!(80:86, "0.2375"); put!(87:93, "0.0072")
    line[96] = 'I'; put!(98:106, dpsi_or_dx); put!(107:115, "0.300")
    put!(117:125, deps_or_dy); put!(126:134, "0.140")
    return String(line)
end

function _finals_file(a, b)
    path = tempname()
    write(path, _finals_line(60310, 1, a, b) * "\n" * _finals_line(60311, 2, a, b) * "\n")
    return path
end

@testset "eop_load refuses a file of the other series" begin
    # The fixed-width files share a layout, so this used to load the wrong numbers silently:
    # finals.all read for IAU2006 gave dX = -110.44 mas where the true value is 0.295.
    finals_1980, finals_2000a = _finals_file("-110.440", "-5.644"), _finals_file("0.295", "-0.095")
    message = _msg(() -> eop_load(finals_1980; theory = IAU2006()))
    @test occursin("IAU 1980", message) && occursin("theory = FK5()", message)
    message = _msg(() -> eop_load(finals_2000a; theory = FK5()))
    @test occursin("IAU 2000A", message) && occursin("theory = IAU2006()", message)

    # The CSV forms used to fail, saying the file was corrupted.
    message = _msg(() -> eop_load(_csv_file(_CSV_1980); theory = IAU2006()))
    @test occursin("IAU 1980", message) && occursin("theory = FK5()", message)
    message = _msg(() -> eop_load(_csv_file(_CSV_2000A); theory = FK5()))
    @test occursin("IAU 2000A", message) && occursin("theory = IAU2006()", message)

    # Each file of the right series loads, with its own values.
    jd = 2460310.5
    saved = (eop(FK5()), eop(IAU2006()))
    try
        @test eop_load(finals_1980; theory = FK5()).δΔψ(jd) ≈ -110.440
        @test eop_load(finals_2000a; theory = IAU2006()).δx(jd) ≈ 0.295
        @test eop_load(_csv_file(_CSV_1980); theory = FK5()).δΔψ(jd) ≈ -110.440
        @test eop_load(_csv_file(_CSV_2000A); theory = IAU2006()).δx(jd) ≈ 0.295
    finally
        foreach(set_eop!, saved)
    end
end

@testset "eop_load runs the checks a download gets" begin
    path = tempname() * ".csv"
    write(path, join((_CSV_HEAD, replace(_CSV_2000A, "0.136912" => "Inf"), _second_day(_CSV_2000A)), "\n") * "\n")
    message = _msg(() -> eop_load(path; theory = IAU2006()))
    @test occursin("not a finite number", message) && occursin(path, message)
end

@testset "a keyword-built body is not an Earth" begin
    # It used to take NAIF 399 and Earth's flattening, so `translate(earth, b, t)` returned zeros.
    @test_throws UndefKeywordError CelestialBody(name = "MyAsteroid", mu = 1e-6, equatorial_radius = 1.0)
    @test_throws UndefKeywordError CelestialBody(mu = 1e-6, equatorial_radius = 1.0, naifid = 2000433)
    @test_throws UndefKeywordError CelestialBody()
    b = CelestialBody(name = "MyAsteroid", mu = 1e-6, equatorial_radius = 1.0, naifid = 2000433)
    @test b.flattening == 0.0
    @test b.texture_file == ""
end

@testset "ephemeris calls take a Time" begin
    utc = Time("2024-01-01T00:00:00", UTC(), ISOT())
    tdb = utc.tdb

    # The call reads the epoch in TDB whatever scale it was given in.
    for t in (utc, utc.tt, utc.tai, tdb)
        @test translate(sun, earth, t) ≈ translate(sun, earth, tdb) atol = 1e-3
        @test translate_state(sun, earth, t) ≈ translate_state(sun, earth, tdb) atol = 1e-3
    end

    # A bare date, the source of a silent 2,095 km error from a UTC date, is refused.
    @test occursin("Time(jd, TDB(), JD())", _msg(() -> translate(sun, earth, utc.jd)))
    @test occursin("Time(jd, TDB(), JD())", _msg(() -> translate_state(sun, earth, utc.jd)))
    @test occursin("Time(jd, TDB(), JD())", _msg(() -> iau2015_orientation(mars, utc.jd)))
    @test occursin("Time(jd, TDB(), JD())", _msg(() -> body_axes_rotation(IAU2015(), 499, utc.jd)))

    # Truth: SPICE at the ephemeris time built from the two parts. A single date is 40 µs coarse here.
    t = Time(2460463.0, 0.0240740741, TDB(), JD())
    et = ((t.jd1 - 2451545.0) + t.jd2) * 86400.0
    @test translate(sun, earth, t) == _AU.spkpos("399", et, "J2000", "NONE", "10")[1]
    @test translate_state(sun, earth, t) == _AU.spkezr("399", et, "J2000", "NONE", "10")[1]

    # Clear messages at the SPICE boundary, where they used to be a conversion error from deep inside.
    @test occursin("not differentiable",
                   _msg(() -> ForwardDiff.derivative(d -> translate(earth, moon, Time(d, TDB(), JD()))[1], 2458849.5)))
    @test occursin("not finite", _msg(() -> translate(earth, moon, Time(NaN, TDB(), JD()))))
end

@testset "ephemeris truth for each body" begin
    # Truth: JPL Horizons, 2020-01-01 TDB, solar system barycenter to body center, as checked in the
    # production audit. Exact agreement with AstroFrames' regression rows is checked there; this is
    # the barycenter against center check an outer planet needs.
    t = Time(2458849.5, TDB(), JD())
    offset(body, barycenter_id) = norm(_AU.spkpos(string(body.naifid), _ephemeris_et(t), "J2000", "NONE", string(barycenter_id))[1])
    _ephemeris_et(t) = (t.jd1 - 2451545.0 + t.jd2) * 86400.0
    @test offset(jupiter, 5) ≈ 110.114 atol = 0.01
    @test offset(saturn, 6) ≈ 302.966 atol = 0.01
    @test offset(pluto, 9) ≈ 2131.033 atol = 0.01
    @test offset(mars, 4) < 0.001

    # Every body resolves at both ends of the default kernel's span.
    for jd in (2433283.0, 2488069.0), body in (mercury, venus, moon, mars, jupiter, saturn, uranus, neptune, pluto, sun)
        @test all(isfinite, translate(earth, body, Time(jd, TDB(), JD())))
    end

    # Velocity is the time derivative of position.
    h = 10.0 / 86400
    for body in (moon, mars, jupiter, neptune)
        ahead = translate(earth, body, Time(t.jd1, t.jd2 + h, TDB(), JD()))
        behind = translate(earth, body, Time(t.jd1, t.jd2 - h, TDB(), JD()))
        state = translate_state(earth, body, t)
        @test state[1:3] == translate(earth, body, t)
        @test norm(state[4:6] - (ahead - behind) / 20.0) < 1e-7
    end
end

@testset "built-in constants against their sources" begin
    # Truth: NAIF gm_de440.tpc (BODY<n>_GM) and pck00011.tpc (BODY<n>_RADII), copied from the files.
    # A planet with moons carries the system μ, since these are set up for Earth orbits.
    table = ((sun, 1.3271244004127942e11, 695700.0, 695700.0),
             (mercury, 2.2031868551400003e4, 2440.53, 2438.26),
             (venus, 3.24858592e5, 6051.8, 6051.8),
             (moon, 4.9028001184575496e3, 1737.4, 1737.4),
             (mars, 4.2828375815756102e4, 3396.19, 3376.20),
             (jupiter, 1.2671276409999998e8, 71492.0, 66854.0),
             (saturn, 3.7940584841799997e7, 60268.0, 54364.0),
             (uranus, 5.7945563999999985e6, 25559.0, 24973.0),
             (neptune, 6.8365271005803989e6, 24764.0, 24341.0),
             (pluto, 9.755e2, 1188.3, 1188.3))
    for (body, mu, equator, pole) in table
        @test body.mu == mu
        @test body.equatorial_radius == equator
        @test body.equatorial_radius * (1 - body.flattening) ≈ pole rtol = 1e-14
    end
    # Truth: WGS84 (NGA TR8350.2).
    @test earth.mu == 398600.4418
    @test earth.equatorial_radius == 6378.137
    @test earth.flattening == 1 / 298.257223563
end

@testset "built-in bodies are constants" begin
    for name in (:sun, :mercury, :venus, :earth, :moon, :mars, :jupiter, :saturn, :uranus, :neptune, :pluto)
        @test isconst(_AU, name)
    end
    # The texture path is the installed package's, found when it loads.
    @test isfile(earth.texture_file)
    @test startswith(earth.texture_file, pkgdir(AstroUniverse))
    @test pluto.texture_file == ""
end

@testset "μ has one setter" begin
    # Four routes used to behave four ways.
    b() = CelestialBody("Body", 1.0, 1.0, 0.0, 1000002)
    @test_throws ArgumentError set_gravparam!(b(), -1.0)
    @test_throws ArgumentError AstroUniverse.set_field!(b(), Mu(), -1.0)
    @test_throws ArgumentError (b().mu = -5.0)
    @test_throws ArgumentError (b().mu = NaN)
    @test_throws ArgumentError (b().equatorial_radius = 0.0)
    @test_throws ArgumentError (b().flattening = 1.0)

    # A Float32 or integer body holds Float64, so a value neither underflows nor is refused as inexact.
    @test CelestialBody("F32", 1f0, 1f0, 0f0, 1000003) isa CelestialBody{Float64}
    integer = CelestialBody("Int", 1, 1, 0, 1000004)
    @test integer isa CelestialBody{Float64}
    @test set_gravparam!(integer, 1.5).mu == 1.5
    @test set_gravparam!(CelestialBody("Tiny", 1f0, 1f0, 0f0, 1000005), 1e-100).mu == 1e-100
end

@testset "a body prints on one line inside anything else" begin
    @test sprint(show, sun) == "CelestialBody(\"Sun\", naifid = 10)"
    @test "$(moon)" == "CelestialBody(\"Moon\", naifid = 301)"
    block = sprint(show, MIME"text/plain"(), earth)
    @test occursin("km³/s²", block) && occursin(" km", block)
end

@testset "derivatives pass through IAU2015 for all nine bodies" begin
    # They used to throw for seven of them.
    for b in (sun, mercury, venus, mars, jupiter, saturn, uranus, neptune, pluto)
        f(d) = body_axes_rotation(IAU2015(), b.naifid, Time(d, TDB(), JD()))[1, 1]
        jd, h = 2458849.5, 1e-3
        @test ForwardDiff.derivative(f, jd) ≈ (f(jd + h) - f(jd - h)) / 2h rtol = 1e-4
    end
    @test iau2015_orientation(sun, Time(big"2458849.5", TDB(), JD())) isa IAU2015Orientation{BigFloat}
end

@testset "an estimated EOP offset carries its derivative" begin
    # A dual-valued offset at a plain date used to throw.
    struct _PoleOffset{T} <: AbstractEopCorrection
        dx::T
    end
    AstroUniverse.eop_correction(m::_PoleOffset, jd_utc, jd_tt) = (x = m.dx, y = 0.0, ut1_utc = 0.0, lod = 0.0)
    try
        @test ForwardDiff.derivative(0.001) do dx
            set_eop_corrections!([_PoleOffset(dx)])
            eop(IAU2006(), 2460000.5, 2460000.5008).x
        end == 1.0
    finally
        set_eop_corrections!([])
    end
end

@testset "corrections reach every field, whichever date carries the derivative" begin
    struct _AllFields <: AbstractEopCorrection end
    AstroUniverse.eop_correction(::_AllFields, jd_utc, jd_tt) =
        (x = 0.001, y = 0.002, ut1_utc = 0.003, lod = 0.004 * (jd_tt - jd_tt + 1))
    jd_utc, jd_tt = 2460000.5, 2460000.5008
    base = eop(IAU2006(), jd_utc, jd_tt)
    try
        set_eop_corrections!([_AllFields()])
        v = eop(IAU2006(), jd_utc, jd_tt)
        @test v.x == base.x + 0.001
        @test v.y == base.y + 0.002
        @test v.ut1_utc == base.ut1_utc + 0.003
        @test v.lod == base.lod + 0.004
        @test v.celestial_pole == base.celestial_pole
        # A derivative in TT alone, with a plain UTC date.
        @test ForwardDiff.derivative(d -> eop(IAU2006(), jd_utc, d).lod, jd_tt) == 0.0
    finally
        set_eop_corrections!([])
    end
end

@testset "eop at a leap second" begin
    # Truth: UT1 − TAI is continuous by definition, and it used to jump a second through 23:59:60.
    ut1_tai(s) = (t = Time(s, UTC(), ISOT()); eop(IAU2006(), t).ut1_utc - tai_minus_utc(t))
    before = ut1_tai("2016-12-31T23:59:59.500")
    @test ut1_tai("2016-12-31T23:59:60.000") ≈ before atol = 1e-3
    @test ut1_tai("2016-12-31T23:59:60.500") ≈ before atol = 1e-3
    @test ut1_tai("2017-01-01T00:00:00.000") ≈ before atol = 1e-3
end

@testset "SPICE from threads" begin
    # CSPICE is not thread-safe, and an unguarded call from threads aborts the process, so this has
    # to run in its own.
    script = """
    using AstroUniverse, AstroEpochs
    bodies = (sun, earth, moon, mars, jupiter)
    results = zeros(2000)
    Threads.@threads :dynamic for i in eachindex(results)
        t = Time(2458849.5 + i, TDB(), JD())
        results[i] = translate(bodies[mod1(i, 5)], bodies[mod1(i + 1, 5)], t)[1]
    end
    expected = [translate(bodies[mod1(i, 5)], bodies[mod1(i + 1, 5)], Time(2458849.5 + i, TDB(), JD()))[1]
                for i in eachindex(results)]
    println(results == expected ? "SAME" : "DIFFERENT")
    """
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) -t 4 -e $script`
    output = IOBuffer()
    process = run(pipeline(ignorestatus(cmd); stdout = output, stderr = devnull))
    @test process.exitcode == 0
    @test occursin("SAME", String(take!(output)))
end

@testset "kernel download" begin
    _file_url(path) = "file:///" * replace(lstrip(replace(path, '\\' => '/'), '/'), " " => "%20")
    mktempdir() do dir
        source = joinpath(dir, "source.bsp")
        write(source, "the kernel's bytes")
        sha = bytes2hex(open(sha256, source))
        size = filesize(source)
        cache = mkpath(joinpath(dir, "cache"))
        url = _file_url(source)

        # A truncated cached file used to be trusted forever.
        write(joinpath(cache, "k.bsp"), "short")
        path = @test_logs (:warn, r"should be") _AU.ensure_kernel_download(cache, "k.bsp", url; sha256 = sha, size = size)
        @test read(path, String) == "the kernel's bytes"

        # The right size is not fetched again.
        @test _AU.ensure_kernel_download(cache, "k.bsp", "file:///no/such/file"; sha256 = sha, size = size) == path

        # Another process's temporary file is not this attempt's, and nothing is left behind.
        write(joinpath(cache, "j.bsp.part"), "someone else's")
        _AU.ensure_kernel_download(cache, "j.bsp", url; sha256 = sha, size = size)
        @test read(joinpath(cache, "j.bsp.part"), String) == "someone else's"
        @test sort(readdir(cache)) == ["j.bsp", "j.bsp.part", "k.bsp"]
    end
end

@testset "processes fetching one kernel at once" begin
    # They shared one temporary file, and all four failed in twelve of twelve trials.
    mktempdir() do dir
        source = joinpath(dir, "source.bsp")
        write(source, rand(UInt8, 5_000_000))
        size = filesize(source)
        sha = bytes2hex(open(sha256, source))
        cache, go = mkpath(joinpath(dir, "cache")), joinpath(dir, "go")
        url = "file:///" * lstrip(replace(source, '\\' => '/'), '/')
        script = """
        using AstroUniverse
        cache, go, url, sha, size = ARGS[1], ARGS[2], ARGS[3], ARGS[4], parse(Int, ARGS[5])
        while !isfile(go); sleep(0.01); end
        AstroUniverse.ensure_kernel_download(cache, "race.bsp", url; sha256 = sha, size = size)
        """
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $script $cache $go $url $sha $size`
        workers = [run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull); wait = false) for _ in 1:4]
        sleep(8.0)                      # long enough for each to load the package and reach the wait
        touch(go)
        foreach(wait, workers)
        @test all(p -> p.exitcode == 0, workers)
        @test readdir(cache) == ["race.bsp"]
        @test read(joinpath(cache, "race.bsp")) == read(source)
    end
end

@testset "loading needs no network, and the failure says what to do" begin
    # Run with an empty depot first on the path and no route out, so the kernels cannot be fetched.
    mktempdir() do depot
        script = """
        using AstroUniverse, AstroEpochs
        println("MU ", earth.mu)
        println("LOADED ", isempty(loaded_spice_kernels()))
        try
            translate(earth, moon, Time(2458849.5, TDB(), JD()))
            println("NO ERROR")
        catch e
            println("MESSAGE ", sprint(showerror, e))
        end
        """
        env = copy(ENV)
        env["JULIA_DEPOT_PATH"] = depot * (Sys.iswindows() ? ";" : ":") * join(DEPOT_PATH, Sys.iswindows() ? ";" : ":")
        env["HTTPS_PROXY"] = "http://127.0.0.1:9"
        env["HTTP_PROXY"] = "http://127.0.0.1:9"
        cmd = setenv(`$(Base.julia_cmd()) --project=$(Base.active_project()) --compiled-modules=existing --startup-file=no -e $script`, env)
        output = IOBuffer()
        process = run(pipeline(ignorestatus(cmd); stdout = output, stderr = devnull))
        text = String(take!(output))
        @test process.exitcode == 0
        @test occursin("MU 398600.4418", text)
        @test occursin("LOADED true", text)
        @test occursin("naif0012.tls", text)
        @test occursin("download each file from the address shown", text)
        @test occursin("spice_kernels", text)
    end
end

@testset "unloading a kernel loaded from elsewhere" begin
    # The name was tested and the cache path unloaded, so this used to hang forever.
    mktempdir() do dir
        name = "pck00011.tpc"
        download_spice_kernel(name)
        copy = cp(joinpath(get_spice_directory(), name), joinpath(dir, name))
        _AU.furnsh(copy)
        @test name in loaded_spice_kernels()
        unload_spice_kernel(name)
        @test !(name in loaded_spice_kernels())

        # A different file that merely shares the name is loaded, not skipped.
        load_spice_kernel(name)
        _AU.furnsh(copy)
        @test count(==(name), loaded_spice_kernels()) == 2
        @test_logs (:warn, r"2 times") unload_spice_kernel(name)
        @test !(name in loaded_spice_kernels())
    end
end

@testset "kernel lists" begin
    @test loaded_spice_kernels() isa Vector{String}
    @test "naif0012.tls" in loaded_spice_kernels()
    cached = cached_spice_kernels()
    @test cached isa Vector{String}
    # Every kernel type, which used to hide .bpc, .tf and .tpc.
    @test "moon_pa_de440_200625.bpc" in cached
    @test "moon_de440_250416.tf" in cached
    @test "naif0012.tls" in cached
end
