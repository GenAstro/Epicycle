# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using LinearAlgebra: norm

using AstroUniverse
using AstroUniverse: ensure_kernel_download, dafopr, dafbfs, daffna, dafgs,
                     dafus, dafcls

# The checksum tests fetch a real file, because what they are checking is the
# download path itself. Skip them rather than fail when there is no network.
have_network() = try
    mktempdir() do dir
        ensure_kernel_download(dir, "naif0012.tls", kernel_source("naif0012.tls").url)
        true
    end
catch
    false
end

@testset "kernel manifest is well formed" begin
    @test !isempty(KERNEL_SOURCES)

    for source in KERNEL_SOURCES
        @test !isempty(source.filename)
        @test !isempty(source.purpose)
        @test startswith(source.url, "https://")
        @test endswith(source.url, source.filename)

        # A checksum that is not 64 lowercase hex digits was mistyped, and would
        # reject every download of that kernel forever.
        @test length(source.sha256) == 64
        @test all(c -> c in "0123456789abcdef", source.sha256)
    end

    names = [s.filename for s in KERNEL_SOURCES]
    @test length(unique(names)) == length(names)

    # A default kernel missing from the manifest fails at `using AstroUniverse`,
    # which is the worst place to find out.
    for filename in AstroUniverse.DEFAULT_KERNELS
        @test kernel_source(filename) !== nothing
    end
end

@testset "kernel_source lookup" begin
    source = kernel_source("naif0012.tls")
    @test source isa KernelSource
    @test source.filename == "naif0012.tls"

    @test kernel_source("no_such_kernel.bsp") === nothing
end

@testset "download_spice_kernel names what it knows" begin
    message = try
        download_spice_kernel("no_such_kernel.bsp")
        ""
    catch e
        sprint(showerror, e)
    end

    @test occursin("no manifest entry", message)
    @test occursin("naif0012.tls", message)          # lists the alternatives
    @test occursin("download_spice_kernel(name, url)", message)
end

@testset "downloads are checked before they are cached" begin
    if !have_network()
        @info "No network; skipping checksum tests."
        @test true
    else
        source = kernel_source("naif0012.tls")   # 5 KB, so this is cheap

        # A good checksum caches the file.
        mktempdir() do dir
            path = ensure_kernel_download(dir, source.filename, source.url;
                                          sha256 = source.sha256)
            @test isfile(path)
            @test filesize(path) > 0
        end

        # A bad one caches nothing, and says both hashes.
        mktempdir() do dir
            wrong = "0"^64
            message = try
                ensure_kernel_download(dir, source.filename, source.url; sha256 = wrong)
                ""
            catch e
                sprint(showerror, e)
            end

            @test occursin("Checksum mismatch", message)
            @test occursin(wrong, message)
            @test occursin(source.sha256, message)

            # Nothing under the real name: a partial or wrong file left there
            # would be trusted by every later session, since the cache check is
            # `isfile` and cannot tell a bad kernel from a good one.
            @test !isfile(joinpath(dir, source.filename))
            @test isempty(readdir(dir))
        end

        # An already-cached kernel is not re-downloaded, so its checksum is not
        # re-examined. Deliberate: rehashing 94 MB on every startup would cost
        # more than it catches.
        mktempdir() do dir
            path = joinpath(dir, source.filename)
            write(path, "not a kernel")
            @test ensure_kernel_download(dir, source.filename, source.url;
                                         sha256 = source.sha256) == path
            @test read(path, String) == "not a kernel"
        end
    end
end

@testset "every planet center resolves from the default kernels" begin
    # The reason the merged kernel exists. Each `CelestialBody` names a body
    # center, and JPL's de440.bsp carries centers only for Mercury and Venus, so
    # against it these six raised a bare SPICE error instead of a position.
    jd = 2458849.5
    for body in (mercury, venus, mars, jupiter, saturn, uranus, neptune, pluto, moon, sun)
        r = translate(earth, body, jd)
        @test length(r) == 3
        @test all(isfinite, r)
        @test norm(r) > 0
    end

    # Both ids reach the same place: Mars center sits well under a metre from
    # the Mars barycenter, and both are in the file.
    mars_barycenter = CelestialBody("MarsBarycenter", mars.mu, mars.equatorial_radius,
                                    mars.flattening, 4)
    @test norm(translate(earth, mars, jd) - translate(earth, mars_barycenter, jd)) < 1e-3
end

@testset "the merged kernel carries what it claims" begin
    # The manifest checksum proves the bytes arrived intact. It cannot prove
    # the file was *built* right: rebuild the kernel with a body missing,
    # update the hash, and every download check still passes.
    #
    # A separate check proves the merged file reproduces its NAIF sources
    # exactly. It needs 3.6 GB of source kernels, so it cannot run here. This
    # is the part that can: walk the
    # segment descriptors and confirm the inventory and the span.
    #
    # It caught a real omission once already. The first build took Uranus from
    # `ura116xl.bsp`, which contains only the irregular satellites 716-724 and
    # no planet centre at all, and the file shipped without body 799.
    path = joinpath(get_spice_directory(), "epicycle_de440_1950-2100.bsp")
    @test isfile(path)

    bodies = Set{Int}()
    earliest, latest = Inf, -Inf
    handle = dafopr(path)
    dafbfs(handle)
    while daffna()
        times, ids = dafus(dafgs(5), 2, 6)
        push!(bodies, Int(ids[1]))
        earliest = min(earliest, times[1])
        latest   = max(latest,   times[2])
    end
    dafcls(handle)

    # Every barycentre, the Sun, Earth, Moon, and every planet centre. The
    # centres are the reason this file exists: `CelestialBody` names a body
    # centre, and JPL's de440.bsp carries none beyond Mercury and Venus.
    for naifid in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                   199, 299, 301, 399, 499, 599, 699, 799, 899, 999)
        naifid == 0 && continue          # the SSB is the root, not a segment
        @test naifid in bodies
    end

    # 1950 to 2100, as an ephemeris time. Checked from the inside so a rebuild
    # with a shorter window fails here rather than at whatever epoch a user
    # happens to try.
    et(year) = (year - 2000.0) * 365.25 * 86_400.0
    @test earliest <= et(1950.02)
    @test latest   >= et(2099.98)
end
