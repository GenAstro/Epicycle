# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# The local file server.
#
# It is a real HTTP server listening on 127.0.0.1, so it is tested by speaking
# HTTP to it rather than by calling `_handle` with a fake connection. That is the
# only way to find out whether the thing a browser talks to actually works.
#
# The case that matters most is directory traversal. The server hands out files
# from a scratch directory, and one line — `basename(path)` — is what stops
# `GET /../../secrets` walking out of it. Nothing else in the package would
# notice if that line were changed, and it is the sort of line that looks
# redundant to someone tidying up.
# =============================================================================

using Test
using EpicycleIO
using Sockets

const IO_ = EpicycleIO

IO_.auto_open!(false)

"""Speak HTTP to the server and return (status, headers, body).

The server closes the connection after responding, so the body is everything up
to EOF and no Content-Length parsing is needed.
"""
function _get(port::Integer, target::AbstractString; host = "127.0.0.1")
    sock = Sockets.connect(Sockets.localhost, port)
    try
        write(sock, "GET $target HTTP/1.1\r\nHost: $host\r\nConnection: close\r\n\r\n")
        raw = read(sock, String)
        head, _, body = partition(raw)
        lines  = split(head, "\r\n")
        status = parse(Int, split(lines[1])[2])
        headers = Dict{String, String}()
        for l in lines[2:end]
            i = findfirst(==(':'), l)
            i === nothing && continue
            headers[lowercase(strip(l[1:i-1]))] = strip(l[i+1:end])
        end
        return status, headers, body
    finally
        close(sock)
    end
end

# Split a raw response at the blank line between headers and body.
function partition(raw::AbstractString)
    i = findfirst("\r\n\r\n", raw)
    i === nothing && return raw, "", ""
    return raw[1:first(i)-1], "\r\n\r\n", raw[last(i)+1:end]
end

@testset "the local server" begin

h    = IO_.ensure_server()
port = h.port
dir  = IO_.assets_dir()

@testset "it serves the dashboard, and / means the dashboard" begin
    status, headers, body = _get(port, "/dashboard.html")
    @test status == 200
    @test occursin("text/html", headers["content-type"])
    @test !isempty(body)

    # A bare / is the dashboard, which is what the URL in the browser bar is.
    status_root, _, body_root = _get(port, "/")
    @test status_root == 200
    @test body_root == body
end

@testset "a missing file is 404, not a crash or a hang" begin
    status, headers, body = _get(port, "/no_such_panel.json")
    @test status == 404
    @test occursin("no_such_panel.json", body)
    @test occursin("text/plain", headers["content-type"])
end

@testset "a request cannot walk out of the served directory" begin
    # Plant a file next to the scratch directory, not inside it. If traversal
    # works, this is what leaks — and on a real machine the same request shapes
    # reach ssh keys and browser profiles.
    outside = joinpath(dirname(dir), "epicycle_io_traversal_probe.txt")
    write(outside, "SECRET")
    try
        @test isfile(outside)                       # the target really is there

        for target in ("/../epicycle_io_traversal_probe.txt",
                       "/../../epicycle_io_traversal_probe.txt",
                       "/..%2Fepicycle_io_traversal_probe.txt",
                       "/subdir/../../epicycle_io_traversal_probe.txt")
            status, _, body = _get(port, target)
            @test status == 404
            @test !occursin("SECRET", body)
        end

        # An absolute path is not a way in either.
        status, _, body = _get(port, "/C:/Windows/win.ini")
        @test status == 404
        @test !occursin("SECRET", body)
    finally
        rm(outside; force = true)
    end
end

@testset "a query string selects nothing but is not part of the filename" begin
    # The page appends ?t=… to defeat caching, and scoped tabs append ?panels=…
    # Both have to be stripped before the file is looked up.
    probe = joinpath(dir, "server_probe.json")
    write(probe, "{\"ok\":true}")
    try
        status, headers, body = _get(port, "/server_probe.json?t=12345")
        @test status == 200
        @test body == "{\"ok\":true}"
        @test occursin("application/json", headers["content-type"])

        status2, _, body2 = _get(port, "/server_probe.json")
        @test status2 == 200 && body2 == body
    finally
        rm(probe; force = true)
    end
end

@testset "content types cover what the dashboard loads" begin
    @test occursin("text/html",       IO_.content_type("dashboard.html"))
    @test occursin("text/javascript", IO_.content_type("app.js"))
    @test occursin("text/css",        IO_.content_type("style.css"))
    @test occursin("application/json", IO_.content_type("panel_A.json"))
    @test occursin("application/json", IO_.content_type("scene.czml"))
    @test IO_.content_type("earth.jpg") == "image/jpeg"
    @test IO_.content_type("earth.PNG") == "image/png"      # case does not matter
    @test IO_.content_type("Aqua.glb")  == "model/gltf-binary"

    # Anything unrecognised still downloads rather than being guessed at.
    @test IO_.content_type("notes.txt") == "application/octet-stream"
    @test IO_.content_type("noextension") == "application/octet-stream"
end

@testset "serving a request is what proves a viewer is alive" begin
    # `wants_browser` asks whether anyone polled recently. That question is only
    # answerable because _handle stamps the clock, so the two are tested together.
    IO_._LAST_REQUEST[] = 0.0
    _get(port, "/dashboard.html")
    @test IO_._LAST_REQUEST[] > 0.0
    @test time() - IO_._LAST_REQUEST[] < IO_.VIEWER_TIMEOUT[]
    @test !IO_.wants_browser()                  # a live viewer: do not open another
end

@testset "the url names the port, and scoping is a query" begin
    url = IO_.dashboard_url()
    @test occursin("127.0.0.1:$port", url)
    @test endswith(url, "/dashboard.html")

    scoped = IO_.dashboard_url("Range residuals", "Porkchop")
    @test occursin("?panels=", scoped)
    # Panel ids, not titles — the same mapping the files on disk use.
    @test occursin(IO_.panel_id("Range residuals"), scoped)
    @test occursin("Porkchop", scoped)
end

@testset "assets are staged, and staging again does not touch them" begin
    # The page polls modification times, so re-copying unchanged files on every
    # call would make the dashboard reload continuously.
    dash = joinpath(dir, "dashboard.html")
    @test isfile(dash)
    before = mtime(dash)
    IO_.stage_assets(dir)
    @test mtime(dash) == before
end

@testset "closing lets the port go, and reopening gets one again" begin
    IO_.close_dashboard()
    @test IO_._SERVER[] === nothing
    @test IO_._LAST_REQUEST[] == 0.0

    # A closed server is not a dead package: the next plot brings one back.
    h2 = IO_.ensure_server()
    @test h2 !== nothing
    status, _, _ = _get(h2.port, "/dashboard.html")
    @test status == 200

    # Calling it twice hands back the same server rather than starting a second.
    @test IO_.ensure_server() === h2
end

IO_.close_dashboard()

end # testset
