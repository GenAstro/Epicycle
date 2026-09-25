# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A local static file server, on Sockets alone so the package takes no HTTP dependency.
#
# The dashboard cannot run from file:// — browsers block fetch there — so the served directory
# is a Scratch space the package writes into and the browser reads from over 127.0.0.1.

mutable struct ServerHandle
    port::Int
    socket::Sockets.TCPServer
    task::Task
end

const _SERVER  = Ref{Union{Nothing, ServerHandle}}(nothing)
const _ASSETS  = Ref{Union{Nothing, String}}(nothing)

# Whether anyone is watching is a question, not something to remember. The page polls several
# times a second, so a recent request means a live viewer and a long silence means the tab was
# closed. Remembering that we once opened one is what left a closed window unreopened.
const _LAST_REQUEST = Ref(0.0)
const _LAST_OPEN    = Ref(0.0)

const VIEWER_TIMEOUT = Ref(4.0)    # no request in this long: nobody is looking
const OPEN_GRACE     = Ref(15.0)   # a tab existed and went away: allow this long before another
const COLD_GRACE     = Ref(120.0)  # nothing has ever connected: a browser can be slow to start

"""Set false to stop plots opening a browser tab. On by default; tests and scripts turn it off."""
const AUTO_OPEN = Ref(true)

"""
    auto_open!(on::Bool)

Whether drawing opens a browser tab when none is watching. On by default.
"""
auto_open!(on::Bool) = (AUTO_OPEN[] = on)

"""The writable directory the dashboard is served from."""
function assets_dir()
    if _ASSETS[] === nothing
        dir = @get_scratch!("epicycle_io")
        _ASSETS[] = dir
    end
    stage_assets(_ASSETS[])
    return _ASSETS[]
end

# Copy packaged assets across when they are missing or older than the source. Skipping unchanged
# files matters: the dashboard polls modification times, so copying unconditionally would make it
# reload on every call.
function stage_assets(dir::AbstractString)
    src_dir = joinpath(@__DIR__, "assets")
    isdir(src_dir) || return dir
    for name in readdir(src_dir)
        src = joinpath(src_dir, name)
        isfile(src) || continue
        dst = joinpath(dir, name)
        if !isfile(dst) || mtime(src) > mtime(dst)
            cp(src, dst; force = true)
        end
    end
    return dir
end

const _MIME = Dict(".html" => "text/html; charset=utf-8",
                   ".js"   => "text/javascript; charset=utf-8",
                   ".css"  => "text/css; charset=utf-8",
                   ".json" => "application/json; charset=utf-8",
                   ".czml" => "application/json; charset=utf-8",
                   ".jpg"  => "image/jpeg",
                   ".png"  => "image/png",
                   ".glb"  => "model/gltf-binary")

content_type(path) = get(_MIME, lowercase(splitext(path)[2]), "application/octet-stream")

"""Start the server if it is not already up, and return its handle."""
function ensure_server()
    h = _SERVER[]
    h === nothing || return h

    dir = assets_dir()

    # The served directory outlives the process. Starting a session onto the last one's panels
    # would show plots this Julia knows nothing about, and they would never go away because
    # nothing here would ever republish them.
    sweep_panel_files()
    atomic_write(joinpath(dir, "dashboard.json"), "{\"panels\":[]}")
    atomic_write(joinpath(dir, NARRATION_FILE), "{}")      # nor an earlier session's caption
    sweep_narration_files()

    socket, port = _listen_on_free_port()
    task = @async _serve(socket, dir)
    h = ServerHandle(port, socket, task)
    _SERVER[] = h
    return h
end

# A fixed range, so a dashboard left open from an earlier session usually finds the same address
# and a bookmarked URL keeps working. Nothing is special about these numbers beyond being high,
# unassigned, and clear of the ports the common local dev servers take.
const PORT_RANGE = 8731:8780

function _listen_on_free_port()
    for port in PORT_RANGE
        try
            return Sockets.listen(Sockets.localhost, port), port
        catch e
            # Already taken — try the next one. Anything else is a real problem with the
            # machine's networking and should be seen rather than scanned past fifty times.
            e isa Base.IOError || rethrow()
        end
    end
    socket = Sockets.listen(Sockets.localhost, 0)      # let the OS pick
    return socket, Int(Sockets.getsockname(socket)[2])
end

function _serve(socket, dir)
    while isopen(socket)
        conn = try
            Sockets.accept(socket)
        catch e
            e isa Base.IOError || rethrow()
            break                       # socket closed while we were waiting
        end
        @async try
            _handle(conn, dir)
        catch e
            # A browser hanging up mid-response is ordinary and says nothing. Anything else is a
            # fault in our own handling, and swallowing it silently would make the page simply
            # stop working with nothing anywhere to say why.
            if !(e isa Base.IOError || e isa EOFError)
                @warn "EpicycleIO: failed to serve a request" exception = (e, catch_backtrace())
            end
        finally
            close(conn)
        end
    end
end

function _handle(conn, dir)
    _LAST_REQUEST[] = time()
    line = readline(conn)
    isempty(line) && return
    parts = split(line)
    length(parts) >= 2 || return
    target = parts[2]

    while true                          # discard headers
        h = readline(conn)
        (isempty(h) || h == "\r") && break
    end

    path = split(target, '?')[1]
    path = path == "/" ? "/dashboard.html" : path
    name = basename(path)               # no directory traversal: filename only

    file = joinpath(dir, name)
    if isfile(file)
        body = read(file)
        write(conn, "HTTP/1.1 200 OK\r\n",
                    "Content-Type: ", content_type(file), "\r\n",
                    "Content-Length: ", string(length(body)), "\r\n",
                    "Cache-Control: no-store\r\n",
                    "Access-Control-Allow-Origin: *\r\n\r\n")
        write(conn, body)
    else
        body = "not found: $name"
        write(conn, "HTTP/1.1 404 Not Found\r\n",
                    "Content-Type: text/plain\r\n",
                    "Content-Length: ", string(length(body)), "\r\n\r\n", body)
    end
end

"""
    dashboard_url(panels...)

The dashboard's address. Naming panels scopes the page to just those, which is how one browser
tab shows a subset without the panels themselves being anywhere different.
"""
function dashboard_url(titles::AbstractString...)
    base = "http://127.0.0.1:$(ensure_server().port)/dashboard.html"
    isempty(titles) && return base
    return base * "?panels=" * join(panel_id.(titles), ",")
end

"""
    open_dashboard()
    open_dashboard(panel, more...)

Open the dashboard in a browser tab.

With no arguments you get every panel. Naming panels opens a tab showing only those, which is
how you spread twelve plots across a few windows without splitting where they live:

```julia
open_dashboard("Range residuals", "Range-rate residuals")   # navigation, in its own tab
open_dashboard("Porkchop")                                  # one plot, full size
```

The panels are still one set on one server. A scoped tab is a view of it, so a panel that
several tabs show updates in all of them.

# Arguments
- `titles`: the panels to show. With none, the tab shows every panel.

# Notes
Opens a tab through the operating system's usual handler. If that fails — no browser, or a
headless machine — you get a warning carrying the URL rather than an error, because failing to
open a window is not a reason to stop a run.

# Returns
`nothing`.
"""
function open_dashboard(titles::AbstractString...)
    url = dashboard_url(titles...)
    try
        if Sys.iswindows()
            run(`cmd /c start "" $url`; wait = false)
        elseif Sys.isapple()
            run(`open $url`; wait = false)
        else
            run(`xdg-open $url`; wait = false)
        end
    catch e
        @warn "EpicycleIO: could not open a browser. Go to $url yourself." exception = e
    end
    _LAST_OPEN[] = time()
    return nothing
end

"""
Open a tab when nothing is watching.

Asks rather than remembers: the page polls several times a second, so a request in the last few
seconds means a viewer is alive. Closing the window therefore brings a new one back on the next
plot, which remembering that a tab was once opened did not.
"""
function wants_browser()
    now = time()
    now - _LAST_REQUEST[] < VIEWER_TIMEOUT[] && return false      # someone is looking

    if _LAST_REQUEST[] == 0.0
        # Nothing has ever connected, so a tab we opened is still starting. A cold browser can
        # take longer than a script does to run, and treating that silence as "no viewer" is
        # what opened a second tab onto the same plots.
        return now - _LAST_OPEN[] >= COLD_GRACE[]
    end

    # Something did connect and then stopped, so a tab existed and was closed.
    return now - _LAST_OPEN[] >= OPEN_GRACE[]
end

function ensure_browser()
    AUTO_OPEN[] || return nothing
    wants_browser() && open_dashboard()
    return nothing
end

"""
    close_dashboard()

Stop the local server.

An open tab stops updating and reports that it is disconnected. Your panels are untouched — they
live in Julia, not in the page — so the next plot starts the server again and a new tab picks up
where this one left off.

# Notes
Calling this when nothing is serving does nothing.

# Returns
`nothing`.

# Examples
```julia
close_dashboard()
```
"""
function close_dashboard()
    h = _SERVER[]
    h === nothing && return nothing
    try
        close(h.socket)
    catch e
        # Already closed is not a failure to close. Anything else is worth knowing about, even
        # though we go on to forget the handle either way.
        e isa Base.IOError || rethrow()
    end
    _SERVER[] = nothing
    _LAST_REQUEST[] = 0.0
    _LAST_OPEN[] = 0.0
    # The next server start sweeps the served directory and writes an empty manifest, so the
    # manifest this session already published is no longer on disk. Forget it, or the next publish
    # compares equal, skips the write, and the page loads an empty panel list.
    _LAST_MANIFEST[] = ""
    return nothing
end
