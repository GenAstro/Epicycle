# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Everything reaching the browser goes through `_publish`. Three properties, all
# of which exist because of live update rather than in spite of it:
#
#   atomic       — write a temporary file and rename it, so the dashboard can never fetch half a
#                  figure. Polling a file mid-write is otherwise a real, intermittent failure.
#   throttled    — a solver callback can publish thousands of times a second. Writes to one id
#                  coalesce, last value wins, and the pending one is flushed at exit.
#   non-blocking — publishing must not slow the run feeding it.
#
# This is also where transport changes. Swapping polled files for a websocket changes this file
# and the dashboard's fetch loop, and nothing a user types.

const MIN_INTERVAL = Ref(0.10)          # seconds between writes to one id

const _LAST    = Dict{String, Float64}()      # id => time of last write
const _PENDING = Dict{String, String}()       # id => payload waiting to go out
const _TIMERS  = Dict{String, Timer}()
const _WRITING = Task[]                       # writes started and possibly not yet on disk
const _LOCK    = ReentrantLock()

"""Write `payload` to `path` so a reader never sees it partly written."""
function atomic_write(path::AbstractString, payload::AbstractString)
    dir = dirname(path)
    tmp = joinpath(dir, ".tmp_" * basename(path) * "_" * string(getpid()) * "_" * string(rand(UInt32)))
    open(tmp, "w") do io
        write(io, payload)
    end
    mv(tmp, path; force = true)          # atomic on Windows and POSIX alike
    return path
end

"""
    _publish(id, payload)

Hand a payload over for delivery and return. Writes to the same id inside `MIN_INTERVAL`
coalesce, so a solver publishing every iteration costs one write per interval rather than one
per iteration.
"""
function _publish(id::AbstractString, payload::AbstractString)
    lock(_LOCK) do
        now  = time()
        last = get(_LAST, id, 0.0)

        if now - last >= MIN_INTERVAL[]
            _LAST[id] = now
            delete!(_PENDING, id)
            _start_write!(id, payload)
        else
            _PENDING[id] = payload                     # last value wins
            if !haskey(_TIMERS, id)
                wait_for = MIN_INTERVAL[] - (now - last)
                _TIMERS[id] = Timer(wait_for) do _
                    lock(_LOCK) do
                        delete!(_TIMERS, id)
                        p = pop!(_PENDING, id, nothing)
                        p === nothing && return
                        _LAST[id] = time()
                        _start_write!(id, p)
                    end
                end
            end
        end
    end
    return nothing
end

# Start a write without waiting for it, and keep the task so `flush_publishes` can. Finished
# tasks are dropped here, so a long run holds only the writes still in flight. Called under _LOCK.
function _start_write!(id::AbstractString, payload::AbstractString)
    filter!(!istaskdone, _WRITING)
    push!(_WRITING, @async _write_now(id, payload))
    return nothing
end

function _write_now(id::AbstractString, payload::AbstractString)
    try
        atomic_write(joinpath(assets_dir(), id), payload)
    catch e
        @warn "EpicycleIO: could not publish $id" exception = e
    end
end

"""
    flush_publishes()

Write out anything still waiting, and return once every publish made so far is on disk. Called
before the process exits so the last state of a live plot is never the one that got dropped.
"""
function flush_publishes()
    writing, pending = lock(_LOCK) do
        for (_, t) in _TIMERS
            close(t)
        end
        empty!(_TIMERS)
        w = copy(_WRITING)
        empty!(_WRITING)
        w, collect(pairs(_PENDING))
    end
    # Writes already started carry older values than anything still pending, so they land first
    # and the pending values overwrite them. Without the wait, a started write could land after
    # the flush and put an older value back.
    foreach(wait, writing)
    for (id, payload) in pending
        _write_now(id, payload)
    end
    lock(_LOCK) do
        empty!(_PENDING)
    end
    return nothing
end

# ─── What gets published ──────────────────────────────────────────────────────────────────────

# A scene is CZML, which Cesium wants under its own extension.
panel_file(p::Panel) = "panel_" * p.id * (p.family === :scene ? ".czml" : ".json")

"""Publish one panel's figure, and the manifest if the set of panels changed."""
function publish_panel(p::Panel)
    ensure_server()
    _publish(panel_file(p), payload(p))
    publish_manifest()
    return nothing
end

const _LAST_MANIFEST = Ref("")

function publish_manifest()
    m = JSON.json(manifest())
    m == _LAST_MANIFEST[] && return nothing
    _LAST_MANIFEST[] = m
    _publish("dashboard.json", m)
    return nothing
end
