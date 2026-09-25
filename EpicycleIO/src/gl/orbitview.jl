# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A trajectory drawn in three dimensions around a central body, with a clock.
#
# This is the one part of EpicycleIO that knows what a spacecraft is, and it has to: a
# trajectory on a globe is meaningless without a frame and an epoch. `xyplot` takes arrays and
# knows nothing; `orbitview` takes the spacecraft.
#
# Earth is the only central body today. Nothing here names it — the body arrives from the
# frame the history was recorded in, and adding the Moon or a planet is a setting rather than
# a rewrite.

"""
    orbitview([panel], sat; kwargs...)

Draw a spacecraft's propagated trajectory in three dimensions, with a clock you can play.

The panel joins the dashboard beside any plots, so a trajectory and the quantities describing
it share one page.

Everything is read from the spacecraft:

- the trajectory, from `sat.history`
- the epoch, from the earliest recorded sample
- one coloured arc per propagation segment, from the same palette the plots use
- the name shown on the entity, from `sat.name`

```julia
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))

orbitview(sat)
orbitview("Mission", sat)          # name the panel
```

# Arguments
- `panel`: optional leading string naming the panel.
- `sat`: the spacecraft. It must have been propagated — the trajectory is read from its history.

# Keyword arguments
- `colors`: palette for the arcs, cycled. Defaults to the shared palette, so a trajectory and
  the plots of it match.
- `linewidth`: arc width in pixels.
- `speed`: initial playback rate, in simulated seconds per real second. `:auto` fits the whole
  trajectory into about two minutes.

# Notes
Positions go to the browser in the inertial frame and Cesium applies the rotation to fixed
itself, so the orbit holds still against the stars and the body turns beneath it. The epoch is
converted to UTC on the way out, whatever scale the history was recorded in, because Cesium
reads times as UTC.

An impulsive burn is marked with a point labelled with its ΔV, differenced from the velocities
either side of it rather than taken from whatever commanded it.

Earth is the only central body today.

# Returns
`nothing`. The view is drawn on the dashboard.

Throws an `ArgumentError` if the spacecraft has no propagated history, or if nothing in that
history is long enough to draw an arc.
"""
function orbitview(args...; kwargs...)
    title, rest = _split_panel(args)
    length(rest) == 1 || throw(ArgumentError(
        "orbitview takes one spacecraft, optionally after a panel name — " *
        "orbitview(sat), or orbitview(\"Mission\", sat). Got $(length(rest)) arguments."))
    set_scene!(title, build_orbitview(only(rest); kwargs...))
    ensure_browser()
    return nothing
end

"""
    orbitview!([panel], sat; kwargs...)

Add a spacecraft to a view that already exists, keeping what is there.

🔴 Not implemented. Several spacecraft in one view needs the CZML writer to merge entity sets,
which it does not yet do.

# Returns
Never returns. Throws an `ArgumentError` naming the way to get the same picture today — give the
second spacecraft its own panel with [`orbitview`](@ref). It refuses rather than quietly drawing
one spacecraft and dropping the other.
"""
function orbitview!(args...; kwargs...)
    throw(ArgumentError(
        "orbitview! is not implemented yet — a view holds one spacecraft. Give the second one " *
        "its own panel: orbitview(\"Chaser\", sat2)."))
end

"""
Build the CZML for one spacecraft: the packets, and the time span the clock runs over.
"""
function build_orbitview(sat; colors = nothing, linewidth::Real = 1.5, speed = :auto)
    epoch_iso, segments, maneuvers = _extract(sat; colors = colors, linewidth = linewidth)

    # Playback: fit the whole trajectory into about two minutes of real time unless told
    # otherwise, clamped so a short arc is not unwatchably slow nor a long one a blur.
    multiplier = if speed === :auto
        t_min = minimum(minimum(s.times_s) for s in segments)
        t_max = maximum(maximum(s.times_s) for s in segments)
        clamp((t_max - t_min) / 120.0, 1.0, 86_400.0)
    else
        Float64(speed)
    end

    return _build_czml(segments, maneuvers, epoch_iso;
                       id                = _entity_id(sat),
                       name              = sat.name,
                       marker_color      = (255, 255, 255, 255),
                       marker_pixel_size = 10,
                       multiplier        = multiplier,
                       inertial          = true)
end

# A CZML entity id has to survive being put in a URL and a JSON key.
_entity_id(sat) = replace(sat.name, r"[^A-Za-z0-9_\-]" => "_")
