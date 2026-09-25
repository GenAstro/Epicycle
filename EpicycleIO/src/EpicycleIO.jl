# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    EpicycleIO

The EpicycleIO module provides plotting, reporting, and three-dimensional trajectory
visualization for Epicycle results. Plots are drawn by Plotly.js and trajectories by Cesium,
both running in a browser page served from the local machine.

Attribute names are Plotly's own (e.g., `name`, `mode`, `line_color`, `marker_size`) and values
are passed through unchanged, so Plotly's reference documentation applies directly.

The plotting and reporting functions operate on arrays and have no knowledge of spacecraft,
states, or coordinate systems. Arrays are typically produced by `history` in AstroCallbacks. The
three-dimensional view is the exception. `orbitview` accepts a Spacecraft because a trajectory
drawn about a central body requires a reference frame and an epoch.

`xyplot` is named for what it draws, like every other wrapper here. That also keeps it clear of
Plots.jl and Makie, which both claim `plot`; the other names still collide, so a user who wants one
of those alongside this renames on import: `using Plots: plot as pplot`.

```julia
using Epicycle
using EpicycleIO

t, alt = history(Calc(epoch, sat), Calc(altitude, sat))
xyplot(t, alt; name = "altitude")
```

Attribute names, trace types and layout settings are documented by Plotly:
<https://plotly.com/javascript/reference/>
"""
module EpicycleIO

using PlotlyBase
using JSON
using Scratch: @get_scratch!
using Sockets

# The 3D view is the one part of this package that knows what a spacecraft is, and it has to:
# a trajectory on a globe means nothing without a frame and an epoch. Confined to src/gl, so a
# change here cannot reach plotting and the plotting tests need no spacecraft.
using AstroModels: Spacecraft, HistorySegment
using AstroEpochs: Time
using Dates: DateTime, Millisecond

# Our wrappers are named after the Plotly traces they draw, which is the point — the name is
# already the documented one. That also means PlotlyBase exports every one of them, so defining
# rival functions would make each name ambiguous the moment a user loaded both, which the
# escape hatch in §6 asks them to do.
#
# So these are the same functions, with methods added. PlotlyBase's build a trace from keywords
# or a Dict; ours take a panel title or arrays and draw. One binding, nothing to disambiguate.
import PlotlyBase: bar, contour, heatmap, histogram, scatter3d, scattergeo, scatterpolar,
                   surface

const TRACE_WRAPPERS = (:scatterpolar, :bar, :histogram, :contour, :heatmap, :scattergeo,
                        :scatter3d, :surface, :band)

# PlotlyBase itself, so the escape hatch works. Plotly has around forty trace types and only
# nine have wrappers here; the rest are reached by building the trace with PlotlyBase and handing
# it to `xyplot`. A user who installs EpicycleIO gets PlotlyBase as an indirect dependency, which
# Julia will not let them load by name, so without this export the documented escape hatch fails
# for every reader who tries it.
export PlotlyBase

export xyplot, xyplot!, panel!, clear!, clear_all!
export scatterpolar, bar, histogram, contour, heatmap, scattergeo, scatter3d, surface, band
export scatterpolar!, bar!, histogram!, contour!, heatmap!, scattergeo!, scatter3d!,
       surface!, band!
export report
export orbitview, orbitview!
export open_dashboard
# Reached rarely enough to be qualified. An export is a name in every user's namespace and, once
# registered, a semver contract, so the bar is whether a user calls it often rather than whether it
# is public. `EpicycleIO.narrate(...)` is public and unexported, which is the usual place for a
# thing a demo drives and an analysis does not.
#     narrate, clear_narration, close_dashboard

include("schema.jl")
include("shapes.jl")
include("panels.jl")
include("server.jl")
include("publish.jl")
include("api.jl")
include("report.jl")
include("narration.jl")
include("gl/palette.jl")
include("gl/czml.jl")
include("gl/segments.jl")
include("gl/orbitview.jl")

# A plot published in the last fraction of a second must not be the one that gets dropped
# because the process ended.
function __init__()
    atexit(flush_publishes)
    return nothing
end

end # module
