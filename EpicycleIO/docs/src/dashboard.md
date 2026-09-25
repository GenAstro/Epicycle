```@meta
CurrentModule = EpicycleIO
```

# The Dashboard

Your plots appear in a browser page served from your own machine. The page stays open while you
keep working and updates as you plot, including from inside a running solve. A
[3D view](plotting_3d.md) is a panel on the same page, so a trajectory and the quantities
describing it sit side by side rather than in separate windows.

The first plot opens a tab. Later plots go to the same one.

## Naming Panels

A leading string names a panel. Without one you get the default panel, which is fine for a quick
look and awkward once you have more than one thing to see.

Plotting into a name that already exists replaces that panel. Run your script twice and you have
one figure, not two. Use `xyplot!` when you mean to add.

```julia
using EpicycleIO

t   = collect(0.0:0.25:12.0)
alt = 400.0 .+ 25.0 .* sin.(t)

xyplot("Altitude", t, alt; name = "Sat A")
xyplot!("Altitude", t, alt .+ 30; name = "Sat B")   # adds
xyplot("Altitude", t, alt .- 30; name = "Sat C")    # replaces both
```

## Titles, Axes and Limits

`panel!` sets the layout. Its keywords go to Plotly's layout, so Plotly's reference applies.
Nested settings are written with underscores rather than as a `Dict`: `yaxis_type` reaches
Plotly as `yaxis.type`, and `legend_orientation = "h"` is how you set a nested legend option.

Settings accumulate, so calling `panel!` twice for the same axis keeps both. Limits are an
ordered pair, so an axis runs backwards if you give them backwards.

```julia
using EpicycleIO

t   = collect(0.0:0.25:12.0)
alt = 400.0 .+ 25.0 .* sin.(t)

xyplot("Altitude", t, alt)

panel!("Altitude";
       title       = "Altitude above the ellipsoid",
       xaxis_title = "hours from epoch",
       yaxis_title = "km",
       yaxis_type  = "log")

xyplot("Decay", t, 500.0 .- 8.0 .* t)
panel!("Decay"; xaxis_range = [12, 0])     # runs backwards
```

Configuring a panel before plotting into it is allowed and creates it empty. If you name a panel
that does not exist and did not mean to create one, you get a warning listing the panels that do.

## Too Many Panels

Twelve panels in one tab is unusable. There are two ways out.

The grid and tabs control at the top right of the page switches between showing every panel at
once and showing one at a time with a strip of names.

Or open a tab holding only the panels you name:

```@raw html
<!-- doc-fragment -->
```
```julia
open_dashboard("Range residuals", "Range-rate residuals")   # navigation, its own window
open_dashboard("Porkchop")                                  # one plot, full size
open_dashboard()                                            # everything, as before
```

The panels are still one set. A tab showing a subset is a view of it, so a panel that two tabs
both show updates in both.

## Watching a Run

`xyplot!` is safe to call from inside a propagation or a solver callback. Writes are throttled and
happen on their own task, so publishing cannot slow the run that feeds it. Zoom into a panel
while this runs and the view stays where you put it, because updates redraw in place rather than
resetting what you set.

```julia
using EpicycleIO

# What a solver callback would do, one iteration at a time
for k in 1:30
    cost      = 10.0 / k
    violation = 1.0 / k^2
    xyplot!("Convergence", [Float64(k)], [cost];      mode = "markers", name = "cost")
    xyplot!("Feasibility", [Float64(k)], [violation]; mode = "markers", name = "max violation")
end
```

One `xyplot!` per point makes one series per call, which is fine for a few hundred and slows the
browser well before a few thousand. If you already hold the whole history, plot it once as a
growing series instead.

## Clearing and Closing

```julia
using EpicycleIO

t = collect(0.0:0.25:12.0)
xyplot("Altitude", t, 400.0 .+ 25.0 .* sin.(t))

clear!("Altitude")     # empty one panel, keeping the layout you set
clear_all!()           # remove every panel
EpicycleIO.close_dashboard()      # stop the server; the page stops updating
```

Closing the browser tab changes nothing in Julia. Your panels live in Julia rather than in the
page, so plotting again opens a new tab and picks up where the last one left off.

To stop plots opening tabs at all, in a script or a test:

```julia
using EpicycleIO

EpicycleIO.auto_open!(false)
```

## What Is Actually Happening

Julia writes each panel to a file and a small local server hands those files to the page, which
polls for changes. Nothing leaves your machine, and no account or key is needed anywhere.

The drawing libraries are the one exception. The page loads Plotly and Cesium from a public CDN
at exact pinned versions, each the first time it is needed. Plotly arrives when the page first
opens. Cesium arrives only when you first render a trajectory, so a session that plots data never
fetches it at all. Together they come to about 10 MB.

Your browser then caches both for a week. The cache is on disk rather than per session, so new
tabs, new Julia sessions, a browser restart and a reboot all reuse it. When the week is up your
browser asks whether the file has changed and is told it has not, which costs a few bytes instead
of another download, and because the versions are pinned that answer does not change.

So you need a network connection the first time you plot on a given browser, and after that you
do not. On a machine that cannot reach the internet at all, the panel stays blank. It will not
tell you why, so if a panel never draws on a restricted network, this is the first thing to
check.
