```@meta
CurrentModule = EpicycleIO
```

```@setup io
# The doc build is headless. This stops each block opening a browser tab. A script or a test
# does the same.
using EpicycleIO
EpicycleIO.auto_open!(false)
```

# EpicycleIO

The EpicycleIO module provides plotting, reporting, and three-dimensional trajectory
visualization for Epicycle results. Plots are drawn by Plotly.js and trajectories by Cesium, both
in a browser page served from the local machine. The page updates while a script is still
running, so a propagation or a solver iteration can be watched as it proceeds.

The plotting and reporting functions take arrays, and have no knowledge of spacecraft, states
or coordinate systems. `history` in AstroCallbacks produces the arrays. Data from two spacecraft
needs no special handling: the two recordings are two pairs of arrays with their own lengths.

The three-dimensional view is the exception. `orbitview` takes a Spacecraft and reads the
reference frame and the epoch from it. Earth is the only central body supported for these views,
and the body comes from the frame the history was recorded in.

Attribute names are Plotly's own, such as `name`, `mode`, `line_color` and `marker_size`, and
values pass through unchanged. Plotly's reference documentation applies directly to this
interface.

References:

- Plotly JavaScript reference, <https://plotly.com/javascript/reference/>
- CesiumJS, <https://cesium.com/platform/cesiumjs/>
- CZML specification, <https://github.com/AnalyticalGraphicsInc/czml-writer/wiki>

## Installation

EpicycleIO is registered in the Gen Astro registry rather than Julia's General registry, because it
is released under the Gen Astro Source Available License. Add the registry once, then install
as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("EpicycleIO")
```

General is still required, since these packages depend on packages registered there.

## Quick Start

The Epicycle umbrella does not re-export `xyplot`; the name collides with Plots.jl and Makie. Load
EpicycleIO explicitly. The example below propagates a spacecraft, draws the trajectory, plots two
quantities describing it, and writes the same data to a file. It runs as written.

```@example io
using Epicycle
using EpicycleIO

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 name  = "Explorer")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))

# The trajectory. orbitview takes the spacecraft; the frame and epoch come from it
orbitview("Trajectory", sat)

# Columns out of the recorded history, one walk for all of them
t, r, rmag = history(Calc(epoch,              sat),
                     Calc(position_vector,    sat, EarthMJ2000Eq),
                     Calc(position_magnitude, sat))

# Time is not a number. Subtract the epoch before plotting against it
hours = [(x - first(t)) * 24 for x in t]

# The leading string names a panel; a browser tab opens on the first plot
xyplot("Radius", hours, rmag; name = "|r|", line_color = "cyan")

# Nested Plotly settings are spelled with underscores, so yaxis_title becomes yaxis.title
panel!("Radius"; xaxis_title = "hours from epoch", yaxis_title = "km")

# A position column is a vector of 3-element vectors and draws as three series
xyplot("Position", hours, r; name = ["x", "y", "z"])

# The same arrays as a file. The keywords become the column headings
report("flight.txt"; time = [x.jd for x in t], position = r, radius = rmag)
```

The API is documented with docstrings; Plotly's own documentation is not duplicated here. In the
REPL, type `?` to enter help mode, then enter a name. `?xyplot` gives the shape rules and how
attributes are distributed across several series, and `?orbitview` gives what is read from a
spacecraft. `names(EpicycleIO)` lists every exported name.

## Data Plotting

`xyplot` draws lines and markers and replaces whatever the panel held, and `xyplot!` adds to a
panel. Other functions draw other Plotly traces and take the name of the trace they draw. Every
one has a `!` form.

Array shapes follow the conventions used by MATLAB and matplotlib. A single array plots against
its index, and a matrix draws one series per column. Several x/y pairs may be given in one
call, and lengths are checked pairwise.

```@example io
using EpicycleIO

t   = collect(0.0:0.25:12.0)
alt = 400.0 .+ 25.0 .* sin.(t)
r   = [[7000.0 + 10i, 100.0i, 1300.0] for i in t]

# Several forms of the same call
xyplot("A", alt)                            # index on the horizontal axis
xyplot("B", t, alt)                         # one series
xyplot("C", t, hcat(alt, alt .+ 10))        # one series per column of a matrix
xyplot("D", t, alt, t[1:20], alt[1:20])     # two series, independent lengths

# A scalar attribute applies to every series; a vector gives one value per series
xyplot("E", t, r; name       = ["radial", "along-track", "cross-track"],
                line_color = ["red", "green", "blue"],
                line_width = [1, 2, 3])

# Other trace types take their data positionally and pass attributes through
histogram("Residuals", randn(200); nbinsx = 40)

# Any Plotly trace can be built directly and handed over, whether or not a wrapper exists
xyplot("Distribution", PlotlyBase.violin(y = randn(200), box_visible = true, name = "range"))
```

A `NaN` or an `Inf` draws as a break in the line. Gaps are ordinary in this domain: a quantity
between ground station passes, argument of periapsis on a circular orbit, a solver step that did
not converge.

## Three-Dimensional Views

`orbitview` draws a propagated trajectory about the central body, with a clock that can be
played, paused and run at a chosen rate. The view is a panel on the same page as the plots, so
a trajectory and the quantities describing it are seen together.

Everything the view needs comes from the spacecraft: the trajectory from its recorded history,
the epoch from the earliest sample, one coloured arc per propagation segment, and the label from
its name. Each impulsive maneuver is marked with a point labelled with its magnitude, computed
from the velocity change across the maneuver.

```@raw html
<!-- doc-fragment -->
```
```julia
# Takes the spacecraft, not arrays. See the Quick Start above for the setup
orbitview("Trajectory", sat)

# Playback defaults to fitting the trajectory into about two minutes of real time
orbitview("Trajectory", sat; speed = 3600)     # one hour of flight per second
```

The camera is fixed in the inertial frame, so the orbit holds still against the stars and the
central body rotates beneath it. A view holds one spacecraft, and a second spacecraft takes its
own panel.

## Reporting

`report` writes named columns to a delimited text file. The keywords become the column
headings, and a vector-valued column expands into one column per component. Columns must be the
same length, so columns from two different recordings cannot share a report.

```@example io
using EpicycleIO

t = collect(0.0:0.25:6.0)
r = [[7000.0 + 10i, 100.0i, 1300.0] for i in t]
v = [[0.0, 7.35 - 0.001i, 1.0] for i in t]

report("flight.txt"; time = t, position = r, velocity = v)
```

`report` does not interpret epochs, since a bare Julian date does not record whether it is TT
or UTC. A `Time` passed directly raises an error, naming `t.jd`, `t.mjd` and `t.isot` as the
conversions available.

## The Dashboard

A leading string names a panel, and plotting into a name that already exists replaces that
panel, so a script run twice leaves one figure rather than two. A control at the top right of the
page switches between showing every panel at once and showing one panel at a time.

```@example io
using EpicycleIO

t = collect(0.0:0.25:12.0)
xyplot("Altitude", t, 400.0 .+ 25.0 .* sin.(t))
xyplot("Radius",   t, 7000.0 .+ 10.0 .* t)

open_dashboard("Altitude")    # a tab showing named panels only
clear!("Altitude")            # empty one panel, leaving its layout
clear_all!()                  # remove every panel
EpicycleIO.close_dashboard()             # stop the server; the page stops updating

EpicycleIO.auto_open!(false)  # stop plots opening tabs, for a script or a test
```

Plotly and Cesium are downloaded the first time they are needed and then cached by the
browser, so the first plot on a machine with no network connection will not draw. Nothing else
leaves the local machine, and no account or key is required.

## Reference

The names a reader works through most. Every exported name is on the [API](api.md) page.

```@docs
xyplot
xyplot!
panel!
orbitview
report
open_dashboard
```

## Where to go next

- [Data plotting](data_plotting.md), what can be drawn, styling, and grids
- [3D plotting](plotting_3d.md), trajectories about a central body, with a clock
- [The dashboard](dashboard.md), panels, live update, and keeping plots where they are wanted
- [Reporting](reporting.md), columns, headings, and what a value may be
- [API](api.md), every exported name
