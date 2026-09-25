```@meta
CurrentModule = EpicycleIO
```

# Data Plotting

Data plotting covers everything drawn against axes: time histories, residuals, ground tracks,
distributions, and scalar fields over a grid.

The functions on this page operate on arrays and have no knowledge of where the arrays came
from. The examples below therefore build their own data, and they run as written. In practice
the arrays usually come from `history` in AstroCallbacks, but nothing on this page requires it.

Attribute names are Plotly's own and values are passed through unchanged, so
[Plotly's reference](https://plotly.com/javascript/reference/) applies directly. This page
describes what EpicycleIO adds: how arrays become series, and how attributes are spread across
them.

## Series and Shapes

`xyplot` replaces the contents of a panel and `xyplot!` adds to a panel. Array shapes follow the
conventions used by MATLAB and matplotlib: a single array plots against its index, a matrix
draws one series per column, and several x/y pairs may be given in one call. Lengths are checked
in pairs, so the pairs need have nothing to do with one another.

```julia
using EpicycleIO

t   = collect(0.0:0.25:12.0)
alt = 400.0 .+ 25.0 .* sin.(t)

# The leading string names a panel. Without one, the default panel is used
xyplot("Altitude", t, alt; name = "altitude", mode = "lines", line_color = "cyan", line_width = 2)

# plot! adds to a panel rather than replacing it. Lengths need not match
t2   = collect(0.0:0.5:9.0)
alt2 = 410.0 .+ 15.0 .* cos.(t2)
xyplot!("Altitude", t2, alt2; name = "second spacecraft")

# One array plots against its index
xyplot("Index", alt)

# A matrix draws one series per column
Y = hcat(alt, alt .+ 10, alt .- 10)
xyplot("Columns", t, Y)

# Several pairs in one call, each with its own length
xyplot("Pairs", t, alt, t2, alt2)
```

A position column arrives from `history` as a vector of 3-element vectors rather than a matrix.
That shape draws as three series without conversion, which is the case the shape handling exists
for.

```julia
using EpicycleIO

t = collect(0.0:0.25:12.0)

# The shape history returns for a vector-valued column: one entry per sample,
# each holding that sample's components
r = [[7000.0 + 10i, 100.0i, 1300.0] for i in t]

# Three series, named position[1], position[2] and position[3]
xyplot("Position", t, r; name = "position")

# Name the components when they mean something
xyplot("Position", t, r; name = ["radial", "along-track", "cross-track"])
```

## Attributes Across Several Series

A scalar attribute applies to every series. A vector gives one value per series.

Plotly reads an array on some attributes as one value per data point rather than per series, and
marks those attributes `arrayOk` in its schema. EpicycleIO reads that schema rather than
maintaining a list, so the rule is correct for every attribute. On an `arrayOk` attribute,
Plotly's per-point meaning is preserved and per-series is expressed by nesting.

```julia
using EpicycleIO

t = collect(0.0:0.25:12.0)
r = [[7000.0 + 10i, 100.0i, 1300.0] for i in t]

# One value per series
xyplot("Styled", t, r; name       = ["x", "y", "z"],
                     line_color = ["red", "green", "blue"],
                     line_width = [1, 2, 3])

# marker_size is arrayOk, so a flat vector keeps Plotly's per-point meaning
xyplot("Per point", t, r; mode = "markers", marker_size = [4, 4, 4])

# Nesting always means per series, whatever the schema says
xyplot("Per series", t, r; mode = "markers", marker_size = [[3], [6], [9]])
```

A vector whose length does not match the series count raises an error rather than recycling.

`name` is the only attribute that rewrites a value rather than placing it. A scalar name across
three series becomes `position[1]`, `position[2]` and `position[3]`, because three identical
legend entries are never what was intended.

## Gaps

A `NaN` or an `Inf` draws as a break in the line rather than raising an error. Gaps are ordinary
in this domain: a quantity between ground station passes, argument of periapsis on a circular
orbit, or a solver step that did not converge. None of those should cost the user the plot.

```julia
using EpicycleIO

t = collect(0.0:0.5:20.0)
el = [x < 6 || x > 14 ? NaN : 30.0 + 10.0 * sin(x) for x in t]

# The pass is drawn; the periods with no visibility are gaps rather than zeros
xyplot("Elevation", t, el; mode = "lines", name = "station 1")
```

## Other Trace Types

Each of the following draws the Plotly trace it is named after and has a `!` form that adds to a
panel rather than replacing it.

| Function | Draws |
|---|---|
| `scatterpolar` | angle and radius, such as a ground station pass or an antenna pattern |
| `contour`, `heatmap` | a scalar field over a grid, such as a porkchop plot or a coverage map |
| `histogram` | how a set of values is distributed |
| `bar` | a value per category, such as ΔV by maneuver |
| `band` | a filled region between two bounds, such as a covariance envelope |
| `scattergeo` | longitude and latitude on a map |
| `scatter3d`, `surface` | three axes of data. For a trajectory about a body see [3D plotting](plotting_3d.md) |

The table is a snapshot. To list the functions the package actually exports, run
`names(EpicycleIO)`.

### Grids

`contour`, `heatmap` and `surface` take two axes and a matrix. The matrix is indexed
`[row, column]`, so rows follow `y` and columns follow `x`. A transposed grid is refused rather
than drawn sideways.

```julia
using EpicycleIO

departure = collect(0.0:2.0:60.0)      # columns
arrival   = collect(100.0:2.0:200.0)   # rows

# C3 is indexed [row, column], so size(C3) == (length(arrival), length(departure))
C3 = [12.0 + 0.02 * (d - 30)^2 + 0.01 * (a - 150)^2 for a in arrival, d in departure]

contour("Porkchop", departure, arrival, C3;
        colorscale = "Viridis", colorbar_title = "C3 km²/s²")

panel!("Porkchop"; xaxis_title = "departure, days from epoch",
                   yaxis_title = "arrival, days from epoch")
```

### Bands

A band is drawn before the data inside it so the points land on top rather than underneath.

```julia
using EpicycleIO

t     = collect(0.0:0.5:24.0)
sigma = 0.05 .+ 0.002 .* t
resid = [0.03 * sin(x) + 0.01 * randn() for x in t]

band("Range residuals", t, -3 .* sigma, 3 .* sigma;
     name = "3 sigma", fillcolor = "rgba(120,120,120,0.25)", line_width = 0)

xyplot!("Range residuals", t, resid; mode = "markers", name = "residual", marker_size = 4)

panel!("Range residuals"; xaxis_title = "hours", yaxis_title = "km")
```

### Polar Plots

A polar plot is drawn the way an operator expects it: north at the top, bearings increasing
clockwise, and the zenith at the centre.

Two of the settings below are correctness rather than preference. Mathematics places zero to the
right and increases counterclockwise, while a compass bearing places zero at the top and
increases clockwise, so an azimuth drawn in the mathematical convention is mirrored rather than
merely rotated. Elevation runs from 90 degrees at the centre to 0 at the rim.

```julia
using EpicycleIO

az = collect(range(0.4, 2.6; length = 40))                 # radians
el = [deg2rad(75.0) * exp(-((x - 1.5) / 0.7)^2) for x in az]

# thetaunit converts the angle. There is no equivalent for the radial axis, so an
# elevation in radians against a degree range would place every point on the rim
scatterpolar("Sky view", az, rad2deg.(el);
             thetaunit = "radians", mode = "lines", name = "Pass 3")

panel!("Sky view";
       polar_angularaxis_direction = "clockwise",   # compass sense
       polar_angularaxis_rotation  = 90,            # north at the top
       polar_radialaxis_range      = [90, 0])       # zenith at the centre
```

## Traces Without a Wrapper

Plotly has around forty trace types and all of them are reachable. A trace built with PlotlyBase
can be handed to `xyplot` directly, whether or not a wrapper exists for it. A wrapper is a
convenience for a common case and never a gate.

```julia
using EpicycleIO

residuals = randn(200)

xyplot("Distribution", PlotlyBase.violin(y = residuals, box_visible = true, name = "range"))
```

## Panel Coordinate Families

A panel has a coordinate family because Plotly does. Polar, geographic and three-dimensional
traces are drawn on their own subplot types and cannot share a panel with a cartesian trace.
Adding a trace from a different family to an existing panel raises an error naming both families
rather than silently drawing nothing. This is Plotly's constraint reported rather than a limit
imposed here.
