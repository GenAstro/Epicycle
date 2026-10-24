```@meta
CurrentModule = EpicycleIO
```

# Data plotting

EpicycleIO displays interactive data plots in a browser dashboard. Use it to
inspect time histories, compare results, explore distributions, and visualize
scalar fields. The plotting functions accept arrays, including data extracted
from a spacecraft's recorded history.

A plot contains one or more **traces**: a line, a set of points, a bar series, or
another representation of data. Traces appear in named **panels**. Supply trace
appearance settings to the plotting function and panel settings, such as axis
labels and legend position, to `panel!`.

Data plots use Plotly. Attribute names follow Plotly's naming, with underscores
for nested properties: `line_color` sets a line's color and `xaxis_title` sets
an axis title. You do not need to construct a Plotly figure to get started.

## Create your first plot

This example compares two synthetic altitude histories:

```julia
using EpicycleIO

hours = collect(0.0:0.1:6.0)
altitude_a = 400.0 .+ 20.0 .* sin.(hours)
altitude_b = 430.0 .+ 15.0 .* sin.(hours .+ 0.5)

xyplot("Altitude", hours, altitude_a;
       mode = "lines", name = "Spacecraft A", line_width = 2)
xyplot!("Altitude", hours, altitude_b;
        mode = "lines", name = "Spacecraft B", line_dash = "dash")
panel!("Altitude";
       title = "Altitude comparison",
       xaxis_title = "Hours from epoch",
       yaxis_title = "Altitude (km)",
       showlegend = true,
       legend_orientation = "h")
```

The first plotting call opens the dashboard automatically. Hover over a plot
to inspect values and use its browser controls to explore the data. See
[The dashboard](dashboard.md) for arranging panels and opening selected views.

The leading string identifies the panel. The `name` keyword labels a trace in
its legend. Multiple named traces can share one panel.

## Plot interfaces

Each convenience plotting function has a replacing form and an adding form:

- A call without `!` replaces the panel's plotted data. Reusing a panel name
  lets you rerun a script without creating another panel.
- A call with `!` adds traces while keeping existing traces. It does not append
  samples to an existing trace or change that trace's data.
- Omitting the panel name uses the default panel, `"Plot"`.
- Plotting functions return `nothing`; their output appears on the dashboard.

Replacing data resets the interactive view. Adding traces preserves the user's
zoom. Layout settings survive data replacement, and successive `panel!` calls
merge settings. You can configure a panel before plotting into it; check the
spelling if `panel!` warns that a panel does not exist yet.

### Accepted data shapes

The paired-series interface accepts the following forms:

| Call | Interpretation |
|:--|:--|
| `xyplot(y)` | Values against sample indices `1:length(y)`. |
| `xyplot(x, y)` | One series; both arrays must have equal lengths. |
| `xyplot(x, Y)` with a matrix | One series per column; rows must match `length(x)`. |
| `xyplot(x, samples)` with a vector of vectors | One trace per component when each outer entry is a sample. |
| `xyplot(x1, y1, x2, y2, ...)` | Independent pairs, each with its own length. |

A matrix supplied without `x` uses its row indices as horizontal coordinates.
A vector of vectors can also contain one full-length vector per series. If both
interpretations are possible, outer entries are treated as samples. Use a matrix
with one column per series to make the orientation explicit.

```julia
x = collect(0.0:0.25:4.0)
Y = hcat(sin.(x), cos.(x))
xyplot("Signals", x, Y;
       mode = "lines", name = ["Sine", "Cosine"],
       line_color = ["cyan", "orange"])
```

The paired-series rules also apply to `bar`, `scatterpolar`, and `scattergeo`.
Grid plots, `scatter3d`, and `band` have distinct interfaces described below.

### Units, frames, and time

Data plots use the supplied values. They do not infer units, transform reference
frames, or convert time systems. Convert data before plotting and include units
in axis labels.

Convert Epicycle `Time` objects before using them as coordinates. For example,
after propagating a spacecraft named `sat`, extract its history and express
elapsed time in hours:

```julia
using Epicycle

times, positions = history(Calc(epoch, sat),
                           Calc(position_vector, sat, EarthMJ2000Eq))
hours = [(t - first(times)) * 24 for t in times]
xyplot("Position", hours, positions;
       mode = "lines", name = ["x", "y", "z"])
panel!("Position";
       xaxis_title = "Hours from epoch",
       yaxis_title = "Position in EarthMJ2000Eq (km)")
```

The position history is a vector of component vectors, so this draws three
traces. For a calendar axis, supply date strings and set `xaxis_type = "date"`.
`NaN` and `Inf` in real-valued series are rendered as gaps.

## Supported plot types

Every function below accepts an optional leading panel name and trace keywords.
Each has an adding form, such as `heatmap!` or `band!`.

| Function | Data arguments | Use |
|:--|:--|:--|
| `xyplot` | `y`, `x, y`, or multiple pairs | Cartesian lines, points, and combinations. |
| `bar` | `categories, values` | Compare discrete categories. |
| `histogram` | `values` | Show a sample distribution. |
| `band` | `x, lower, upper` | Shade an interval or uncertainty envelope. |
| `contour` | `x, y, Z` or `Z` | Show levels of a scalar field. |
| `heatmap` | `x, y, Z` or `Z` | Show values as colored cells. |
| `scatterpolar` | `angle, radius` | Plot angular and radial data. |
| `scattergeo` | `longitude, latitude` | Plot geographic locations or tracks. |
| `scatter3d` | `x, y, z` | Plot points or lines in three dimensions. |
| `surface` | `x, y, Z` or `Z` | Draw a scalar field as a 3D surface. |

### Lines, markers, and filled curves

Use an explicit `mode` to choose the representation:

```julia
x = collect(0.0:0.2:4.0)
y = sin.(x)

xyplot("Line", x, y;
       mode = "lines", name = "Prediction",
       line_color = "cyan", line_width = 2, line_dash = "dash")
xyplot("Points", x, y;
       mode = "markers", name = "Measurements",
       marker_color = "orange", marker_size = 6, marker_symbol = "circle")
xyplot("Line and points", x, y;
       mode = "lines+markers", name = "Samples")
xyplot("Filled curve", x, abs.(y);
       mode = "lines", fill = "tozeroy",
       fillcolor = "rgba(0,220,255,0.2)", name = "Magnitude")
```

Use `line_shape = "hv"` for a step curve or `line_shape = "spline"` for a smooth
curve. Choose interpolation that represents your data appropriately. Error bars
can be specified with `error_y_type = "data"` and
`error_y_array = uncertainties`, where `uncertainties` contains one magnitude
per point.

`xyplot` and `xyplot!` use WebGL for every trace in the call when any series in
it exceeds 4,000 points. Rendering and supported attributes can differ from
ordinary scatter; other plot types do not use this automatic switch.

### Bar charts

```julia
bar("Maneuver budget", ["Departure", "Correction", "Arrival"],
    [2.82, 1.14, 0.48]; name = "Nominal", marker_color = "cyan")
bar!("Maneuver budget", ["Departure", "Correction", "Arrival"],
     [2.95, 1.20, 0.51]; name = "Reserve", marker_color = "orange")
panel!("Maneuver budget";
       barmode = "group", xaxis_title = "Maneuver", yaxis_title = "Delta-V (km/s)")
```

Use `barmode = "stack"` to stack series. For horizontal bars, supply numerical
values as the first coordinate and category labels as the second, with
`orientation = "h"`.

### Histograms

```julia
residuals = [0.1, -0.2, 0.05, 0.3, -0.1, 0.08, 0.15, -0.12]
histogram("Residual distribution", residuals;
          name = "Range", xbins_size = 0.1, marker_color = "cyan")
panel!("Residual distribution";
       xaxis_title = "Range residual (km)", yaxis_title = "Count")
```

Supply samples rather than precomputed bin counts. Plotly chooses bins
unless you configure them. `nbinsx` sets an upper bound on the number of bins;
`xbins_size` sets the width. Use `histnorm = "probability"` or
`"probability density"` for normalized displays and label the vertical axis
accordingly. For overlapping distributions, add traces with `histogram!`, set
panel `barmode = "overlay"`, and adjust each trace's `opacity`.

### Uncertainty bands

```julia
x = collect(0.0:0.1:4.0)
prediction = sin.(x)
sigma = fill(0.15, length(x))

band("Envelope", x, prediction .- 3 .* sigma, prediction .+ 3 .* sigma;
     name = "3 sigma", fillcolor = "rgba(120,120,120,0.25)")
xyplot!("Envelope", x, prediction;
        mode = "lines", name = "Prediction", line_color = "cyan")
panel!("Envelope"; xaxis_title = "Time (s)", yaxis_title = "Value")
```

All three arrays must have equal lengths. Draw the band first so overlaid data
appear on top. `band` constructs a filled polygon from the bounds; it is a
convenience helper rather than a separate Plotly trace type.

### Contours and heatmaps

For grid plots, rows of `Z` follow `y` and columns follow `x`:
`size(Z)` must equal `(length(y), length(x))`.

```julia
x = collect(range(-2.0, 2.0; length = 41))
y = collect(range(-1.0, 1.0; length = 21))
Z = [xi^2 + 2yi^2 for yi in y, xi in x]

contour("Cost contours", x, y, Z;
        colorscale = "Viridis", contours_showlabels = true,
        colorbar_title = "Cost")
heatmap("Cost map", x, y, Z;
        colorscale = "Viridis", colorbar_title = "Cost")
panel!("Cost contours"; xaxis_title = "x", yaxis_title = "y")
panel!("Cost map"; xaxis_title = "x", yaxis_title = "y")
```

Supplying only `Z` uses grid indices for the coordinates. Use `zmin` and `zmax`
to set color limits. The color bar explains scalar values; configure it with
trace attributes such as `colorbar_title`, separately from the panel legend.

### Polar plots

`scatterpolar` takes angles in degrees by default. Set `thetaunit = "radians"`
for angular inputs in radians. Radius values are used as supplied.

```julia
azimuth_deg = collect(0.0:15.0:180.0)
elevation_deg = 60.0 .* sin.(deg2rad.(azimuth_deg))
scatterpolar("Sky view", azimuth_deg, elevation_deg;
             mode = "lines+markers", name = "Pass")
panel!("Sky view";
       polar_angularaxis_direction = "clockwise",
       polar_angularaxis_rotation = 90,
       polar_radialaxis_range = [90, 0])
```

This convention places north at the top and zenith at the center. Convert
elevation to degrees yourself: `thetaunit` affects only the angular coordinate.
Polar layout settings begin with `polar_`.

### Geographic plots

`scattergeo` takes longitude in degrees east and latitude in degrees north.
Use `mode = "lines"` for tracks or `"markers"` for locations:

```julia
scattergeo("Stations", [-116.9, -4.2, 148.98], [35.4, 40.4, -35.4];
           mode = "markers", name = "Stations", marker_size = 8)
panel!("Stations";
       geo_projection_type = "equirectangular", geo_showland = true)
```

Geographic layout settings begin with `geo_`. Convert radians before plotting
and use geodetic latitude when comparing positions with map coastlines. The
plot does not perform coordinate-frame conversions.

### Three-dimensional points, lines, and surfaces

`scatter3d` takes three equal-length coordinate arrays. Its `mode` can be
`"markers"`, `"lines"`, or `"lines+markers"`. `surface` uses the same grid
orientation as `contour` and `heatmap`.

```julia
angle = collect(range(0.0, 4pi; length = 200))
scatter3d("Helix", cos.(angle), sin.(angle), angle ./ (2pi);
          mode = "lines", name = "Curve")
panel!("Helix";
       scene_xaxis_title = "x", scene_yaxis_title = "y",
       scene_zaxis_title = "Turns", scene_aspectmode = "data")

x = collect(range(-2.0, 2.0; length = 41))
y = collect(range(-1.0, 1.0; length = 21))
Z = [xi^2 + 2yi^2 for yi in y, xi in x]
surface("Cost surface", x, y, Z;
        colorscale = "Viridis", colorbar_title = "Cost")
panel!("Cost surface";
       scene_xaxis_title = "x", scene_yaxis_title = "y",
       scene_zaxis_title = "Cost")
```

Three-dimensional axis settings begin with `scene_`. For coordinates with the
same units, `scene_aspectmode = "data"` preserves relative scale. Use
`scatter3d!` to overlay points or a curve on a surface.

For animated spacecraft trajectories around Earth, use `orbitview`. See
[3D plotting](plotting_3d.md) for that interface and its playback controls.

## Legends and trace appearance

Set legend text with trace `name`. For multi-component data, a scalar name such
as `"Position"` produces `Position[1]`, `Position[2]`, and so on. Supply a vector
of names to label components individually:

```julia
x = collect(0.0:0.25:4.0)
Y = hcat(sin.(x), cos.(x))
xyplot("Signals", x, Y;
       mode = "lines", name = ["Sine", "Cosine"],
       line_color = ["cyan", "orange"], line_dash = ["solid", "dash"])
panel!("Signals";
       showlegend = true, legend_title_text = "Signal",
       legend_orientation = "h", legend_x = 0, legend_y = 1.15)
```

Panel `showlegend = true` explicitly displays a legend, including for a single
trace. Set `showlegend = false` on a trace to omit its entry, or on the panel
to hide the whole legend. Give related traces the same `legendgroup` string
to group their legend interactions.

When a paired-series call expands into multiple traces:

- A scalar style value applies to every trace.
- A vector supplies one value per trace for attributes that do not accept
  per-point arrays, such as `line_color`, `line_width`, and `line_dash`.
- Attributes that accept per-point arrays, such as `marker_size` and
  `marker_color`, retain that meaning. Nest arrays to supply separate values
  to separate traces: `marker_size = [sizes_a, sizes_b]`.

For example, `marker_size = [[3], [6]]` supplies a different marker-size array
to each of two traces. The distribution rule does not apply to the single-trace
grid, spatial, and band interfaces, or to directly constructed PlotlyBase traces.

## Axis labels, ranges, and scales

Use `panel!` for layout settings and include units in labels:

```julia
panel!("Altitude";
       xaxis_title = "Hours from epoch",
       yaxis_title = "Altitude (km)",
       xaxis_range = [0, 6], yaxis_range = [350, 500],
       xaxis_showgrid = true, yaxis_showgrid = true)
```

Ranges are ordered pairs; `[6, 0]` reverses an axis. Set
`yaxis_autorange = true` to restore automatic vertical limits.
For a logarithmic Cartesian axis, use `yaxis_type = "log"`. Explicit range
endpoints then use base-10 exponents: `[0, 3]` displays values from 1 to 1,000.
Only positive values can appear on a log axis.

Set `xaxis_tickvals` and `xaxis_ticktext` to control tick positions and displayed
labels. Choose layout settings for the coordinate system you are using:

| Coordinates | Layout prefixes |
|:--|:--|
| Cartesian | `xaxis_`, `yaxis_` |
| Polar | `polar_angularaxis_`, `polar_radialaxis_` |
| Geographic | `geo_` |
| Three-dimensional data | `scene_xaxis_`, `scene_yaxis_`, `scene_zaxis_` |

Use underscore keywords for nested properties, such as
`legend_orientation = "h"`, rather than nested string-keyed dictionaries.

## Combine plots and update results

Traces with compatible coordinates can share a panel, such as a band and a
Cartesian line, or a surface and 3D scatter. Keep Cartesian, polar, geographic,
and three-dimensional plots in separate panels. Use separate panels for orbit
views as well.

```julia
open_dashboard("Altitude", "Signals")
clear!("Signals")   # Empty the panel, retaining its layout.
```

`clear_all!()` removes all panels and their layouts from the dashboard.

Plotting calls can be made during propagation or solver callbacks. Dashboard
publishing is throttled and runs on a separate task. Adding calls create more
traces; use the replacing form with updated arrays to refresh a complete history.
Avoid creating a new trace for every sample in a long run.

## Use PlotlyBase traces directly

For plot types without convenience functions, construct PlotlyBase traces and
pass them to `xyplot` or `xyplot!`. EpicycleIO exports `PlotlyBase`, so it is
available after `using EpicycleIO`:

```julia
using EpicycleIO

samples = [398.0, 401.0, 399.5, 402.0, 400.2]
xyplot("Insertion altitude",
       PlotlyBase.box(y = samples, name = "Monte Carlo", boxmean = "sd"))
panel!("Insertion altitude"; yaxis_title = "Altitude (km)", showlegend = true)

xyplot("Residuals",
       PlotlyBase.violin(y = [-0.2, 0.1, 0.05, -0.1, 0.3],
                         name = "Range", box_visible = true))
```

Supply several trace objects in one call, or combine directly constructed traces
with convenience-function traces in compatible panels. Array expansion and
attribute distribution do not apply to trace objects; configure each trace when
constructing it.

You can also supply a layout object:

```julia
panel!("Insertion altitude",
       PlotlyBase.Layout(yaxis_title = "Altitude (km)", showlegend = true))
```

Additional styles are available through trace constructors supported by the
installed PlotlyBase and dashboard renderer. The panel interface distinguishes
Cartesian, polar, geographic, and three-dimensional families; specialized Plotly
subplot types may require configuration beyond the convenience interface.

## Troubleshooting

| Symptom | What to check |
|:--|:--|
| Labels or legend settings do not change | Supply layout settings to `panel!` using the plot's panel name. |
| A single trace has no legend | Provide a trace `name` and set panel `showlegend = true`. |
| A call removes earlier data | Use the adding form with `!` to keep existing traces. |
| A grid is rejected or transposed | Rows must follow `y` and columns must follow `x`; check `size(Z)`. |
| Polar or geographic data look incorrect | Check degrees versus radians and the coordinate-system layout settings. |
| Traces cannot share a panel | Use separate panels for different coordinate families. |
| An attribute-name warning appears | Check spelling and whether the property belongs to a trace or the layout. The bundled schema can lag newer renderer attributes, so a warning does not always mean a property is unsupported. |

See the [API reference](api.md) for signatures and the repository's `examples`
directory for complete plotting scripts.
