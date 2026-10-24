```@meta
CurrentModule = EpicycleIO
```

# 3D Plotting

EpicycleIO displays animated spacecraft trajectories and three-dimensional data
plots in the browser dashboard. Choose the function that matches your input:

| Function | Use | Input |
|:--|:--|:--|
| `orbitview` | An animated trajectory around Earth | A spacecraft with propagated history |
| `scatter3d` | Points or lines in three data axes | Three equal-length coordinate arrays |
| `surface` | A scalar field drawn as a surface | A matrix, optionally with coordinate arrays |

The first plotting call opens the dashboard automatically. See
[The dashboard](dashboard.md) for panel layout and browser requirements.

## Display a spacecraft trajectory

Use `orbitview` after propagating a spacecraft. The view reads positions and
times from its recorded history and uses the spacecraft's name as its label.
There is no need to extract position arrays or supply a separate epoch.

This example propagates a low Earth orbit for six hours:

```julia
using Epicycle

sat = Spacecraft(
    state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
    time = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
    name = "Explorer",
)

prop = OrbitPropagator(
    ForceModel(PointMassGravity(earth, ())),
    IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0),
)

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))
orbitview("Trajectory", sat)
```

The Cartesian state contains position in kilometers followed by velocity in
kilometers per second. Each recorded propagation segment appears as a separate
colored arc, with the spacecraft animated along the trajectory.

A leading string selects the panel. Without it, `orbitview(sat)` uses the default
panel, `"Plot"`. Calling `orbitview` again with the same panel name replaces its
contents and resets the displayed time and camera framing. Call it again after
further propagation to display the updated history.

### Frames and time

The camera follows an inertial frame, so Earth rotates beneath the trajectory
during playback. Recorded times are converted to UTC for display. Histories
containing different time scales are supported.

Earth is currently the only supported central body. The viewer treats recorded
positions as Earth-centered inertial coordinates and converts kilometers to
meters for the browser. It does not transform positions from the history's
coordinate system. Use an Earth-centered inertial history; positions recorded in
a body-fixed frame or another coordinate system require conversion before viewing.

### Playback controls

The controls beneath each view let you:

- **Play / pause** the animation.
- **Reset** the time to the start of the recorded trajectory.
- Run the animation **slower** or **faster**.
- Toggle the **stars** background.

The time readout shows the current UTC epoch. Each trajectory panel has its own
clock and playback controls.

By default, `speed = :auto` chooses a rate that plays the recorded time span in
about two minutes. To set an initial rate, supply the number of simulated seconds
per real second:

```julia
orbitview("Trajectory", sat; speed = 3600)  # One simulated hour per real second.
```

Playback rates in the browser are limited to between 1 and 86,400 simulated
seconds per real second. Automatic playback uses the same limits, so very short
or very long trajectories may take a different amount of time to play.

### Arc colors and width

The default arc colors use the same palette as the data plots. Use `colors` to
supply a different palette and `linewidth` to set the arc width in pixels:

```julia
orbitview("Trajectory", sat;
          colors = [(55, 200, 255, 255), (255, 165, 0, 255)],
          linewidth = 2.5,
          speed = 600)
```

Each color is an `(r, g, b, a)` tuple of integers from 0 to 255; an alpha value
of 255 is opaque. Colors follow propagation-segment order and cycle when there
are more arcs than colors. A single tuple, such as
`colors = (55, 200, 255, 255)`, applies one color to every arc.
The default line width is 1.5 pixels.

### Maneuver markers

When the history contains a single-sample impulsive-maneuver marker between
propagation arcs, the view marks the burn and labels it with its name and ΔV
magnitude in kilometers per second. The ΔV comes from the recorded velocity
change across the marker; no separate maneuver list is needed.

A velocity discontinuity without a recorded marker does not produce a burn
label. Markers without a preceding propagation arc, or with a velocity change
below 1 mm/s, are omitted.

## Show a trajectory beside its data

Use separate named panels to display an orbit view beside plots from the same
run. Continuing the spacecraft example above:

```julia
t, r = history(Calc(epoch, sat), Calc(position_vector, sat, EarthMJ2000Eq))
hours = [(epoch_i - first(t)) * 24 for epoch_i in t]
radius_km = [sqrt(sum(abs2, position_i)) for position_i in r]

xyplot("Radius", hours, radius_km; name = "|r|")
panel!("Radius";
       xaxis_title = "Hours from epoch",
       yaxis_title = "Radius (km)")

open_dashboard("Trajectory", "Radius")
```

The radius plot provides a quantitative view of the same history. Its horizontal
axis does not move with the orbit view's playback clock.

Each orbit view currently holds one spacecraft. `orbitview!` is not implemented
and raises an `ArgumentError`. To compare spacecraft, give each its own panel:

```julia
# sat1 and sat2 must each have propagated history.
orbitview("Chaser", sat1)
orbitview("Target", sat2)
open_dashboard("Chaser", "Target")
```

## Plot points and lines in three axes

Use `scatter3d` for a three-dimensional data plot. It takes one `x`, `y`, and `z`
value per point and uses the supplied values directly, without interpreting
units, reference frames, or epochs.

```julia
using EpicycleIO

angle = collect(range(0.0, 4π; length = 200))
x = cos.(angle)
y = sin.(angle)
z = angle ./ (2π)

scatter3d("Helix", x, y, z;
          mode = "lines", name = "Curve", line_color = "cyan")
panel!("Helix";
       scene_xaxis_title = "x",
       scene_yaxis_title = "y",
       scene_zaxis_title = "Turns",
       scene_aspectmode = "data")

scatter3d!("Helix", [first(x)], [first(y)], [first(z)];
           mode = "markers", name = "Start", marker_size = 5)
```

Choose `mode = "markers"`, `"lines"`, or `"lines+markers"` for the display.
`scatter3d` replaces a panel's plotted data; `scatter3d!` adds a trace. Axis layout
properties begin with `scene_`. For coordinates with the same units,
`scene_aspectmode = "data"` preserves their relative scale.

## Plot a surface

Use `surface` for values on a rectangular grid. Rows of `Z` correspond to `y`,
and columns correspond to `x`, so `size(Z)` must equal
`(length(y), length(x))`.

```julia
using EpicycleIO

x = collect(range(-2.0, 2.0; length = 41))
y = collect(range(-1.0, 1.0; length = 21))
Z = [xi^2 + 2yi^2 for yi in y, xi in x]

surface("Cost", x, y, Z; colorscale = "Viridis", colorbar_title = "Cost")
panel!("Cost";
       scene_xaxis_title = "x",
       scene_yaxis_title = "y",
       scene_zaxis_title = "Cost")
```

With `surface("Cost", Z)`, the coordinate axes use the matrix indices.
`surface!` adds a surface to an existing 3D data panel.

Keep orbit views and data plots in separate panels. An orbit view uses a globe
renderer, while `scatter3d` and `surface` use Plotly's three-axis scenes. Plotly
layout settings supplied through `panel!` do not configure the orbit viewer.

## Troubleshooting

| Symptom | What to check |
|:--|:--|
| `orbitview` reports missing history | Call `propagate!` before displaying the spacecraft. An initial state alone is insufficient. |
| The history has no drawable arcs | At least one propagation segment must contain two or more samples. Empty segments and single-sample markers cannot form an arc. |
| The trajectory has an unexpected orientation or location | Check that recorded positions are Earth-centered inertial coordinates in kilometers. The view does not transform other frames. |
| A maneuver has no label | Check that the history records a single-sample marker and a measurable velocity change. |
| The globe panel remains blank | Check that the browser can load CesiumJS and Earth imagery. See [The dashboard](dashboard.md) for network requirements. |

See the [API reference](api.md) for function signatures and
[Data plotting](data_plotting.md) for other plot types and trace attributes.
