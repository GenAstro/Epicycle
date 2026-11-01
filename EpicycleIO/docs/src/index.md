```@meta
CurrentModule = EpicycleIO
```

```@setup io
using EpicycleIO
EpicycleIO.auto_open!(false)
```

# EpicycleIO

EpicycleIO provides tools for visualizing and exporting Epicycle results. Use it
to plot data, animate a spacecraft trajectory, write tabular reports, and exchange
ephemeris files with other software.

Plots and animated trajectories appear together in a browser dashboard served from
your machine. The dashboard updates as your Julia session produces results, so you
can inspect a propagation or follow a solver while it runs.

| Task | Function | Input |
|:--|:--|:--|
| Plot data | `xyplot` and other plotting functions | Arrays of values |
| Animate a trajectory | `orbitview` | A spacecraft with a recorded propagation history |
| Write a tabular report | `report` | Named columns of values |
| Read an ephemeris | `read_ephemeris` | A CCSDS OEM, STK `.e`, or JPL SPK file |
| Export a trajectory | `write_ephemeris` | A spacecraft or spacecraft history |

Plotting and reporting work with arrays, including data from sources other than
Epicycle. For spacecraft results, use `history` and `Calc` to extract the quantities
you need. Trajectory views and ephemeris functions work directly with spacecraft
or history objects so they can use the associated epochs and coordinate systems.

## Installation

EpicycleIO is distributed through the Gen Astro registry under the Gen Astro
Source Available License. Add the registry, then install the package:

```julia
using Pkg
pkg"registry add General https://github.com/GenAstro/GenAstroRegistry.git"
pkg"add EpicycleIO"
```

The registry command includes General because EpicycleIO's dependencies are
registered there. You only need to add the registries once for a Julia installation.

Load `EpicycleIO` to plot arrays or write reports. The spacecraft examples on this
page also use the `Epicycle` umbrella package, which must be installed in your
active Julia environment. Epicycle re-exports EpicycleIO's public functions.

## Quick start: plot data

Create a named panel with `xyplot`, then set its axis labels with `panel!`:

```@example io
using EpicycleIO

hours = collect(0.0:0.25:12.0)
altitude_km = 400.0 .+ 25.0 .* sin.(hours)

xyplot("Altitude", hours, altitude_km; name = "Spacecraft A")
panel!("Altitude"; xaxis_title = "Hours from epoch", yaxis_title = "Altitude (km)")
```

The first plot opens the dashboard in your browser.
Subsequent plots update the same dashboard. A leading string selects the panel;
plotting again with the same name replaces that panel's plotted data. Use the `!`
form to add a series:

```@example io
xyplot!("Altitude", hours, altitude_km .+ 30; name = "Spacecraft B")
```

Trace options such as `name`, `mode`, and `line_color` follow Plotly conventions.
Use `panel!` for layout options such as titles, axis labels, and ranges. Nested
options use underscores: `yaxis_title` sets Plotly's `yaxis.title` property.

A matrix draws one series per column, and a vector of vectors draws one series per
component. See [Data plotting](data_plotting.md) for array shapes, trace types,
styling, and subplot layouts.

## Work with a propagated trajectory

The following example propagates an Earth-orbiting spacecraft for one hour. It
then displays the trajectory and plots the recorded distance from Earth's center.
Initial position and velocity components are in km and km/s.

```@example io
using Epicycle

sat = Spacecraft(
    state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
    time = Time("2026-10-21T00:00:00", UTC(), ISOT()),
    name = "Explorer",
)

prop = OrbitPropagator(
    ForceModel(PointMassGravity(earth, ())),
    IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12),
)
propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 3600.0))

orbitview("Trajectory", sat)

times, radius_km = history(Calc(epoch, sat), Calc(position_magnitude, sat))
hours_from_epoch = [(t - first(times)) * 24 for t in times]

xyplot("Radius", hours_from_epoch, radius_km; name = "Explorer")
panel!("Radius"; xaxis_title = "Hours from epoch", yaxis_title = "Radius (km)")
```

Subtracting two `Time` values gives an interval in days; multiplying by 24 produces
the elapsed hours used on the horizontal axis.

`orbitview` reads the trajectory, epochs, and reference frame from the spacecraft's
history. Its playback controls let you animate the orbit, and impulsive maneuvers
are marked with their velocity changes. Earth is currently the only supported
central body, and each trajectory panel holds one spacecraft. See
[3D plotting](plotting_3d.md) for playback and display options.

## Export results

Choose a tabular report when you want columns of calculated quantities. Choose an
ephemeris when you want a trajectory with its epochs, coordinate system, and
interpolation information.

### Tabular reports

`report` writes one column for each keyword. All columns must have the same number
of rows. Include units in the column names, and convert epochs to the representation
you want readers to use.

```@example io
output_dir = mktempdir()
report_path = joinpath(output_dir, "radius.txt")
report(report_path; hours_from_epoch = hours_from_epoch, radius_km = radius_km)
```

An existing file is overwritten. Vector-valued columns expand into component
columns, such as `position_1`, `position_2`, and `position_3`. For absolute epochs,
convert `Time` values explicitly, for example with `t.isot` or `t.jd`, and identify
the time system in the column heading. See [Reporting](reporting.md) for formatting
and column requirements.

### Ephemeris files

Write the spacecraft's recorded history, then read it back as a `SpacecraftHistory`:

```@example io
ephemeris_path = joinpath(output_dir, "explorer.oem")
write_ephemeris(ephemeris_path, sat; step = 60.0)

eph = read_ephemeris(ephemeris_path)
coverage(eph)
```

The extension selects CCSDS OEM (`.oem`), STK ephemeris (`.e`), or JPL SPK (`.bsp`);
SPK also needs the object's `naif_id`.
The writer's `step` specifies the output sampling interval in seconds. Without it,
the writer keeps the source epochs. An existing output file is replaced.

Use `state_at` to evaluate a loaded history within its coverage. You can also pass
the history to `write_ephemeris` to convert it to another supported format. See
[Ephemeris files](ephemeris.md) for examples, format requirements, metadata,
interpolation, and diagnostics.

## Manage the dashboard

Named panels let you keep related plots and trajectory views together. Use
`open_dashboard` to show selected panels, and clear panels when you no longer
need their contents:

```julia
open_dashboard("Altitude", "Radius")
clear!("Altitude")  # Remove the contents of one panel.
clear_all!()        # Remove all panels.
```

For scripts or automated runs, disable automatic browser opening before plotting:

```julia
EpicycleIO.auto_open!(false)
```

Plots use Plotly.js, and trajectory views use CesiumJS. The browser loads these
libraries from external servers, so an internet connection is needed if they are
not already cached. The dashboard itself runs locally and does not require an
account or API key. See [The dashboard](dashboard.md) for panel layout, live
updates, and session controls.

## Guides and reference

- [Data plotting](data_plotting.md): array shapes, plot types, styling, and subplots.
- [3D plotting](plotting_3d.md): trajectory animation and display options.
- [The dashboard](dashboard.md): panels, browser views, and live updates.
- [Reporting](reporting.md): column names, units, epoch formatting, and text output.
- [Ephemeris files](ephemeris.md): reading, evaluating, comparing, and exporting trajectories.
- [API reference](api.md): function signatures, keywords, and detailed behavior.

In the Julia REPL, enter help mode with `?`, then type a function name such as
`xyplot`, `orbitview`, or `write_ephemeris` to read its documentation.
