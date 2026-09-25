```@meta
CurrentModule = EpicycleIO
```

# 3D Plotting

`orbitview` draws a propagated trajectory about the central body with a clock that can be played,
paused and run at a chosen rate. The view is a panel on the same dashboard as the plots, so a
trajectory and the quantities describing it are seen on one page.

This is the one part of EpicycleIO that takes a spacecraft rather than arrays. A trajectory drawn
about a body requires a reference frame and an epoch, and neither is styling. The example below
propagates a spacecraft and draws it, and it runs as written.

```julia
using Epicycle
using EpicycleIO

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 name  = "Explorer")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))

# Everything the view needs is read from the spacecraft
orbitview("Trajectory", sat)
```

## What Is Read From the Spacecraft

The trajectory comes from `sat.history`, the epoch from the earliest recorded sample, one
coloured arc from each propagation segment, and the label on the entity from `sat.name`. The
palette is the one the plots use, so a trajectory and a plot of the same run match.

Segments are what make a maneuver visible. A coast, a burn and a coast leave two propagation
segments in the history, and they draw in different colours, so the same maneuver appears as a
colour change in the view and as a step in the plot beside it.

## Maneuvers

An impulsive maneuver is recorded as a one-sample segment in the history and is marked with a
point labelled with its magnitude. Nothing is passed in for it. The magnitude is the velocity
change between the arcs either side of the marker, so the number on the label is read back out of
the trajectory rather than taken from whatever commanded it.

```julia
using Epicycle
using EpicycleIO

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 1300.0, 0.0, 7.35, 1.0]),
                 time  = Time("2020-01-01T00:00:00.000", UTC(), ISOT()),
                 name  = "Explorer")

prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))

# Coast, burn, coast: two arcs in two colours with a labelled point between them
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.08))
maneuver!(sat, ImpulsiveManeuver(axes = VNB(), element1 = 0.20,
                                 element2 = 0.0, element3 = 0.0))
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 0.25))

orbitview("Transfer", sat)
```

Because the magnitude is differenced from the recorded velocities, a label that disagrees with
the value commanded means the wrong pair of arcs was found rather than that the label is
cosmetic.

## Playback and Appearance

Playback defaults to fitting the whole trajectory into about two minutes of real time, clamped so
that a short arc is not unwatchably slow and a long one is not a blur.

```@raw html
<!-- doc-fragment -->
```
```julia
orbitview("Trajectory", sat; speed = 3600)              # one hour of flight per second
orbitview("Trajectory", sat; linewidth = 3)             # arc width in pixels
orbitview("Trajectory", sat; colors = [(255, 80, 80, 255), (55, 200, 255, 255)])
```

Play, pause, reset, slower, faster and stars sit in a bar under the view, with the clock time on
the right. Holding shift while clicking slower or faster changes the rate in fives rather than
doubles.

## The Camera

The camera is fixed in the inertial frame. The orbit therefore holds still against the stars and
the central body rotates beneath it, which is almost always what is wanted when looking at a
trajectory, and it is why the view does not follow the spacecraft by default.

Positions are sent to the browser in the inertial frame and Cesium applies the rotation to the
fixed frame itself, once per rendered frame. The epoch is converted to UTC on the way out
whatever scale the history was recorded in, because Cesium reads CZML times as UTC.

## One Spacecraft per View

A view holds one spacecraft. A second spacecraft is given its own panel.

```@raw html
<!-- doc-fragment -->
```
```julia
orbitview("Chief", chief)
orbitview("Deputy", deputy)
```

`orbitview!` exists and raises an error. Drawing several spacecraft in one view requires the CZML
writer to merge entity sets, which it does not yet do, so it says so rather than silently drawing
one and dropping the other.

## Central Bodies

Earth is currently the only supported central body. The body is taken from the frame the history
was recorded in rather than assumed, and nothing in the naming or the design refers to Earth, so
adding the Moon or a planet is a configuration change rather than a rewrite. The name
`orbitview` was chosen over `globe` for the same reason: the first interplanetary cruise makes
"globe" the wrong word.

## Implementation

The trajectory is written as [CZML](https://github.com/AnalyticalGraphicsInc/czml-writer/wiki), a
time-dynamic scene format, and the page hands it to
[Cesium](https://cesium.com/platform/cesiumjs/) to render.

Cesium is downloaded the first time a 3D panel appears and not before, so a dashboard of plots
does not pull it. It is about 5 MB and your browser keeps it, so you pay that once rather than per
session. The first 3D view on a machine with no network connection will not draw; see
[the dashboard page](dashboard.md) for what is cached and for how long.
