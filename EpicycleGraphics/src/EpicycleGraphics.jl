# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
The EpicycleGraphics module renders spacecraft trajectories, celestial bodies and star fields
in an interactive native 3D window, built on GLMakie.

A view is a [`View3D`](@ref), which holds a coordinate system and the spacecraft drawn in it.
Spacecraft are added with [`add_spacecraft!`](@ref) and the window opens with
[`display_view`](@ref); the trajectory comes from each spacecraft's recorded history, so
propagate first and render afterwards.

```julia
using Epicycle, EpicycleGraphics

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]),
                 time  = Time("2015-09-21T00:00:00", TAI(), ISOT()))
prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                       IntegratorConfig(DP8(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))

view = View3D()                      # Earth-centred ICRF by default
add_spacecraft!(view, sat)
display_view(view)
```

Epicycle draws in two places and they are not interchangeable. This package opens a native
window through GLMakie, which is the heavier dependency and is why it is a separate package
rather than part of the `Epicycle` umbrella. `EpicycleIO` draws in a browser instead — Plotly
for data and Cesium for 3D through `orbitview` — and needs no graphics stack. Choose by where
the picture has to appear.
"""
module EpicycleGraphics

using GLMakie
using FileIO
using MeshIO
using GeometryBasics
using Random
using LinearAlgebra
using Colors: RGB, N0f8

using AstroModels: Spacecraft
using AstroUniverse: CelestialBody, earth, moon
using AstroFrames: CoordinateSystem, ICRF

export View3D, add_spacecraft!, display_view

include("view3d.jl")
include("trajectories.jl")
include("bodies.jl")
include("stars.jl")
include("planes.jl")
include("labels.jl")
include("models.jl")

end # module EpicycleGraphics
