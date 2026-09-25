# Coordinate System

The `coord_sys` field specifies the coordinate system (origin body and axes) for the spacecraft's state.

## Basic Usage

```julia
using AstroModels, AstroFrames, AstroUniverse

# Earth-centered inertial (default)
sc = Spacecraft(
    coord_sys = CoordinateSystem(earth, ICRF())
)

# Mars-centered inertial
sc = Spacecraft(
    coord_sys = CoordinateSystem(mars, ICRF())
)

# Moon-centered inertial
sc = Spacecraft(
    coord_sys = CoordinateSystem(moon, ICRF())
)
```

## Accessing Coordinate System

<!-- doc-continue -->
```julia
# Get origin body
origin = sc.coord_sys.origin  # CelestialBody

# Get axes type
axes = sc.coord_sys.axes  # AbstractAxes
```
