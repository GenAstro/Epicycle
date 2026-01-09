# Coordinate System

The `coord_sys` field specifies the coordinate system (origin body and axes) for the spacecraft's state.

## Basic Usage

```julia
using AstroFrames, AstroUniverse

# Earth-centered inertial (default)
sc = Spacecraft(
    coord_sys = CoordinateSystem(earth, ICRFAxes())
)

# Mars-centered inertial
sc = Spacecraft(
    coord_sys = CoordinateSystem(mars, ICRFAxes())
)

# Moon-centered inertial
sc = Spacecraft(
    coord_sys = CoordinateSystem(moon, ICRFAxes())
)
```

## Accessing Coordinate System

```julia
# Get origin body
origin = sc.coord_sys.origin  # CelestialBody

# Get axes type
axes = sc.coord_sys.axes  # AbstractAxes
```
