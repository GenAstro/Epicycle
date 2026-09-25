```@meta
CurrentModule = AstroModels
```

# AstroModels

The AstroModels package provides models of physical objects in space systems such as spacecraft and ground-stations for example. 

## Installation

Versions through 0.2.0 are in Julia's General registry. From the next version AstroModels is
released under the Gen Astro Source Available License, which General does not carry, so later
versions come from the Gen Astro registry. Add it once, then install as usual:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url = "https://github.com/GenAstro/GenAstro.git"))
Pkg.add("AstroModels")
```

General is still required, since these packages depend on packages registered there. Installing
without the Gen Astro registry resolves to 0.2.0, the last version General carries, and
reports nothing about the newer ones.

## Quick Start

```julia
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse

# A spacecraft carries its state, epoch, coordinate system, and physical data.
epoch = Time("2015-09-21T12:23:12", TAI(), ISOT())
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time = epoch,
    mass = 1000.0,
    coord_sys = CoordinateSystem(earth, ICRF()),
    name = "MySat"
)

# A ground station uses geodetic degrees, altitude in km, and an elevation mask.
station = GroundStation(
    name = "DSS-14",
    body = earth,
    latitude = 35.4267,
    longitude = -116.89,
    altitude = 1.0,
    min_elevation = 10.0
)

# The station state defaults to GCRF axes; visibility uses the same inertial axes.
r_station, v_station = get_state(station, epoch)
visible = is_visible(station, to_posvel(sc)[1:3], epoch)
```
`Spacecraft` stores its state as an `OrbitState` and accepts concrete state representations such as `CartesianState` and `KeplerianState` at construction. See [State](spacecraft_state.md) for the representation and conversion rules.

## Spacecraft

- [State](spacecraft_state.md) - Orbital state representations
- [Time](spacecraft_time.md) - Epoch and time scales
- [Mass](spacecraft_mass.md) - Spacecraft mass
- [Coordinate System](spacecraft_coord_sys.md) - Reference frames
- [CAD Model](spacecraft_cad_model.md) - 3D visualization
- [History](history.md) - Trajectory data
- [Reference](reference.md) - Complete API documentation

## Ground Stations

A `GroundStation` is a fixed site on Earth. Its location is geodetic latitude and longitude in degrees and altitude in km above the reference ellipsoid, which takes Earth's equatorial radius and flattening from AstroUniverse. `get_state` returns the station's position and velocity in km and km/s, in GCRF by default or in any axes named by the `axes` keyword. `is_visible` takes a spacecraft position in GCRF and compares its elevation, measured from the geodetic vertical at the station, with the station's `min_elevation` in degrees. The default cutoff of -90° makes every position visible. The constructor raises an `ArgumentError` for a latitude outside [-90°, 90°] or a cutoff outside [-90°, 90°). Stations on bodies other than Earth are not supported yet.

```julia
using AstroModels, AstroUniverse, AstroEpochs, AstroFrames

# Geodetic degrees, altitude in km, elevation cutoff in degrees
madrid = GroundStation(name = "DSS-63", body = earth,
                       latitude = 40.43, longitude = -4.25, altitude = 0.8,
                       min_elevation = 10.0)

epoch = Time("2020-06-01T00:00:00", UTC(), ISOT())

# Position and velocity in km and km/s, in GCRF unless other axes are named
r, v = get_state(madrid, epoch)
r_fixed, v_fixed = get_state(madrid, epoch; axes = ITRF())   # v_fixed is zero

# Visibility takes a GCRF position
is_visible(madrid, 1.05 .* r, epoch)    # overhead: true
is_visible(madrid, -r, epoch)           # the far side of Earth: false
```
