```@meta
CurrentModule = AstroUniverse
```

# AstroUniverse

The AstroUniverse package provides celestial-body definitions, physical constants,
SPICE ephemerides, and orientation models for astrodynamics applications. It ships
with definitions for the Sun, Moon, planets, and Pluto, and supports additional
bodies through `CelestialBody`.

Planet and Moon positions come from the DE440 ephemeris, downloaded the first time AstroUniverse loads. Planetary
orientations use the IAU 2015 model by default. Earth orientation is handled
separately through IERS Earth orientation parameters and a selectable frame theory.

## Quick Start

Built-in bodies expose their physical properties directly. Distances are in
kilometers, time intervals are in seconds, and gravitational parameters are in
km^3/s^2.

```julia
using AstroUniverse

earth.mu
mars.equatorial_radius
venus.naifid

phobos = CelestialBody(
    "Phobos",
    7.0875e-4,  # gravitational parameter [km^3/s^2]
    11.1,       # equatorial radius [km]
    0.0,        # flattening
    401,        # NAIF ID
)
```

The built-in bodies are `sun`, `mercury`, `venus`, `earth`, `moon`, `mars`,
`jupiter`, `saturn`, `uranus`, `neptune`, and `pluto`.

## Ephemeris Translation

`translate` and `translate_state` evaluate relative positions and states from
loaded SPICE kernels. Epochs are TDB Julian dates.

```julia
jd_tdb = 2458849.5

r_mars_from_earth = translate(earth, mars, jd_tdb)
x_moon_from_earth = translate_state(earth, moon, jd_tdb)
```

## Body Orientation

The Sun, planets other than Earth, and Pluto use published IAU orientation
polynomials. A custom body can use a Julia orientation model or a frame supplied
by a loaded SPICE kernel.

```julia
model = orientation_model(mars)
rotation = body_axes_rotation(model, mars.naifid, 2458849.5)
```

Earth and the Moon use dedicated frame models instead of the planetary
pole-and-prime-meridian model. AstroFrames provides those frame transformations.

## Earth Orientation

The active frame theory selects the default Earth precession-nutation chain.
`IAU2006()` is the default; `FK5()` selects the classical IAU-76/80 chain.

```julia
frame_theory()
set_frame_theory!(FK5())
```

Earth orientation parameters (EOP) are the measured corrections to Earth's rotation that the
IERS publishes: UT1−UTC and polar motion. They enter every transformation to or from an
Earth-fixed frame. Each frame theory reads its own table, `IAU2006()` the IAU 2000A series and
`FK5()` the IAU 1980 series.

A table loads the first time a transformation needs it and stays loaded for the session. Loading
goes through SatelliteToolboxTransformations, which keeps the IERS files in its own on-disk cache
and downloads a new copy when that cache is missing or out of date. `eop_refresh!` downloads the
latest series immediately. A run that must not reach the network, or must give the same numbers
every time, loads a file with `eop_load` or installs a table with `set_eop!` before its first
transformation.

```@raw html
<!-- doc-fragment -->
```
```julia
eop()                                              # the active theory's table, loaded on first use
eop_refresh!()                                     # download the latest IERS series now
eop_refresh!(; theory = FK5())                     # the same for the FK5 table

eop_load("finals2000A.all"; theory = IAU2006())    # a local IERS file, no network
set_eop!(table)                                    # a table already in hand
```

## SPICE Kernels

AstroUniverse downloads its default kernels once, verifies their SHA-256
checksums, and stores them between Julia sessions with Scratch.jl. The default
set contains a leap-second kernel, the 1950-2100 DE440 ephemeris, and the lunar
orientation kernels used by Moon frames.

The download happens the first time AstroUniverse loads, directly or through a
package that uses it, and is about 110 MB, most of it the ephemeris. That first
load needs a network connection and waits for the download; later sessions read
the stored kernels and work offline.

Additional kernels can be downloaded and loaded independently:

```julia
download_spice_kernel(
    "de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
)
load_spice_kernel("de440s.bsp")
```

`list_downloaded_spice_kernels` reports files stored on disk, while
`list_cached_spice_kernels` reports kernels loaded in the current process.
`unload_spice_kernel` and `unload_all_spice_kernels` remove kernels from the
current SPICE session without deleting downloaded files.

## Texture Maps

Small texture maps for the Sun and planets are distributed in the package's
`data` directory. The images are provided by Solar System Scope under the
Creative Commons Attribution 4.0 license.

## API Reference

```@index
```

```@autodocs
Modules = [AstroUniverse]
Public  = true
Private = false
Order = [:type, :function, :macro, :constant]
```
