```@meta
CurrentModule = AstroUniverse
```

# AstroUniverse

The AstroUniverse package provides celestial-body definitions, physical constants, SPICE
ephemerides, and orientation models. Other packages such as AstroFrames, AstroProp, and
AstroSolve use these to build coordinate systems, evaluate forces, and estimate trajectories.

The `CelestialBody` type holds a body's gravitational parameter, radius, flattening, and NAIF ID.
The package includes the Sun, the Moon, all eight planets, and Pluto, and you can define your own
bodies. Ephemeris functions such as `translate` and the lunar orientation models read SPICE kernels,
which are downloaded on first use and cached for later sessions. Body orientation is a model set on
the body. Earth orientation uses data published by the International Earth Rotation and Reference
Systems Service (IERS).

AstroUniverse defines the interface for corrections to Earth orientation. Precision corrections
are available in the `EpicycleEnterprise` package.

## Quick Start

Start by reading a built-in body's properties or defining a body of your own. Distances are in
kilometers, and gravitational parameters are in km³/s².

```julia
using AstroUniverse

# Read the constants of a built-in body
earth.mu                    # km³/s²
mars.equatorial_radius      # km
venus.naifid

# Define a body
phobos = CelestialBody(name = "Phobos", mu = 7.0875e-4, equatorial_radius = 13.0, naifid = 401)
```

To query an ephemeris or evaluate an orientation model, supply an AstroEpochs `Time`, which states
its time scale. The first ephemeris call may download about 111 MB of kernel data.

```julia
using AstroEpochs, AstroUniverse

t = Time("2020-01-01T00:00:00", UTC(), ISOT())

# The Moon's position relative to Earth, km
translate(earth, moon, t)

# The Moon's position and velocity relative to Earth, km and km/s
translate_state(earth, moon, t)

# Transform a state from ICRF to Mars-fixed axes
body_axes_rotation(orientation_model(mars), mars.naifid, t)
```

## Celestial Bodies

The built-in bodies are `sun`, `mercury`, `venus`, `earth`, `moon`, `mars`, `jupiter`, `saturn`,
`uranus`, `neptune`, and `pluto`.

Each `CelestialBody` stores a name, gravitational parameter, equatorial radius, flattening, and NAIF
ID. To create a body, provide `name`, `mu`, `equatorial_radius`, and `naifid`. Flattening defaults to
zero, a sphere, so nothing of Earth's is assumed.

You can change a body's gravitational parameter with `set_gravparam!` or by assigning its `mu` field.
The gravitational parameter and radius must be finite and positive, and flattening must be finite and
in [0, 1). These checks apply when you create a body and whenever you change one of these fields.

```julia
using AstroUniverse

# Read a body's value
get_gravparam(moon)

# Use a different value for the rest of the session
set_gravparam!(moon, 4902.8)
```

Changing a built-in body changes the shared object for the current Julia session, so other code using
that body sees the new value.

The table lists the built-in constants and their sources.

| Body | μ, km³/s² | Equatorial radius, km | μ source | Shape source |
|:-----|--:|--:|:--|:--|
| Sun | 1.3271244004127942e11 | 695700 | `gm_de440.tpc` | `pck00011.tpc` |
| Mercury | 22031.868551400003 | 2440.53 | `gm_de440.tpc` | `pck00011.tpc` |
| Venus | 324858.592 | 6051.8 | `gm_de440.tpc` | `pck00011.tpc` |
| Earth | 398600.4418 | 6378.137 | WGS84 | WGS84 |
| Moon | 4902.8001184575496 | 1737.4 | `gm_de440.tpc` | `pck00011.tpc` |
| Mars | 42828.375815756102 | 3396.19 | `gm_de440.tpc`, system | `pck00011.tpc` |
| Jupiter | 1.2671276409999998e8 | 71492 | `gm_de440.tpc`, system | `pck00011.tpc` |
| Saturn | 3.7940584841799997e7 | 60268 | `gm_de440.tpc`, system | `pck00011.tpc` |
| Uranus | 5.7945563999999985e6 | 25559 | `gm_de440.tpc`, system | `pck00011.tpc` |
| Neptune | 6.8365271005803989e6 | 24764 | `gm_de440.tpc`, system | `pck00011.tpc` |
| Pluto | 975.5 | 1188.3 | `gm_de440.tpc`, system | `pck00011.tpc` |

A body's flattening comes from the polar and equatorial radii of the same source. Earth's is
1/298.257223563.

## Choosing Gravitational Parameters

The built-in gravitational parameters are set up for Earth orbits. Earth uses the WGS84 value. Mars,
Jupiter, Saturn, Uranus, Neptune, and Pluto use system values that include their moons, since from
Earth a planet and its moons act as one mass. Two other kinds of analysis need different values.

If your analysis needs constants consistent with DE440, such as work centered on the solar system
barycenter, use Earth's DE440 value.

```julia
using AstroUniverse

# Earth's μ from DE440
set_gravparam!(earth, 398600.43550702266)
```

This choice matters. Starting from the same state with point-mass forces from Earth and ten third
bodies, moving Earth's μ from 398600.4418 to 398600.43551 changes the position after 7 days by 145 m
for an orbit with a semimajor axis of 7000 km, and by 59 m at 42164 km. Choose one value and use it
throughout an analysis.

If you model a planet's moons as separate bodies, as you usually would in orbit about that planet, use
the planet's own value.

```julia
using AstroUniverse

# Mars's own μ, without Phobos and Deimos
set_gravparam!(mars, 42828.37362069909)
```

| Body | Built-in μ, km³/s² | Own μ, `gm_de440.tpc`, km³/s² |
|:-----|--:|--:|
| Earth | 398600.4418 | 398600.43550702266 |
| Mars | 42828.375815756102 | 42828.37362069909 |
| Jupiter | 1.2671276409999998e8 | 1.266865319003704e8 |
| Saturn | 3.7940584841799997e7 | 3.793120623436167e7 |
| Uranus | 5.7945563999999985e6 | 5.793951256527211e6 |
| Neptune | 6.8365271005803989e6 | 6.835103145462294e6 |
| Pluto | 975.5 | 869.6138177608748 |

The Sun, Mercury, Venus, and the Moon keep their built-in values in both cases.

## Ephemeris Translation

`translate(from, to, t)` returns the position of `to` relative to `from`, and
`translate_state(from, to, t)` returns its position and velocity. Results use ICRF axes, with position
in kilometers and velocity in kilometers per second. They are geometric: no light-time or aberration
correction is applied.

Supply the epoch as an AstroEpochs `Time`, in any time scale. A bare Julian date is an error, because
a number does not say which scale it is on.

The default ephemeris covers 1950 to 2100, and a date outside the loaded kernels' coverage raises an
error from SPICE. A custom body needs a loaded kernel that carries its NAIF ID. The Uranus and Neptune
centers differ from JPL Horizons by 60 km to 4190 km at the dates compared, as described under
[Validation](#Validation).

```julia
using AstroEpochs, AstroUniverse

t = Time("2020-01-01T00:00:00", UTC(), ISOT())

# Mars's position relative to Earth, km
translate(earth, mars, t)

# Mars's position and velocity relative to Earth, km and km/s
translate_state(earth, mars, t)
```

When querying several bodies at one epoch, convert the epoch once and reuse it, as in `t_tdb = t.tdb`.
Both functions convert to TDB on each call otherwise. They are safe to call from several threads, and
the SPICE queries run one at a time.

## Body Orientation

An orientation model defines a body's rotating axes as a function of time, and each model is named for
the publication that defines it. Except for Earth's, a model takes the epoch as an AstroEpochs `Time`
and returns the 6×6 matrix that transforms position and velocity from ICRF to the body's axes,
including the contribution from the axes' rotation.

Use `orientation_model(body)` to read a body's model and `set_orientation!(body, model)` to change it.
The selection is stored by NAIF ID and shared throughout the Julia session. For built-in models,
`set_orientation!` checks that the model supports the body. Custom models are accepted without being
evaluated, and a SPICE frame is checked when the rotation is evaluated.

| Body | Default model | Other included models |
|:-----|:--------|:-------------|
| Sun, planets other than Earth, Pluto | `IAU2015()` | `IAU1991()` for Mars |
| Moon | `LunarPA()`, principal axes (DE440) | `LunarME()`, mean-Earth axes |
| Earth | the frame theory, `IAU2006()` | `FK5()` |

```julia
using AstroEpochs, AstroUniverse

t = Time(2458849.5, TDB(), JD())

# Read Mars's model and its rotation
model = orientation_model(mars)
body_axes_rotation(model, mars.naifid, t)

# Use the mean-Earth axes for the Moon
set_orientation!(moon, LunarME())
```

Earth's orientation is the frame theory, so `set_orientation!(earth, FK5())` and
`set_frame_theory!(FK5())` are the same call. Its rotation is the ITRF chain, which needs UT1 and TT
as well as TDB. AstroFrames evaluates it, and its `body_fixed_rotation(model, naifid, epoch)`
evaluates any model, Earth's included.

Some data products, such as gravity fields, are defined in axes that differ from the body's default.
Specify the orientation the data product requires rather than changing the body's default for
unrelated calculations.

A custom body can use a frame from a loaded SPICE kernel, `SpiceOrientation("…")`, or a model written
in Julia. A subtype of `AbstractOrientationModel` represents a way of orienting a body, and it
implements one method. Writing the fields with a type parameter lets an estimator store dual numbers
in them.

```julia
using AstroEpochs, AstroUniverse

struct SimpleSpin{T<:Real} <: AbstractOrientationModel
    pole_ra::T      # deg
    pole_dec::T     # deg
    pm0::T          # deg at J2000
    spin_rate::T    # deg/day
end

function AstroUniverse.body_axes_rotation(m::SimpleSpin, naifid, t::Time)
    d = t.tdb.jd - 2451545.0
    Ẇ = deg2rad(m.spin_rate)
    return pole_axes_rotation(deg2rad(m.pole_ra), deg2rad(m.pole_dec),
                              deg2rad(m.pm0) + Ẇ * d, zero(Ẇ), zero(Ẇ), Ẇ / 86_400)
end

# Set the model on a body
bennu = CelestialBody("Bennu", 4.892e-9, 0.2825, 0.0, 2101955)
set_orientation!(bennu, SimpleSpin(85.46, -60.36, 89.6, 2011.145))
```

Once the model is set, coordinate systems and force models that read the body's orientation use it.

## Earth Orientation

The active frame theory selects the default Earth precession-nutation chain. `IAU2006()` is the
default, and `FK5()` selects the classical IAU-76/80 chain.

```julia
using AstroUniverse

# Read the active theory
frame_theory()

# Select FK5
set_frame_theory!(FK5())
```

Earth orientation parameters (EOP) describe measured variations in Earth's rotation: UT1−UTC, polar
motion, length of day, and celestial-pole offsets. They enter every transformation to or from an
Earth-fixed frame. `IAU2006()` uses the IERS IAU 2000A series and `FK5()` the IAU 1980 series. Each
series has its own table, so both can be loaded in one session.

A table loads on first use and stays in memory for the session. AstroUniverse reads it from a local
cache, and downloads a new copy if the cached file is missing, unreadable, or a week old. It tries the
IERS data center first, then the US Naval Observatory, which publishes the same series in a fixed-width
format.

A download replaces the cached file only after it passes validation, which includes a check that its
UT1−UTC values reach today. If no source supplies usable data, AstroUniverse uses a readable cached
copy and warns about its age. With no usable copy, the error names the sources and how to load a file
by hand.

`eop_refresh!` downloads and installs the latest table immediately. It throws if no source supplies
usable data, and leaves the cached file unchanged.

`eop_load` reads either format from a local file. It applies the format and finite-value checks a
download gets and refuses a file of the other theory's series, but it does not require the file to
reach today, so a historical file loads.

For repeatable results, load a saved file with `eop_load` or install a table with `set_eop!` before
the first transformation, for every series your calculation uses. This prevents EOP downloads. For a
run with no network at all, the SPICE kernels described below must also be in the cache.

A finals file carries the rapid values of IERS Bulletin A and, for past dates, the final values of
Bulletin B. `eop_load` reads Bulletin A by default. With `bulletin = :B`, it uses Bulletin B where
available and Bulletin A elsewhere. Length of day and the uncertainties always come from Bulletin A.

```@raw html
<!-- doc-fragment -->
```
```julia
# Read the active theory's table, loaded on first use
eop()

# Download the latest series now
eop_refresh!()
eop_refresh!(; theory = FK5())

# Read a local file, CSV or fixed width, with no network
eop_load("finals2000A.all"; theory = IAU2006())
eop_load("finals2000A.all"; bulletin = :B)

# Install a table already in hand
set_eop!(table)
```

Outside a table's span, `eop` holds the value at the boundary. A query outside the polar-motion span
warns once per session. Length of day and the celestial-pole offsets can end earlier, and are held
without a warning, so check their coverage if your analysis depends on them.

An `AbstractEopCorrection` adds terms the daily tables do not carry. A subtype implements
`eop_correction(model, jd_utc, jd_tt)`, and `set_eop_corrections!` installs a list of them. A custom
correction is responsible for returning valid values in the documented units, because a non-finite
value propagates into the computed Earth orientation.

!!! note "Enterprise Feature"
    `IersSubdaily`, the ocean-tide, libration, and zonal corrections, requires the `EpicycleEnterprise` package. The open package defines `AbstractEopCorrection`, the interface it plugs into.

## SPICE Kernels

AstroUniverse reads ephemerides and lunar orientation from SPICE kernels. The default set is a
leap-second kernel from NAIF, the 1950–2100 DE440 ephemeris published with Epicycle, and NAIF's lunar
orientation kernels. They are stored between Julia sessions with Scratch.jl, in the folder
`get_spice_directory()` returns.

The download happens at the first call that needs a kernel, directly or through a package that uses
AstroUniverse, and is about 111 MB, most of it the ephemeris.

Downloads from registered kernel sources are verified against their SHA-256 checksums, and a cached
registered file is checked by size and fetched again if the size is wrong. A download from any other
URL is not checked.

`using AstroUniverse` itself never reaches the network. Later sessions can use the cached kernels
without downloading them again, as long as all the files are present at the right size. If a download
fails, the error names each file, where it comes from, and the folder to place it in by hand.

```julia
using AstroUniverse

# Find the stored kernels
get_spice_directory()
cached_spice_kernels()

# List the kernels loaded in this session
loaded_spice_kernels()

# Download another kernel and load it
download_spice_kernel(
    "de440s.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
)
load_spice_kernel("de440s.bsp")
```

`unload_spice_kernel` and `unload_all_spice_kernels` remove kernels from the current session without
deleting downloaded files. `unload_all_spice_kernels` also turns off automatic loading of the default
set for the rest of the session, so afterward you download and load every kernel your calculation needs.

`KERNEL_SOURCES` lists the kernels the package knows by name, for
`download_spice_kernel("pck00011.tpc")` and the like. The kernel pool belongs to the process, and a
kernel loaded later takes precedence over earlier ones for the bodies and times it covers.

## Texture Maps

Texture maps for the Sun, Moon, and planets are distributed in the package's `data` directory, and
Pluto has none. Earth's image is from the NASA Earth Observatory Blue Marble collection. The others
are provided by Solar System Scope under the Creative Commons Attribution 4.0 license.

## Conventions

- Distances are in km, velocities in km/s, and gravitational parameters in km³/s².
- Ephemeris queries and orientation models take AstroEpochs `Time` values and evaluate them in TDB. Earth's transformation chain also uses UT1 and TT.
- Inertial axes are ICRF. SPICE labels them `J2000`, and they are not FK5 mean equator and equinox of J2000.
- `translate(from, to, t)` is the position of `to` relative to `from`.
- Orientation model angles are in radians and rates in rad/s when passed to `pole_axes_rotation`.
- Polar motion is in arcseconds, UT1−UTC in seconds, length of day in milliseconds, and celestial-pole offsets in milliarcseconds.
- Bodies, orientation selection, frame theory, EOP corrections, and the SPICE kernel pool are shared across a Julia session.

## Differentiability

`IAU2015`, `IAU1991`, and EOP values are differentiable with ForwardDiff, and tests compare the derivatives with central differences. A custom model written in Julia is differentiable when its code and field types carry derivative information, and installed EOP corrections need to do the same.

Calls that read SPICE, which are `translate`, `translate_state`, `LunarPA`, `LunarME`, and `SpiceOrientation`, have no derivative through the date and raise an error when given one.

```julia
using AstroEpochs, AstroUniverse, ForwardDiff

# How one element of Mars's rotation changes per day
ForwardDiff.derivative(2458849.5) do jd
    body_axes_rotation(IAU2015(), mars.naifid, Time(jd, TDB(), JD()))[1, 1]
end
```

## Validation

The test suite compares positions with GMAT R2022a at two stored cases, Sun to Earth in 2015 and Moon to Earth in 1995, within 1 part in 10¹⁰.

Barycentric positions were compared once, outside the test suite, with JPL Horizons at 2020-01-01 TDB. Mercury, Venus, Mars, Jupiter, Saturn, Pluto, Earth, and the Sun agree to the 0.0001 km the comparison reported, and the Moon to 0.0017 km.

The Uranus and Neptune centers differ from JPL Horizons by 700 km and 240 km at 2020-01-01, by 360 km and 60 km at J2000, and for Uranus by 4190 km at 2050. The difference comes from how each source relates a planet's center to the solar system barycenter. Horizons follows the barycenter chain of its satellite solutions. The default kernel keeps DE440's barycenters, which agree with Horizons to the 0.001 km the comparison reported, and adds the offset to the center. The offsets themselves agree to 1.0 km and 0.12 km.

The default kernel was compared once with JPL's `de440.bsp` at 5000 epochs for each barycenter, the Sun, Mercury, Venus, the Moon, and Earth, and every difference is exactly zero. The tests check the Jupiter, Saturn, and Pluto center offsets against JPL Horizons within 0.01 km at 2020-01-01, and that every body resolves at both ends of the span.

`translate_state`'s velocity is compared with a central difference of position for the Moon, Mars, Jupiter, and Neptune, within 1e-7 km/s.

`IAU2015` is compared with NAIF's `pck00011.tpc` for eight bodies, the planets other than Earth and Pluto, at one epoch, within 1e-5 arcsec in rotation and 1 part in 10⁹ in rate. This checks the implementation of the same publication and does not establish how well the publication describes the body. `IAU1991` is compared with NAIF's `preliminary_pck00004.tpc` for Mars at three dates, within 1e-3 arcsec. `LunarPA` and `LunarME` are equal to the named SPICE frames.

Every built-in constant is compared with the source in the table under [Celestial Bodies](#Celestial-Bodies). EOP values are compared with published IERS rows within 1e-12 of each field's units, and the fixed-width and CSV files give identical values at the dates tested.

## References

- Acton, C. H. (1996), Ancillary data services of NASA's Navigation and Ancillary Information Facility, Planetary and Space Science 44, 65.
- Archinal, B. A., et al. (2018), Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2015, Celestial Mechanics and Dynamical Astronomy 130, 22.
- Davies, M. E., et al. (1992), Report of the IAU/IAG Working Group on Cartographic Coordinates and Rotational Elements of the Planets and Satellites: 1991, Celestial Mechanics and Dynamical Astronomy 53, 377.
- International Earth Rotation and Reference Systems Service, Earth orientation data, <https://www.iers.org>.
- National Geospatial-Intelligence Agency (2014), Department of Defense World Geodetic System 1984, NGA.STND.0036_1.0.0_WGS84.
- NAIF, generic kernels `gm_de440.tpc` and `pck00011.tpc`, <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/>.
- Park, R. S., et al. (2021), The JPL Planetary and Lunar Ephemerides DE440 and DE441, The Astronomical Journal 161, 105.

## API Reference

```@index
```

```@autodocs
Modules = [AstroUniverse]
Public  = true
Private = false
Order = [:type, :function, :macro, :constant]
```
