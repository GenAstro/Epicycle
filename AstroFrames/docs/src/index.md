```@meta
CurrentModule = AstroFrames
```

# AstroFrames

AstroFrames transforms position and velocity between coordinate systems. Each coordinate system pairs an origin with an axes orientation. The package provides IAU 2006 and FK5 Earth frames, planetary body-fixed frames, lunar frames, and orbit-relative frames.

For a state conversion, use a `Coordinate` to carry the state, coordinate system, and epoch together. The same conversion interface accepts a `Spacecraft` or another object that supplies this information. For calculations on arrays, use the rotation and translation operators directly.

AstroFrames is tested against [pyerfa](https://github.com/liberfa/pyerfa) for Earth-frame transformations, and against SPICE for planetary body-fixed and lunar transformations. pyerfa is the Python interface to ERFA, the open-licence release of the IAU's SOFA library.

AstroFrames uses the Earth-frame models in [SatelliteToolboxTransformations.jl](https://juliaspace.github.io/SatelliteToolboxTransformations.jl/stable/). We thank the authors of that package for their high-quality work.

The implementation follows these references:

- Chagas, R. A. J., *SatelliteToolboxTransformations.jl*, v1.2.0, JuliaSpace. https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
- Petit, G. and Luzum, B., eds. (2010), *IERS Conventions (2010)*, IERS Technical Note No. 36, Verlag des Bundesamts für Kartographie und Geodäsie.
- Vallado, D. A. (2013), *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm Press / Springer.
- Acton, C. H. (1996), *Ancillary Data Services of NASA's Navigation and Ancillary Information Facility*, Planetary and Space Science, 44(1).

## Quick Start

Create a coordinate and convert it to Earth-fixed axes. Cartesian states contain position followed by velocity, `[x, y, z, vx, vy, vz]`, in km and km/s.

```@example frames1
using AstroFrames, AstroUniverse, AstroEpochs

# Define the source and target coordinate systems
cs_j2000 = CoordinateSystem(earth, MJ2000Eq())
cs_fixed = CoordinateSystem(earth, ITRF())

# Associate the state with its frame and epoch
epoch = Time(2458849.5, 0.0, TDB(), JD())
c     = Coordinate([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0], cs_j2000, epoch)

# Convert to Earth-fixed axes at the same epoch
c_fixed = Coordinate(c, cs_fixed)
```

In the REPL, `using InteractiveUtils; subtypes(AbstractAxes)` lists the axes types. Help mode entries such as `?MoonME` and `?axes_rotation` describe individual calls.

## Coordinate Systems

A `CoordinateSystem` pairs an origin with an axes orientation. Some axes require a specific origin: Earth axes require Earth, lunar axes require the Moon, and body-fixed axes require the body encoded by the axes. The constructor raises an `ArgumentError` for an invalid pairing.

```@example frames2
using AstroFrames, AstroUniverse, AstroModels

# Any body from AstroUniverse can be an origin.
cs_earth_fixed   = CoordinateSystem(earth, ITRF())
cs_moon_ecliptic = CoordinateSystem(moon, MJ2000Ec())

# Define Mars body-fixed coordinates
cs_mars_fixed = CoordinateSystem(mars, CelestialBodyFixed())

# Choose mean-Earth axes for lunar mapping or principal axes for gravity
cs_moon_me = CoordinateSystem(moon, MoonME())
cs_moon_pa = CoordinateSystem(moon, MoonPA())

# Use a spacecraft as both origin and reference orbit
chief   = Spacecraft(name = "chief")
cs_ric  = CoordinateSystem(chief, RIC())
cs_vnb  = CoordinateSystem(chief, VNB())
cs_lvlh = CoordinateSystem(chief, LVLH())

# Use named systems for common origin-and-axes pairings
EarthMJ2000Eq
MoonFixed

# Reject an origin incompatible with the axes
# CoordinateSystem(mars, GCRF())   # GCRF axes require an Earth origin
```

`CelestialBodyFixed` uses the body's selected orientation model. For Earth, it gives `ITRF` under the active frame theory. For the Moon, the default is `MoonPA`; calling `set_orientation!(moon, LunarME())` selects mean-Earth axes instead. Use `ITRF`, `MoonPA`, or `MoonME` explicitly when the frame must be independent of the body's model setting.

## Coordinates and Conversions

`Coordinate(c, target)` returns a new coordinate in the target frame at the original epoch. The conversion rotates the axes and translates the origin when needed. To obtain only the Cartesian state, use `CartesianState(c, target)`.

The same calls accept a `Spacecraft`. Other types can participate by implementing [`state_of`](@ref), [`frame_of`](@ref), and [`epoch_of`](@ref).

```@example frames3
using AstroFrames, AstroUniverse, AstroEpochs, AstroModels, AstroStates, AstroCallbacks

# Define the state in Earth-centered J2000 equatorial coordinates
epoch = Time(2458849.5, 0.0, TDB(), JD())
system = CoordinateSystem(earth, MJ2000Eq())
state = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
c = Coordinate(state, system, epoch)

# Convert to Earth-centered J2000 ecliptic coordinates
c_ecliptic = Coordinate(c, CoordinateSystem(earth, MJ2000Ec()))

# Convert the Coordinate to Moon-centered MoonME axes
c_moon = Coordinate(c, CoordinateSystem(moon, MoonME()))

# Convert a spacecraft state to Earth-fixed coordinates
sc = Spacecraft()
Coordinate(sc, CoordinateSystem(earth, ITRF()))

# Obtain a position vector or Cartesian state in the requested frame
position_vector(sc, CoordinateSystem(earth, ITRF()))
CartesianState(sc, CoordinateSystem(moon, MoonME()))

# Evaluate right ascension of the ascending node in J2000 ecliptic coordinates
raan(sc, CoordinateSystem(earth, MJ2000Ec()))

# Compute the rotation matrix from ICRF to ITRF
M = axes_rotation(ICRF(), ITRF(), epoch)
```

`axes_rotation` returns a 6×6 state transform with block form `[R 0; Ṙ R]`, also used by SPICE `sxform`. Multiplying it by a Cartesian state rotates position and velocity and includes the velocity contribution from moving axes. For Earth frames, `Ṙ` includes sidereal rotation but treats precession, nutation, polar motion, and obliquity of date as constant during the transformation.

## Earth Frame Theories

AstroFrames supports two Earth frame theories, selected with `set_frame_theory!`:

- IAU 2006 precession with IAU 2000A nutation uses `GCRF → CIRS → TIRS → ITRF`. It is the default and is recommended for most work.
- FK5 uses IAU 1976 precession and IAU 1980 nutation along `ICRF → MODEq → TODEq → PEF → ITRF`. Use it for consistency with legacy data and models.

Both apply IERS celestial pole offsets. Their ICRF-to-ITRF results agree within about 0.3 mas at the test epochs.

The IAU 2006 precession-nutation step evaluates the celestial intermediate pole (CIP) coordinates X, Y and the CIO locator s from a series of over a thousand terms. By default, AstroFrames interpolates these values from a cached 30-minute grid. The measured interpolation difference is under 6.3e-5 µas.

Interpolation is faster for closely spaced epochs, such as propagation steps. For scattered epochs, such as hourly output, direct series evaluation was about twice as fast in the performance checks. Select it with:

```@example frames4
using AstroFrames

set_cip_interpolation!(false)     # evaluate the series at every epoch
cip_interpolation()               # false
set_cip_interpolation!(true)      # the default
```

`MJ2000Eq` is the FK5 mean equator and equinox of J2000. It differs from the ICRF by the IERS frame bias, a constant rotation of about 23 mas (0.75 m at LEO, 4.7 m at GEO), and is used for data that arrives labelled as FK5 J2000. It is not on the path of the FK5 chain: the IERS offsets applied at the nutation step already absorb the bias, so the chain starts at the ICRF. `MJ2000Ec` is the mean ecliptic of J2000, the same frame rotated by the J2000 obliquity.

```text
            MJ2000Ec
               │ (obliquity)
            MJ2000Eq
               │ (frame bias, 23 mas)
               ICRF ─────────────── GCRF          IAU 2006/2000A
                │       (identity)    │
  FK5 / IAU 76-80                    CIRS
                │                     │
              MODEq ── MODEc         TIRS
                │                     │
              TODEq ── TODEc         ITRF
              │   │                   │
            PEF  TEME                 └──── shared endpoint
              │
            ITRF
```

## Orbit-Relative Axes

Orbit-relative axes include RIC, LVLH and VNB. When the origin carries a state, AstroFrames uses it as the reference orbit for conversions into and out of that frame. Otherwise the reference state is supplied explicitly, in ICRF, in km and km/s.

At the operator level, `reference_accel` is optional. Omitting it leaves the orbit plane fixed
for RIC and LVLH and gives VNB zero frame rate. Supplying acceleration in km/s² includes its
contribution to the axes rates. When a state-carrying origin supplies the reference orbit,
AstroFrames also supplies two-body acceleration when its central body's gravitational parameter
is available. Explicit reference data take precedence.

A zero position, zero velocity or radial reference orbit leaves these axes undefined. The call
warns and returns NaNs so an optimizer can continue exploring other states. Malformed reference
data throw an `ArgumentError` naming the field.

```@example frames5
using AstroFrames, AstroUniverse, AstroEpochs, AstroModels, AstroStates

epoch = Time(2458849.5, 0.0, TDB(), JD())

# Create two spacecraft and compute the deputy's state in the chief's RIC frame
chief  = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]),
                    time = epoch, name = "chief")
deputy = Spacecraft(state = CartesianState([7000.1, 0.0, 0.0, 0.0, 7.546, 0.0]),
                    time = epoch, name = "deputy")
Coordinate(deputy, CoordinateSystem(chief, RIC()))

# Rotate from ICRF to RIC with an explicit reference state
p = (; reference_state = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0])
M = axes_rotation(ICRF(), RIC(), epoch, p)
```

## Lunar Frames

AstroFrames provides two lunar body-fixed frames. `MoonME` is the mean Earth/mean rotation axis frame used for published lunar coordinates, including surface latitude and longitude. `MoonPA` is aligned with the Moon's principal axes of inertia and is the frame used for lunar gravity coefficients. Both include the full physical, forced, and free libration. They differ by a fixed rotation of about 104 arcseconds, or about 875 metres at the surface.

```@example frames6
using AstroFrames, AstroUniverse, AstroEpochs

epoch = Time(2458849.5, 0.0, TDB(), JD())

# Express a surface site in lunar mapping axes
site = Coordinate([1737.4, 0.0, 0.0, 0.0, 0.0, 0.0],
                  CoordinateSystem(moon, MoonME()), epoch)

# Express an orbital state in lunar principal axes for dynamics
orb = Coordinate([1837.4, 0.0, 0.0, 0.0, 1.63, 0.0],
                 CoordinateSystem(moon, MoonPA()), epoch)
```

AstroUniverse loads the default lunar orientation kernels and solar system ephemeris at startup. Using either lunar frame does not require a kernel path.

## Axes Reference

| Family | Members | Origin |
|---|---|---|
| Origin-agnostic inertial | [`ICRF`](@ref), [`MJ2000Eq`](@ref), [`MJ2000Ec`](@ref) | any |
| Earth inertial (IAU 2006) | [`GCRF`](@ref), [`CIRS`](@ref) | Earth |
| Earth inertial (FK5) | [`MODEq`](@ref), [`TODEq`](@ref), [`MODEc`](@ref), [`TODEc`](@ref) | Earth |
| Earth rotating | [`TIRS`](@ref), [`ITRF`](@ref), [`PEF`](@ref), [`TEME`](@ref) | Earth |
| Moon rotating | [`MoonPA`](@ref), [`MoonME`](@ref) | Moon |
| Body-fixed | [`CelestialBodyFixed{NAIFID}`](@ref) | matching body |
| Orbit-relative | [`RIC`](@ref), [`LVLH`](@ref), [`VNB`](@ref) | any |
| Ambient inertial | [`Inertial`](@ref) | any |

## EOP Data

AstroUniverse manages Earth orientation parameters (EOP), including data loading and corrections. Earth-frame transformations use those data at the requested epoch. Outside the loaded range, AstroUniverse warns and holds the endpoint values; those values are not predictions.

## Origins Reference

An origin is any `AbstractPoint`. Celestial bodies from AstroUniverse serve as origins for body-fixed and body-centered inertial frames. A spacecraft, or another point that carries a state, can serve as the origin and reference orbit for orbit-relative axes.

AstroFrames provides named coordinate systems for common pairings: `EarthMJ2000Eq`, `EarthMJ2000Ec`, `EarthFixed`, `EarthTODEq`, `EarthICRF`, `MoonFixed`, `MoonPrincipalAxes`, and `SunMJ2000Ec`.

## Validation

AstroFrames is tested against external references where they are available.

| What | Checked against | Agreement |
|---|---|---|
| Earth chains, both theories | pyerfa (ERFA, the IAU's SOFA) | under 2.5 µas |
| Body-fixed axes for the Sun, seven planets (excluding Earth), and Pluto | SPICE `sxform` with `pck00011` | under 10 µas |
| IAU orientation polynomials | NAIF's own constants | under 2 nas |
| Lunar frames | SPICE `sxform` with the DE440 lunar kernels, plus synchronous rotation | under 0.03 µas; Earth within 15° of +X |
| Ephemeris positions | Earth-Moon barycentre closure | 4.8e-10 relative |
| Ephemeris velocities | Central difference of position over ±10 s | 1.4e-6 relative |
