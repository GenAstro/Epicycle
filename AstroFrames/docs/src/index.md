```@meta
CurrentModule = AstroFrames
```

# AstroFrames

AstroFrames defines coordinate systems and transforms states between them. A coordinate system pairs an origin with an axes orientation. The package includes IAU 2006 and FK5 Earth frames, lunar frames, body-fixed frames for bodies with IAU orientation data, and orbit-relative frames. Transformations accept either a `Coordinate`, which holds a state, coordinate system, and epoch, or any object that provides the same information, including a `Spacecraft`.

AstroFrames is tested against [pyerfa](https://github.com/liberfa/pyerfa) for the Earth systems, and against SPICE for planetary body-fixed axes and lunar systems. pyerfa is the Python interface to ERFA, the open-licence release of the IAU's SOFA library.

AstroFrames leverages extensive work on Earth-based systems available in [SatelliteToolboxTransformations.jl](https://juliaspace.github.io/SatelliteToolboxTransformations.jl/stable/). We thank the authors of that package for their high-quality work.

The implementation follows these references:

- Chagas, R. A. J., *SatelliteToolboxTransformations.jl*, v1.2.0, JuliaSpace. https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
- Petit, G. and Luzum, B., eds. (2010), *IERS Conventions (2010)*, IERS Technical Note No. 36, Verlag des Bundesamts für Kartographie und Geodäsie.
- Vallado, D. A. (2013), *Fundamentals of Astrodynamics and Applications*, 4th ed., Microcosm Press / Springer.
- Acton, C. H. (1996), *Ancillary Data Services of NASA's Navigation and Ancillary Information Facility*, Planetary and Space Science, 44(1).

## Quick Start

The examples below illustrate how to create coordinate systems, convert a state, and evaluate an orbit property in another frame.

```julia
using AstroFrames, AstroUniverse, AstroEpochs, AstroModels, AstroStates, AstroCallbacks

# Create some coordinate systems providing the origin and axes
cs_j2000 = CoordinateSystem(earth, MJ2000Eq())
cs_fixed = CoordinateSystem(earth, ITRF())

# Define a coordinate: state + coordinate system + time
epoch = Time(2458849.5, 0.0, TDB(), JD())
c     = Coordinate([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0], cs_j2000, epoch)

# Convert the coordinate to the cs_fixed system
c_fixed = Coordinate(c, cs_fixed)

# Create a spacecraft and convert its state to the cs_fixed system
sc = Spacecraft()
Coordinate(sc, cs_fixed)

# Calculate an orbit property in ecliptic coordinates
cs_j2000Ec = CoordinateSystem(earth, MJ2000Ec())
raan(sc, cs_j2000Ec)
```

In the REPL, `subtypes(AbstractAxes)` lists every supported axes type, and `?MoonME`, for example, shows the help for one.

## Coordinate Systems

A `CoordinateSystem` pairs an origin with an axes. Some axes require a specific origin: Earth axes require Earth, lunar axes require the Moon, and body-fixed axes require the body encoded by the axes. The constructor raises an `ArgumentError` for an invalid pairing.

```julia
using AstroFrames, AstroUniverse, AstroModels

# Any body from AstroUniverse can be an origin.
cs_earth_fixed   = CoordinateSystem(earth, ITRF())
cs_moon_ecliptic = CoordinateSystem(moon, MJ2000Ec())

# A planet fixed system
cs_mars_fixed = CoordinateSystem(mars, CelestialBodyFixed())

# Moon fixed systems: mean Earth axes for surface coordinates, principal axes for gravity
cs_moon_me = CoordinateSystem(moon, MoonME())
cs_moon_pa = CoordinateSystem(moon, MoonPA())

# Spacecraft-based systems: a spacecraft is an origin, and orbit-relative axes take their
# reference orbit from it
chief   = Spacecraft(name = "chief")
cs_ric  = CoordinateSystem(chief, RIC())
cs_vnb  = CoordinateSystem(chief, VNB())
cs_lvlh = CoordinateSystem(chief, LVLH())

# The common pairings already have names.
EarthMJ2000Eq
MoonFixed

# An illegal pairing raises here, naming what the axes are for.
# CoordinateSystem(mars, GCRF())   # GCRF axes require an Earth origin
```

Earth and the Moon do not use `CelestialBodyFixed`. Their higher-fidelity frames are `ITRF`, `MoonME`, and `MoonPA`; requesting `CelestialBodyFixed` for either body raises an error.

## Coordinates and Conversions

A `Coordinate` holds a state, its coordinate system, and its epoch, so a conversion only needs the target coordinate system. The same conversion works with a `Spacecraft` or another type that contains state, coordinate system, and time.  

The examples below show how to convert a `Coordinate` to a new frame, convert a `Spacecraft' state to a new frame, and how to request the rotation matrix between systems. 

```julia
using AstroFrames, AstroUniverse, AstroEpochs, AstroModels, AstroStates, AstroCallbacks

# Define a coordinate providing time, state, and system
epoch = Time(2458849.5, 0.0, TDB(), JD())
system = CoordinateSystem(earth, MJ2000Eq())
state = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]
c = Coordinate(state,system, epoch)

# Convert the coordinate to Earth MJ2000Ec, this returns a new Coordinate
c_ecliptic = Coordinate(c, CoordinateSystem(earth, MJ2000Ec()))

# Convert the Coordinate to Moon centere MoonME axes
c_moon = Coordinate(c, CoordinateSystem(moon, MoonME()))

# Define a spacecraft centered Coordinate
sc = Spacecraft()
Coordinate(sc, CoordinateSystem(earth, ITRF()))

# Compute some quantities in different coordinates
position_vector(sc, CoordinateSystem(earth, ITRF()))
CartesianState(sc, CoordinateSystem(moon, MoonME()))

# Compute the rotation matrix from ICRF to ITRF 
M = axes_rotation(ICRF(), ITRF(), epoch)
```

`axes_rotation` returns a 6×6 transform in SPICE `sxform` block form. The upper-left and lower-right blocks contain the rotation, and the lower-left block contains its time derivative. For Earth frames, the rate includes sidereal rotation but treats precession, nutation, polar motion, and obliquity of date as constant during the transformation. 

## Earth Frame Theories

AstroFrames supports two Earth frame theories including FK5 IAU 1980, and IAU 2006/2010.  
IAU 2006/2010 uses `GCRF → CIRS → TIRS → ITRF` and we recommend that system for most work. FK5 / IAU 76-80 uses `MJ2000Eq → MODEq → TODEq → PEF → ITRF` and is supported for legacy system consistency.

```text
                    ICRF ─────────────── GCRF          IAU 2006/2000A
                     │      (identity)     │
             (frame bias)                 CIRS
                     │                     │
  FK5 / IAU 76-80  MJ2000Eq              TIRS
                     │                     │
                   MODEq ── MODEc         ITRF
                     │                     │
                   TODEq ── TODEc          └──── shared endpoint
                   │   │
                 PEF  TEME
                   │
                 ITRF
```

## Orbit-Relative Axes

There are several orbit relative frames supported such as RIC, LVLH, and VNB. If the origin carries a state, as a spacecraft does, AstroFrames uses that state as the reference. Otherwise, you must provide a state in the conversion interface as shown in the examples below. 

```julia
using AstroFrames, AstroUniverse, AstroEpochs, AstroModels, AstroStates

epoch = Time(2458849.5, 0.0, TDB(), JD())

# Create two spacecraft and compute the deputy's state in the chief's RIC frame
chief  = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]),
                    time = epoch, name = "chief")
deputy = Spacecraft(state = CartesianState([7000.1, 0.0, 0.0, 0.0, 7.546, 0.0]),
                    time = epoch, name = "deputy")
Coordinate(deputy, CoordinateSystem(chief, RIC()))

# Compute the rotation matrix between orbit relative systems by providing the 
# state explicity.
p = (; reference_state = [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0])
M = axes_rotation(ICRF(), RIC(), epoch, p)
```

## Lunar Frames

AstroFrames provides two lunar body-fixed frames. `MoonME` is the mean Earth/mean rotation axis frame used for published lunar coordinates, including surface latitude and longitude. `MoonPA` is aligned with the Moon's principal axes of inertia and is the frame used for lunar gravity coefficients. Both include the full physical, forced, and free libration. They differ by a fixed rotation of about 104 arcseconds, or about 875 metres at the surface.

```julia
using AstroFrames, AstroUniverse, AstroEpochs

epoch = Time(2458849.5, 0.0, TDB(), JD())

# Mapping: landing sites, surface features, latitude and longitude.
site = Coordinate([1737.4, 0.0, 0.0, 0.0, 0.0, 0.0],
                  CoordinateSystem(moon, MoonME()), epoch)

# Dynamics: gravity field evaluation, propagation about the Moon.
orb = Coordinate([1837.4, 0.0, 0.0, 0.0, 1.63, 0.0],
                 CoordinateSystem(moon, MoonPA()), epoch)
```

AstroFrames loads the lunar orientation kernels, leap seconds, and solar system ephemeris at startup. Using either lunar frame does not require a kernel path.

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

See the documentation for AstroUniverse for EOP data management.

## Origins Reference

An origin is any `AbstractPoint`. Celestial bodies from AstroUniverse serve as origins for body-fixed and body-centered inertial frames. A spacecraft, or another point that carries a state, can serve as the origin and reference orbit for orbit-relative axes.

AstroFrames provides named coordinate systems for common pairings: `EarthMJ2000Eq`, `EarthMJ2000Ec`, `EarthFixed`, `EarthTODEq`, `EarthICRF`, `MoonFixed`, `MoonPrincipalAxes`, and `SunMJ2000Ec`.

## Validation

AstroFrames is tested against external references where they are available.

| What | Checked against | Agreement |
|---|---|---|
| Earth chains, both theories | pyerfa (ERFA, the IAU's SOFA) | under 2.5 µas |
| Body-fixed axes, Sun and 9 planets | SPICE `sxform` with `pck00011` | under 10 µas |
| IAU orientation polynomials | NAIF's own constants | under 2 nas |
| Lunar frames | SPICE `sxform` with the DE440 lunar kernels, plus synchronous rotation | under 0.03 µas; Earth within 15° of +X |
| Ephemeris positions | Earth-Moon barycentre closure | 4.8e-10 relative |
| Ephemeris velocities | Central difference of position over ±10 s | 1.4e-6 relative |

## Deprecated Names

`ICRFAxes` and `MJ2000Axes` remain as deprecated aliases for [`ICRF`](@ref) and [`MJ2000Eq`](@ref). They warn on use.
