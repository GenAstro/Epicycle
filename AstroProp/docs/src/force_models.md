# Force Models

A force model is a set of forces, each an `OrbitODE` that contributes an acceleration, which
`ForceModel` sums for the propagator. This page documents the full capability, open and Enterprise.
Every Enterprise force is named as Enterprise below and needs the `EpicycleEnterprise` package,
which is licensed commercially; everything else is open.

Fidelity is chosen by the *model* a force is given rather than by swapping the force type, so moving
from open to Enterprise fidelity changes one argument and leaves the rest of the script unchanged.

## Composing forces

```julia
using AstroProp, AstroModels, AstroUniverse

forces = ForceModel(
    HarmonicGravity(earth; degree = 5, order = 0, model = Zonal()),
    AtmosphericDrag(earth; model = Exponential()),
    SolarRadiationPressure(earth; shadow = DualCone()),
)
```

## Central-body gravity

```@docs
PointMassGravity
HarmonicGravity
AbstractGeopotential
Zonal
```

**Enterprise.** `EGM96` and `EGM2008` provide the full spherical-harmonic field, to high degree
and order, through the same `AbstractGeopotential` extension interface. They ship in the
`EpicycleEnterprise` package, which is licensed commercially. The call is the one above with a
different model:

```@raw html
<!-- doc-fragment -->
```
```julia
using EpicycleEnterprise
HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())
```

Naming `EGM96()` without `EpicycleEnterprise` loaded raises an `UndefVarError`.

Enterprise also provides fields for other bodies: `GL0660B` for the Moon, `JGM85F01` for Mars,
and `IcgemGravity(url_or_path)` for a field in any ICGEM file.

**Which axes a field is evaluated in.** A field's coefficients describe the body's mass in the
body-fixed axes they were estimated in, so `HarmonicGravity` rotates into those axes, evaluates
there, and rotates back. It takes them from the `orientation` keyword if you give one, otherwise
from the field ([`AstroProp.field_orientation`](@ref)), otherwise from the body's orientation
model: the frame theory's ITRF for the Earth, `LunarPA()` for the Moon, `IAU2015()` for the
planets. `JGM85F01` declares `IAU1991()`, the Mars axes of the 1991 report it was estimated in.
The choice is made when the force is built, so set the frame theory first.

```julia
HarmonicGravity(mars; degree = 50, order = 50, model = JGM85F01())                # its IAU 1991 axes
HarmonicGravity(mars; degree = 80, order = 80, model = IcgemGravity(url),
                orientation = IAU1991())                                        # a file in those axes
```

Epicycle's Earth axes are the frame theory's ITRF chain from ICRF, with the IERS celestial-pole
corrections. GMAT reaches its Earth-fixed frame from FK5 mean J2000 without the corrections, which
over a day of LEO moves the state about 12 cm under `IAU2006()`. A comparison with GMAT passes
GMAT's axes with `orientation`; the AstroProp tests do this with a small orientation model.

```@docs
AstroProp.field_orientation
```

## Atmospheric drag

```@docs
AtmosphericDrag
AbstractDensityModel
Exponential
```

**Enterprise.** Six atmospheres of the Earth, reached through the same `AbstractDensityModel`
extension interface. They ship in `EpicycleEnterprise`, which is licensed commercially.

| Model | Atmosphere | Altitude | Space weather | Checked against |
|---|---|---|---|---|
| `MSISE00` | NRLMSISE-00 | surface to 1000 km | F10.7, Ap | Orekit, GMAT |
| `JB2008` | Jacchia-Bowman 2008 | 90 to 3000 km | F10.7, S10, M10, Y10, Dst | Orekit |
| `JR1971` | Jacchia-Roberts 1971 | 90 to 3000 km | F10.7, Kp | GMAT's published densities |
| `Jacchia1977` | Jacchia 1977 (SAO Special Report 375) | 90 to 2000 km | F10.7, Kp | — |
| `HarrisPriester` | Harris-Priester, mean solar activity | 100 to 1000 km | none | Orekit |
| `HarrisPriesterModified` | Harris-Priester, smooth, scaled by F10.7 (Hatten and Russell 2017) | 100 to 1000 km | F10.7 81-day average | — |

Space weather comes from `SpaceIndices` tables, loaded on first use; an epoch outside them is an
error naming the span available. Below a model's lower limit the propagation stops with an error
naming the limit, and above its upper limit the density is zero. The two Harris-Priester models
take the cosine exponent of the diurnal bulge as `n`, 4 by default. The atmosphere is the Earth's
in every case, so `AtmosphericDrag` refuses another central body.

```@raw html
<!-- doc-fragment -->
```
```julia
using EpicycleEnterprise
AtmosphericDrag(earth; model = MSISE00())
AtmosphericDrag(earth; model = HarrisPriesterModified(n = 6))    # a polar orbit
```

## Solar radiation pressure

```@docs
SolarRadiationPressure
```

Two shadow models ship. `DualCone`, the default, models both umbra and penumbra exactly.
`SmoothedConical` replaces the penumbra with a smooth curve that stays within 0.04 of `DualCone`, so
the SRP acceleration has continuous derivatives through eclipse entry and exit; use it where SRP
feeds a gradient.

```@docs
DualCone
SmoothedConical
```

## Spacecraft geometry

Drag and SRP read their geometry from the spacecraft rather than from the force, so moving from a
spherical geometry to a higher-fidelity one changes the type of a spacecraft field and leaves the
force untouched.

```@docs
SphericalDrag
SphericalSRP
```
