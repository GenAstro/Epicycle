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

## Atmospheric drag

```@docs
AtmosphericDrag
AbstractDensityModel
Exponential
```

**Enterprise.** `MSISE00` is the NRLMSISE-00 empirical density model, reached through the same
`AbstractDensityModel` extension interface. It ships in `EpicycleEnterprise`, which is licensed
commercially:

```@raw html
<!-- doc-fragment -->
```
```julia
using EpicycleEnterprise
AtmosphericDrag(earth; model = MSISE00())
```

## Solar radiation pressure

```@docs
SolarRadiationPressure
```

One shadow model ships and is the default: `DualCone`, which models both umbra and penumbra.

## Spacecraft geometry

Drag and SRP read their geometry from the spacecraft rather than from the force, so moving from a
spherical geometry to a higher-fidelity one changes the type of a spacecraft field and leaves the
force untouched.

```@docs
SphericalDrag
SphericalSRP
```
