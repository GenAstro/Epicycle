# Force Models

Each Epicycle force model is an `OrbitODE` that contributes an acceleration, and `ForceModel` sums
any mix of them for the propagator. This page documents the full capability, open and Enterprise.
Enterprise features are marked and require the `EpicycleEnterprise` package under a commercial
licence; everything else is open.

Fidelity is chosen by a `model` tag on the force rather than by swapping the force type. Moving
from open to Enterprise fidelity changes that one tag and leaves the surrounding script unchanged.

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

!!! note "Enterprise — full-field geopotential"
    `EGM96` and `EGM2008` provide the full spherical-harmonic field at high degree and order,
    behind the same `AbstractGeopotential` extension interface (abstract type). They ship in the
    `EpicycleEnterprise` package under a commercial licence. The call is identical to the open
    path, and only the `model` tag changes:

    ```julia
    using EpicycleEnterprise
    HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())
    ```

    Without `EpicycleEnterprise` loaded, naming `EGM96()` raises an `UndefVarError`.

## Atmospheric drag

```@docs
AtmosphericDrag
AbstractDensityModel
Exponential
```

!!! note "Enterprise — NRLMSISE-00 density"
    `MSISE00` is the NRLMSISE-00 empirical density model, reaching `AtmosphericDrag` through the
    `AbstractDensityModel` extension interface (abstract type). It ships in `EpicycleEnterprise`
    under a commercial licence:

    ```julia
    using EpicycleEnterprise
    AtmosphericDrag(earth; model = MSISE00(),
                     space_weather = ConstantSpaceWeather(f107 = 150, f107a = 150, magnetic_index = 3))
    ```

## Solar radiation pressure

```@docs
SolarRadiationPressure
```

Shadow (eclipse) model, open: `DualCone`, which models both umbra and penumbra. It is the
default and the only one.

## Spacecraft geometry

Drag and SRP read their geometry from the spacecraft, not from the force — switching cannonball to a
higher-fidelity geometry is a type change on the spacecraft field, leaving the force untouched.

```@docs
SphericalDrag
SphericalSRP
```
