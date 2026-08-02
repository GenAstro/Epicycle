# Force Models

Epicycle force models are composable: each is an `OrbitODE` that contributes an acceleration, and
`ForceModel` sums any mix of them for the propagator. This page documents the **full** capability —
open and Enterprise. Enterprise features are marked with an admonition and require the
`EpicycleEnterprise` package (commercial license); everything else is open.

Fidelity is chosen by a *model* tag on the force, not by swapping the force type. Moving from open to
Enterprise fidelity changes one tag — the surrounding script is unchanged. (See the model-seam
pattern in the architecture spec, §11.10.)

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
    `EGM96` and `EGM2008` provide the full spherical-harmonic field (high degree and order) through
    the same `AbstractGeopotential` seam. They ship in the `EpicycleEnterprise` package (commercial
    license). The call is identical to the open path — only the `model` tag changes:

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
    `MSISE00` is the NRLMSISE-00 empirical density model, plugged into `AtmosphericDrag` through the
    `AbstractDensityModel` seam. It ships in `EpicycleEnterprise` (commercial license):

    ```julia
    using EpicycleEnterprise
    AtmosphericDrag(earth; model = MSISE00())
    ```

## Solar radiation pressure

```@docs
SolarRadiationPressure
```

Shadow (eclipse) models, all open: `NoShadow`, `Cylindrical`, and `DualCone` — the dual-cone
(umbra + penumbra) model.

## Spacecraft geometry

Drag and SRP read their geometry from the spacecraft, not from the force — switching from a spherical
to a higher-fidelity geometry is a type change on the spacecraft field, leaving the force untouched.

```@docs
SphericalDrag
SphericalSRP
```
