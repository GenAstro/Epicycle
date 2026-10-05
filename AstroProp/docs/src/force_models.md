# Force Models

AstroProp combines force accelerations in a `ForceModel` for orbit propagation. Forces include
gravity, atmospheric drag and solar radiation pressure (SRP). Enterprise forces and models require
the commercially licensed `EpicycleEnterprise` package and are identified below.

A force's underlying model sets its fidelity. Moving from the open J2–J5 field to EGM96, or from
the exponential atmosphere to NRLMSISE-00, changes the `model` argument of the same force type.

## Composing a force model

`ForceModel` combines forces about a common central body. All forces must name the same body.
Exactly one force supplies central gravity: `PointMassGravity` includes it by default; when
`HarmonicGravity` supplies it instead, use `include_center = false` on the point-mass force to
add only third-body gravity. `ForceModel` raises an error if central bodies differ or central
gravity is counted twice.

Drag and SRP read geometry, coefficients and mass from the spacecraft. These properties are
configured under [Spacecraft geometry](@ref force-models-spacecraft-geometry).

The open package supports point-mass gravity and the Earth's J2–J5 zonal field:

```julia
using AstroProp, AstroUniverse

ForceModel(PointMassGravity(earth))                              # two-body
ForceModel(PointMassGravity(earth, (moon, sun)))                 # with the Moon and Sun
ForceModel(PointMassGravity(sun, (earth, mars, jupiter)))        # Sun-centred, interplanetary
ForceModel(HarmonicGravity(earth; degree = 5, order = 0, model = Zonal()),
           PointMassGravity(earth, (moon, sun); include_center = false))   # J2–J5 and third bodies
```

!!! note "Enterprise Feature"
    This example requires `EpicycleEnterprise` for EGM96, solid tides, NRLMSISE-00 and relativity.

An Earth configuration with harmonic gravity, solid tides, third-body gravity, drag, SRP and
relativity uses Enterprise models. This example selects degree and order 70 and uses the default
tide, SRP and relativity settings:

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse, EpicycleEnterprise

forces = ForceModel(
    HarmonicGravity(earth; degree = 70, order = 70, model = EGM96(), tides = SolidTides()),
    PointMassGravity(earth, (moon, sun); include_center = false),
    AtmosphericDrag(earth; model = MSISE00()),
    SolarRadiationPressure(earth),
    Relativity(earth),
)
```

The following sections describe orientation, model choices and limits.
[Additional body configurations](@ref) gives examples for other central bodies.

## Body orientation

The body's orientation model supplies the fixed axes used by gravity, atmospheric drag, tides
and relativity. Forces read the orientation when they are built, so configure it before
constructing the force model:

```@raw html
<!-- doc-fragment -->
```
```julia
set_frame_theory!(FK5())                       # the Earth: IAU2006() (default) or FK5()
set_orientation!(bennu, MySpinModel(...))      # a body Epicycle ships no model for
```

The Earth's orientation is its frame theory. `IAU2006()` and `FK5()` both follow the IERS
conventions, including the observed celestial-pole offsets, and give the same Earth-fixed axes to
within 0.25 mas. The Moon defaults to its principal axes, `LunarPA()`, and the Sun, the planets
and Pluto to `IAU2015()`. A gravity field can specify its own coefficient axes, as described under
[Harmonic gravity](@ref).

## Point-mass gravity

`PointMassGravity` provides central-body gravity and gravity from specified third bodies.
Set `include_center = false` when harmonic gravity already supplies the central contribution.
Gravitational parameters come from the bodies, such as `earth.mu` and `moon.mu`; third-body
positions come from their ephemerides, with DE440 used for the built-in planetary ephemeris.
Each third-body contribution includes the direct acceleration of the spacecraft and subtracts
the acceleration of the central body.

The full signature is in the reference: [`PointMassGravity`](@ref).

## Harmonic gravity

`HarmonicGravity` evaluates a gravity field at the selected degree and order. The open package
provides `Zonal()`, the Earth's J2–J5 zonal field. Enterprise provides fields including:

!!! note "Enterprise Feature"
    The full gravity fields listed below, including ICGEM file support, require the commercially licensed `EpicycleEnterprise` package. `Zonal()` is available in the open package.

| Model | Body | Degree and order | Axes | Tide system |
|:---|:---|:---|:---|:---|
| `EGM96()` | Earth | 360 | the frame theory | tide-free |
| `EGM2008()` | Earth | 2190 × 2159 | the frame theory | tide-free |
| `GL0660B()` | Moon | 660 | lunar principal axes, `LunarPA()` | tide-free |
| `JGM85F01()` | Mars | 85 | IAU 1991 Mars axes, `IAU1991()` | tide-free |
| `IcgemGravity(file)` | any | the file's | the body's orientation model, or as given | from the file's header, or as given |

A field brings its own gravitational parameter and reference radius, its axes and its tide system,
so `HarmonicGravity` needs only the body, the degree and order, and the field:

```@raw html
<!-- doc-fragment -->
```
```julia
using EpicycleEnterprise
HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())
HarmonicGravity(mars;  degree = 50, order = 50, model = JGM85F01())
HarmonicGravity(moon;  degree = 60, order = 60, model = GL0660B())
```

An ICGEM file states its gravitational parameter, radius and tide system in its header. Where the
file's coefficients are defined in axes other than the body's default, or its header lacks the tide
system, give them when the field is built:

```@raw html
<!-- doc-fragment -->
```
```julia
ggm2b = IcgemGravity("https://icgem.gfz-potsdam.de/getmodel/gfc/.../ggm2bc80.gfc";
                     orientation = IAU1991())
HarmonicGravity(mars; degree = 80, order = 80, model = ggm2b)
```

`IcgemGravity` refuses a file whose gravitational parameter differs from the central body's by
more than 10 %, which is what a file for the wrong body looks like. Coefficients that vary with
time are evaluated at the propagation epoch. Degree and order have no default, and are checked
against what the field supports.

The full signatures are in the reference: [`HarmonicGravity`](@ref), [`Zonal`](@ref),
[`AbstractGeopotential`](@ref) and [`AstroProp.field_orientation`](@ref), which a new field
defines to declare its axes.

### Solid tides

!!! note "Enterprise Feature"
    Solid-tide corrections require the commercially licensed `EpicycleEnterprise` package.

`SolidTides` adds solid-tide corrections to a `HarmonicGravity` field through its `tides` keyword.
It requires Enterprise:

```@raw html
<!-- doc-fragment -->
```
```julia
HarmonicGravity(earth; degree = 70, order = 70, model = EGM96(), tides = SolidTides())
```

`SolidTides(; raising, k2, k3, degree3 = true)` specifies the raising bodies and Love numbers:

| Input | Default | Meaning |
|:---|:---|:---|
| `raising` | `(moon, sun)` for the Earth | the bodies raising the tide |
| `k2`, `k3` | 0.30190 and 0.093 for the Earth, IERS Conventions (2010) Table 6.3 | the Love numbers of degree 2 and 3 |
| `degree3` | `true` | include the degree-3 term |

The model is IERS Conventions (2010) §6.2 step 1, with one Love number for every order of a degree.
It leaves out the frequency-dependent step-2 corrections and the pole tide, a few percent of the
tidal acceleration in low Earth orbit.

The gravity field supplies the reference radius, spin axis and tide system. For a tide-free field,
such as EGM96 or EGM2008, the correction includes the permanent tide. For a zero-tide field,
it removes that contribution from the correction because the field already contains it.
Mean-tide fields are rejected.

The default Love numbers are the Earth's and apply only to the Earth. Tides on another body need
its own:

```@raw html
<!-- doc-fragment -->
```
```julia
HarmonicGravity(moon; degree = 60, order = 60, model = GL0660B(),
                tides = SolidTides(raising = (earth, sun), k2 = 0.02405, k3 = 0.0089))
```

## Atmospheric drag

`AtmosphericDrag` computes drag using a density model and the spacecraft's drag geometry,
coefficient and mass. The open package provides `Exponential`. Each density model describes one
body's atmosphere; `AtmosphericDrag` rejects a model for a different body. The atmosphere rotates
with the body, so density, altitude and wind use the body's orientation.

Enterprise provides Earth atmosphere models including:

!!! note "Enterprise Feature"
    The atmosphere models listed below require the commercially licensed `EpicycleEnterprise` package. `Exponential` is available in the open package.

| Model | Atmosphere | Altitude | Space weather | Checked against |
|:---|:---|:---|:---|:---|
| `MSISE00()` | NRLMSISE-00 | surface to 1000 km | F10.7, Ap | Orekit, GMAT |
| `JB2008()` | Jacchia-Bowman 2008 | 90 to 3000 km | F10.7, S10, M10, Y10, Dst | Orekit |
| `JR1971()` | Jacchia-Roberts 1971 | 90 to 3000 km | F10.7, Kp | GMAT's published densities |
| `Jacchia1977()` | Jacchia 1977 (SAO Special Report 375) | 90 to 2000 km | F10.7, Kp | — |
| `HarrisPriester(; n = 4)` | Harris-Priester, mean solar activity | 100 to 1000 km | none | Orekit |
| `HarrisPriesterModified(; n = 4)` | Harris-Priester, smooth, scaled by F10.7 (Hatten and Russell 2017) | 100 to 1000 km | F10.7 81-day average | — |

The Harris-Priester models take the cosine exponent of the diurnal bulge as `n`, from 2 to 6.
Below a model's lower limit the propagation stops with an error naming the limit. Above its upper
limit the density is zero; Orekit, by contrast, evaluates NRLMSISE-00 and JB2008 above theirs.
`HarrisPriesterModified` is smooth in its first derivatives, which suits an optimizer or a state
transition matrix.

Space weather is read from the published index tables, downloaded and cached on first use. An epoch
outside them is an error naming the span available.

```@raw html
<!-- doc-fragment -->
```
```julia
using EpicycleEnterprise
AtmosphericDrag(earth; model = MSISE00())
AtmosphericDrag(earth; model = HarrisPriesterModified(n = 6))    # a polar orbit
```

The full signatures are in the reference: [`AtmosphericDrag`](@ref) and [`Exponential`](@ref),
with [`AbstractDensityModel`](@ref) and [`atmosphere_body`](@ref), the interface a new atmosphere
implements.

## Solar radiation pressure

`SolarRadiationPressure` computes SRP using the spacecraft's SRP geometry, coefficient and mass,
with eclipse attenuation from a shadow model. Only the central body's shadow is included. The Sun
and occulting body are treated as spheres with their equatorial radii; the Sun's radius is
`sun.equatorial_radius`.

Supported shadow models include `DualCone`, the default, and `SmoothedConical`. `DualCone` computes
umbra and penumbra from the overlap of the apparent disks under these spherical assumptions.
`SmoothedConical` replaces the penumbra transition with a smooth curve. At its default sharpness,
its lighting factor differs from `DualCone` by at most 0.04 at LEO and GEO. Use it when SRP feeds
a gradient, since its acceleration has continuous derivatives through eclipse entry and exit.

`DualCone` has a kink at each penumbra boundary. A step across a boundary can introduce an error
that the integrator's error estimate misses, so `propagate!` predicts the boundaries and ends a
step at each one. `SmoothedConical` needs no boundary steps. [Verification](@ref) summarizes
the measured effect of boundary handling.

The full signatures are in the reference: [`SolarRadiationPressure`](@ref), [`DualCone`](@ref)
and [`SmoothedConical`](@ref).

## Relativity

!!! note "Enterprise Feature"
    Relativistic corrections require the commercially licensed `EpicycleEnterprise` package.

`Relativity` provides the post-Newtonian correction of IERS Conventions (2010) eq. 10.12 and
requires Enterprise. For a body other than Earth, supply its angular momentum in km²/s or disable
the Lense-Thirring term:

```@raw html
<!-- doc-fragment -->
```
```julia
Relativity(earth)
Relativity(mars; angular_momentum = 297.0)        # km²/s, (C/MR²) R² ω
Relativity(mars; lense_thirring = false)
```

`Relativity(body; schwarzschild = true, lense_thirring = true, de_sitter = true, angular_momentum)`
has one switch per term, all on by default:

| Term | What it is |
|:---|:---|
| Schwarzschild | the central body's mass, the largest term |
| Lense-Thirring | frame dragging by the body's spin |
| de Sitter | geodesic precession from the body's motion about the Sun; left out when the Sun is the central body |

The Lense-Thirring term needs the body's spin angular momentum per unit mass, `angular_momentum`
in km²/s. For Earth it defaults to the IERS value, 980 km²/s; for any other body it must be given
or the term turned off. `Relativity` raises an error otherwise. The spin axis is the z axis of
the body's orientation.

Angular momentum per unit mass is (C/MR²) R² ω. The examples below use 297 km²/s for Mars with
C/MR² = 0.3644, and 3.16 km²/s for the Moon with C/MR² = 0.3929. The same formula gives Earth's
980 km²/s. Cross-tool comparisons must use the same angular momentum; see [Verification](@ref).

## [Spacecraft geometry](@id force-models-spacecraft-geometry)

The spacecraft's `drag` and `srp` fields supply the geometry and coefficients used by drag and
SRP. Both forces use the spacecraft's total mass. Changing a geometry model changes the
corresponding spacecraft field:

```@raw html
<!-- doc-fragment -->
```
```julia
sat = Spacecraft(; state = ..., time = ..., mass = 1000.0,
                 drag = SphericalDrag(c_d = 2.2, drag_area = 10.0),
                 srp  = SphericalSRP(c_r = 1.8, srp_area = 10.0))
```

The full signatures are in the reference: [`SphericalDrag`](@ref) and [`SphericalSRP`](@ref).

## Additional body configurations

These examples configure forces for Mars, the Moon, the Sun and Venus, then define bodies and
orientations for Phobos and Bennu. Body-specific inputs and missing models are noted beside each
example. Angular-momentum values are described under [Relativity](@ref).

### Mars

!!! note "Enterprise Feature"
    This configuration requires `EpicycleEnterprise` for JGM85F01, solid tides and relativity.

This Mars configuration includes harmonic gravity, tides, third-body gravity, SRP and relativity.
Supply Mars's Love numbers and angular momentum explicitly. The example uses k₂ = 0.169
(Konopliv et al. 2016) and disables degree 3 because k₃ is not well determined. No Mars atmosphere
model ships, so drag is omitted.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse, EpicycleEnterprise

forces = ForceModel(
    HarmonicGravity(mars; degree = 85, order = 85, model = JGM85F01(),
                    tides = SolidTides(raising = (sun,), k2 = 0.169, degree3 = false)),
    PointMassGravity(mars, (sun, earth, jupiter); include_center = false),
    SolarRadiationPressure(mars; shadow = DualCone(), solar_flux = 1367.0,
                           nominal_sun = 149597870.691),
    Relativity(mars; schwarzschild = true, lense_thirring = true, de_sitter = true,
               angular_momentum = 297.0),
)
```

### Moon

!!! note "Enterprise Feature"
    This configuration requires `EpicycleEnterprise` for GL0660B, solid tides and relativity.

This lunar configuration uses GRAIL Love numbers k₂ = 0.02405 and k₃ = 0.0089
(Konopliv et al. 2013) and explicit lunar angular momentum. No lunar atmosphere model ships,
so drag is omitted. SRP includes only the Moon's shadow; the Earth's shadow on a lunar orbiter
is not modelled.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse, EpicycleEnterprise

forces = ForceModel(
    HarmonicGravity(moon; degree = 100, order = 100, model = GL0660B(),
                    tides = SolidTides(raising = (earth, sun), k2 = 0.02405, k3 = 0.0089,
                                       degree3 = true)),
    PointMassGravity(moon, (earth, sun); include_center = false),
    SolarRadiationPressure(moon; shadow = DualCone(), solar_flux = 1367.0,
                           nominal_sun = 149597870.691),
    Relativity(moon; schwarzschild = true, lense_thirring = true, de_sitter = true,
               angular_momentum = 3.16),
)
```

### Sun

!!! note "Enterprise Feature"
    This configuration requires `EpicycleEnterprise` for relativity. Point-mass gravity and SRP are available in the open package.

The planets provide third-body gravity in this Sun-centred configuration. SRP remains in full
sunlight because the Sun does not occult its own light. `Relativity` omits the de Sitter term
about the Sun. This example disables Lense-Thirring rather than supplying solar angular momentum.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse, EpicycleEnterprise

forces = ForceModel(
    PointMassGravity(sun, (mercury, venus, earth, mars, jupiter, saturn, uranus, neptune)),
    SolarRadiationPressure(sun),
    Relativity(sun; lense_thirring = false),
)
```

### Venus with an ICGEM field

!!! note "Enterprise Feature"
    This configuration requires `EpicycleEnterprise` for ICGEM gravity fields and relativity.

This Venus configuration loads a gravity field from an ICGEM file. Its coefficients are assumed
to use the body's default IAU 2015 axes unless an orientation is supplied with the field.
Lense-Thirring is disabled because the example does not supply Venus's angular momentum.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse, EpicycleEnterprise

forces = ForceModel(
    HarmonicGravity(venus; degree = 60, order = 60, model = IcgemGravity("venus_field.gfc")),
    PointMassGravity(venus, (sun, earth); include_center = false),
    SolarRadiationPressure(venus),
    Relativity(venus; lense_thirring = false),
)
```

Without a field file, `PointMassGravity(venus, (sun, earth))` supplies Venus's central gravity
and the Sun's and Earth's third-body contributions.

### Phobos with a SPICE orientation

Epicycle provides no built-in Phobos body. This example defines it from its gravitational
parameter, radius and NAIF ID, loads its ephemeris from a Mars satellite kernel, and sets its
orientation from the IAU frame in a planetary constants kernel.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse

download_spice_kernel("mar097.bsp",
    "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar097.bsp")
load_spice_kernel("mar097.bsp")
load_spice_kernel("pck00011.tpc")

phobos = CelestialBody("Phobos", 7.11e-4, 11.08, 0.0, 401)     # μ km³/s², radius km, flattening, NAIF ID
set_orientation!(phobos, SpiceOrientation("IAU_PHOBOS"))

forces = ForceModel(
    PointMassGravity(phobos, (mars, sun)),
    SolarRadiationPressure(phobos),
)
```

### Bennu with a custom spin model

This example defines Bennu and its spin model with `pole_axes_rotation`. It requires a local
`bennu.bsp` ephemeris kernel, obtained from the mission or JPL Horizons. The spin model is set
before the forces are built.

```@raw html
<!-- doc-fragment -->
```
```julia
using AstroProp, AstroModels, AstroUniverse

load_spice_kernel("bennu.bsp")                                  # the asteroid's ephemeris

struct SimpleSpin <: AbstractOrientationModel
    pole_ra::Float64       # deg
    pole_dec::Float64      # deg
    pm0::Float64           # prime meridian at J2000, deg
    spin_rate::Float64     # deg/day
end

function AstroUniverse.body_axes_rotation(m::SimpleSpin, naifid, jd_tdb)
    d = jd_tdb - 2451545.0
    Ẇ = deg2rad(m.spin_rate)
    return pole_axes_rotation(deg2rad(m.pole_ra), deg2rad(m.pole_dec),
                              deg2rad(m.pm0) + Ẇ * d, 0.0, 0.0, Ẇ / 86_400)
end

bennu = CelestialBody("Bennu", 4.892e-9, 0.2825, 0.0, 2101955)
set_orientation!(bennu, SimpleSpin(85.46, -60.36, 89.6, 2011.145))

forces = ForceModel(
    PointMassGravity(bennu, (sun, earth, jupiter)),
    SolarRadiationPressure(bennu),
)
```

The next section describes the interface for adding a custom acceleration to these force models.

## Writing your own force

To add a custom force, define a subtype of `OrbitODE` and implement `accel_eval!`. `ForceModel`
sums your acceleration with the built-in forces.

Write acceleration in km/s² to rows 4 to 6 of `dy`. Each force receives a zero-filled `dy`, and
its contribution is summed afterwards, so assignment and addition give the same result. Leave
rows 1 to 3 to the propagator, and keep the method's element types open so automatic
differentiation can compute the Jacobian. The example adds a constant acceleration along the
velocity, such as a low-thrust engine held prograde, to point-mass gravity.

```julia
using AstroEpochs, AstroStates, AstroUniverse, AstroModels, AstroProp
using LinearAlgebra: norm
import AstroProp: accel_eval!

# A constant acceleration along the velocity vector, in km/s².
struct AlongTrackThrust <: OrbitODE
    accel::Float64
end

# Rows 4 to 6 only. Leave the element types open so the Jacobian can be taken through it.
function accel_eval!(f::AlongTrackThrust, t, y, dy, sc, params)
    v = y[4:6]
    dy[4:6] .= f.accel .* v ./ norm(v)
    return dy
end

sat = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0]),
                 time  = Time("2020-01-01T00:00:00", UTC(), ISOT()))
forces = ForceModel(PointMassGravity(earth, ()), AlongTrackThrust(1.0e-7))
prop   = OrbitPropagator(forces, IntegratorConfig(Tsit5(); dt = 60.0, reltol = 1e-10, abstol = 1e-10))

# A day of thrusting raises the orbit.
propagate!(prop, sat, StopAt(sat, PropDurationDays(), 1.0))
```

If your force needs time conversions, body rotations or ephemeris states, use
`AstroProp.force_epoch`, `force_rotation`, `force_position` and `force_state` through `params`.
These helpers share cached results across forces in an evaluation. When you call `accel_eval!`
outside a propagation, they compute the results directly.

```@raw html
<!-- doc-fragment -->
```
```julia
# Inside accel_eval!(f, t, y, dy, sc, params):
jd_utc = AstroProp.force_epoch(params, t).utc
R      = AstroProp.force_rotation(params, orientation_model(earth), 399, t)   # 6×6, [R 0; Ṙ R]
r_sun  = AstroProp.force_position(params, earth, sun, t)                      # km, ICRF
```

## Verification

Force verification includes comparisons with GMAT and Orekit on Earth orbits from LEO to Molniya,
and combined-model comparisons about the Sun, Moon and Mars with Orekit and about Phobos and
Bennu with GMAT. The atmosphere table identifies the comparison source for each density model;
a dash indicates that no comparison is listed. The
[Verification](https://genastro.github.io/Epicycle/Epicycle/dev/verification/) page gives the cases,
configurations and measured agreement.

Cross-tool comparisons also depend on orientation and physical parameters. GMAT's Earth-fixed
frame follows the FK5 chain but omits celestial-pole offsets; in 2020 its pole differed from the
IERS pole by 47 mas. GMAT also approximates angular momentum per unit mass as 2/5 R²ω,
assuming a uniform sphere. This gives 1186.6 km²/s for Earth, compared with the IERS value of
980 km²/s used by `Relativity`. The verification page documents the effect of the orientation
difference and the agreement after it is removed.

Eclipse-boundary handling was checked over one day in LEO with the full force model and Vern9
at a tolerance of 1e-11. Ending steps at the penumbra boundaries reduced the difference
from the converged solution from 15 cm to 2.2 mm, with 1 % fewer force evaluations.
