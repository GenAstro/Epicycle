# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
#
# Atmospheric drag — API force type plus the extension seams pluggable atmospheres and space-weather
# providers implement. Concrete atmospheres (e.g. Exponential, MSISE00) live in their own files.

# ─────────────────────────── atmosphere model seam ───────────────────────────
"""
    AbstractDensityModel

The atmosphere used by [`AtmosphericDrag`](@ref), chosen with its `model` keyword — it gives the air
density at the spacecraft.

# Available models
- [`Exponential`](@ref) — analytic exponential atmosphere. Open.
- `MSISE00` (NRLMSISE-00), `JB2008` (Jacchia-Bowman 2008), `JR1971` (Jacchia-Roberts 1971),
  `Jacchia1977`, `HarrisPriester` and `HarrisPriesterModified` — **Enterprise**. The Force Models
  guide gives each one's altitude range, space-weather inputs and validation.

# Writing your own
Define a type that subtypes `AbstractDensityModel` and give it a
`density(model, jd, r_eci, earth_axes)` method returning density in kg/m³. `AtmosphericDrag` then
works with it unchanged.

!!! note "Enterprise"
    The Enterprise models come from the `EpicycleEnterprise` package (commercial license). See the
    Force Models guide for details.
"""
abstract type AbstractDensityModel end

"""
    density(model::AbstractDensityModel, jd, r_eci, earth_axes) -> ρ   [kg/m³]

Air density at the spacecraft: `jd` the UTC Julian date, `r_eci` the ICRF position [km], and
`earth_axes` what takes it to the Earth-fixed frame. Each atmosphere provides this method;
`AtmosphericDrag` calls it, and the open-source version implements it for [`Exponential`](@ref).

`AtmosphericDrag` passes the ICRF-to-ITRF rotation, a 3×3 matrix, from the axes it chose at
construction (by default the frame theory's). A direct call may instead pass an EOP table, as
before, which takes SatelliteToolbox's FK5 route from mean J2000. A model reaches geodetic
coordinates through a helper that accepts either.
"""
function density end

# ────────────────────────────── atmospheric drag ─────────────────────────────
"""
    AtmosphericDrag(body; model = Exponential(), orientation = nothing)

Atmospheric drag on the spacecraft, from its velocity relative to a rigidly rotating atmosphere.

The cannonball drag acceleration is

```math
\\vec{a}_{\\mathrm{drag}} = -\\tfrac{1}{2}\\, \\rho\\, \\frac{C_d A}{m}\\,
                              |\\vec{v}_{\\mathrm{rel}}|\\, \\vec{v}_{\\mathrm{rel}},
\\qquad \\vec{v}_{\\mathrm{rel}} = \\vec{v} - \\vec{\\omega}_\\oplus \\times \\vec{r},
```

with ``\\rho`` from `model`, ``\\vec{\\omega}_\\oplus`` the Earth's rotation about its own pole,
taken from the same Earth-fixed axes the density is evaluated in (for the frame theories, the
sidereal rate corrected for the length of day), and ``\\vec{r}, \\vec{v}`` the spacecraft's ICRF
position and velocity.

# Arguments
- `body::CelestialBody`: the central body whose atmosphere acts on the spacecraft (positional,
  required — the atmosphere is that body's, so there is no default).
- `model::AbstractDensityModel`: which atmosphere gives the air density in kg/m³.
  Open: [`Exponential`](@ref). Enterprise: `MSISE00`, `JB2008`, `JR1971`, `Jacchia1977`,
  `HarrisPriester`, `HarrisPriesterModified`.
- `orientation::AbstractOrientationModel`: the Earth-fixed axes the atmosphere is evaluated in.
  Leave it out for the Earth's orientation model, the frame theory.

The central body must be the Earth: the density models and the geodetic altitude they take are
the Earth's.

# Notes
Reads the drag coefficient ``C_d`` and area ``A`` from the spacecraft's `SphericalDrag`, and
total mass from the spacecraft, at each integration step — set `sc.drag` before propagating.

The atmosphere is placed on the Earth through the axes chosen at construction: the `orientation`
keyword if given, otherwise `orientation_model(body)`, which for the Earth is the frame theory in
force then. Set the frame theory before building the force; changing it afterwards does not change
a force already built. This is the same rule `HarmonicGravity` follows.

# Examples
```julia
using AstroProp, AstroUniverse
drag = AtmosphericDrag(earth; model = Exponential())
```
"""
struct AtmosphericDrag{OT<:AbstractOrientationModel, DM<:AbstractDensityModel} <: OrbitODE
    central_body::CelestialBody
    model::DM
    orientation::OT
    dependencies::Vector{Type{<:AbstractVarTag}}
    num_funs::Int
end

function AtmosphericDrag(body::CelestialBody;
                        model::AbstractDensityModel = Exponential(),
                        orientation::Union{Nothing,AbstractOrientationModel} = nothing)
    # The density models and the geodetic conversion under them are the Earth's.
    body.naifid == 399 || throw(ArgumentError(
        "AtmosphericDrag models the Earth's atmosphere; got central body $(body.name). " *
        "Its density models and the geodetic altitude they take are the Earth's."))
    axes = orientation === nothing ? orientation_model(body) : orientation
    _check_axes(axes, body)
    return AtmosphericDrag(body, model, axes,
                           Type{<:AbstractVarTag}[PosVel], 6)
end

function Base.show(io::IO, ::MIME"text/plain", f::AtmosphericDrag)
    println(io, "AtmosphericDrag:")
    println(io, "  central_body = ", f.central_body.name)
    println(io, "  model        = ", nameof(typeof(f.model)), "()")
end
Base.show(io::IO, f::AtmosphericDrag) = show(io, MIME"text/plain"(), f)

function accel_eval!(force::AtmosphericDrag, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    geom = sc.drag
    geom === nothing && throw(ArgumentError(
        "sc.drag must be a SphericalDrag for AtmosphericDrag; got nothing. " *
        "Set sc.drag = SphericalDrag(; c_d, drag_area) before propagating."))
    jd  = force_epoch(params, t).utc
    r   = SVector{3}(x̄[1], x̄[2], x̄[3])
    v   = SVector{3}(x̄[4], x̄[5], x̄[6])
    R, Ṙ = _rotation_blocks(force_rotation(params, force.orientation, force.central_body.naifid, t))
    ρ   = density(force.model, jd, x̄, R)
    BC  = geom.c_d * geom.drag_area / total_mass(sc)         # Cd·A/m [m²/kg]
    # Velocity relative to the co-rotating atmosphere, from the same rotation the density used:
    # in the Earth-fixed axes it is R v + Ṙ r, and back in ICRF Rᵀ(R v + Ṙ r). That is
    # v − ω × r with ω the Earth's spin about its own pole, not the ICRF z axis.
    v_app = R' * (R * v + Ṙ * r)
    dfac  = -0.5 * BC * ρ * norm(v_app) * 1.0e3               # → km/s²
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = dfac * v_app[1]; x̄̇[5] = dfac * v_app[2]; x̄̇[6] = dfac * v_app[3]
    return x̄̇
end
