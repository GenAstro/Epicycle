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
`density(model, jd, r_eci, body_axes)` method returning density in kg/m³. A model of an atmosphere
other than the Earth's also gives [`atmosphere_body`](@ref), the NAIF ID of its body.
`AtmosphericDrag` then works with it unchanged.

!!! note "Enterprise"
    The Enterprise models come from the `EpicycleEnterprise` package (commercial license). See the
    Force Models guide for details.
"""
abstract type AbstractDensityModel end

"""
    density(model::AbstractDensityModel, jd, r_eci, body_axes) -> ρ   [kg/m³]

Air density at the spacecraft: `jd` the UTC Julian date, `r_eci` the ICRF position [km], and
`body_axes` what takes it to the body-fixed frame. Each atmosphere provides this method;
`AtmosphericDrag` calls it, and the open-source version implements it for [`Exponential`](@ref).

`AtmosphericDrag` passes the rotation from ICRF to the body's fixed axes, a 3×3 matrix, from the
body's orientation model; for the Earth that is the frame theory's ITRF. The model takes the
altitude from it in its own way, the Earth's models from the WGS84 ellipsoid. A direct call to an
Earth model may instead pass an EOP table, which takes SatelliteToolbox's FK5 route from mean
J2000. A model reaches geodetic coordinates through a helper that accepts either.
"""
function density end

"""
    atmosphere_body(model::AbstractDensityModel) -> Int

The NAIF ID of the body whose atmosphere `model` describes. [`AtmosphericDrag`](@ref) refuses a
model whose body is not its own. Every shipped model is the Earth's, which is the default; a model
of another body's atmosphere defines this method.

# Returns
- `Int`: the NAIF ID, 399 for every shipped model.

# Example
```julia
struct MarsExponential <: AbstractDensityModel end
AstroProp.atmosphere_body(::MarsExponential) = 499
```
"""
atmosphere_body(::AbstractDensityModel) = 399

# ────────────────────────────── atmospheric drag ─────────────────────────────
"""
    AtmosphericDrag(body; model = Exponential())

Atmospheric drag on the spacecraft, from its velocity relative to an atmosphere turning with its
body.

The cannonball drag acceleration is

```math
\\vec{a}_{\\mathrm{drag}} = -\\tfrac{1}{2}\\, \\rho\\, \\frac{C_d A}{m}\\,
                              |\\vec{v}_{\\mathrm{rel}}|\\, \\vec{v}_{\\mathrm{rel}},
\\qquad \\vec{v}_{\\mathrm{rel}} = \\vec{v} - \\vec{\\omega} \\times \\vec{r},
```

with ``\\rho`` from `model`, ``\\vec{\\omega}`` the body's rotation from its orientation model (for
the Earth, the frame theory's: the sidereal rate corrected for the length of day, about the
Earth's own pole), and ``\\vec{r}, \\vec{v}`` the spacecraft's ICRF position and velocity relative
to the body.

# Arguments
- `body::CelestialBody`: the central body whose atmosphere acts on the spacecraft (positional,
  required — the atmosphere is that body's, so there is no default).
- `model::AbstractDensityModel`: which atmosphere gives the air density in kg/m³. It must be
  `body`'s atmosphere, [`atmosphere_body`](@ref). Open: [`Exponential`](@ref). Enterprise:
  `MSISE00`, `JB2008`, `JR1971`, `Jacchia1977`, `HarrisPriester`, `HarrisPriesterModified`. All of
  these are the Earth's.

# Fields
- `central_body::CelestialBody`: the body whose atmosphere acts, `body`.
- `model::AbstractDensityModel`: the density model, `model`.
- `orientation::AbstractOrientationModel`: the axes the atmosphere turns with,
  `orientation_model(body)` when the force was built.

# Notes
Reads the drag coefficient ``C_d`` and area ``A`` from the spacecraft's `SphericalDrag`, and
total mass from the spacecraft, at each integration step — set `sc.drag` before propagating.

The atmosphere turns with the body in the body's own axes, `orientation_model(body)`, read when the
force is built; for the Earth that is the frame theory in force then. Set the frame theory before
building the force; changing it afterwards does not change a force already built. The density,
the altitude and the wind all come from that one rotation, the one the body's gravity uses.

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

function AtmosphericDrag(body::CelestialBody; model::AbstractDensityModel = Exponential())
    n = atmosphere_body(model)
    n == body.naifid || throw(ArgumentError(
        "$(nameof(typeof(model)))() is the atmosphere of NAIF $n" *
        (n == 399 ? " (the Earth)" : "") * ", not of $(body.name) (NAIF $(body.naifid)). " *
        "Choose a model of $(body.name)'s atmosphere; a model of another body's defines " *
        "`AstroProp.atmosphere_body`."))
    return AtmosphericDrag(body, model, orientation_model(body),
                           Type{<:AbstractVarTag}[PosVel], 6)
end

function Base.show(io::IO, ::MIME"text/plain", f::AtmosphericDrag)
    println(io, "AtmosphericDrag:")
    println(io, "  central_body = ", f.central_body.name)
    println(io, "  model        = ", nameof(typeof(f.model)), "()")
end
Base.show(io::IO, f::AtmosphericDrag) = show(io, MIME"text/plain"(), f)

function accel_eval!(force::AtmosphericDrag, t::Time, x̄::AbstractVector, x̄̇::AbstractVector,
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
    # in the body-fixed axes it is R v + Ṙ r, and back in ICRF Rᵀ(R v + Ṙ r). That is
    # v − ω × r with ω the body's spin, which for the Earth is about its own pole.
    v_app = R' * (R * v + Ṙ * r)
    dfac  = -0.5 * BC * ρ * norm(v_app) * 1.0e3               # → km/s²
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = dfac * v_app[1]; x̄̇[5] = dfac * v_app[2]; x̄̇[6] = dfac * v_app[3]
    return x̄̇
end
