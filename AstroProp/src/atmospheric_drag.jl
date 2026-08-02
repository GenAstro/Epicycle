# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
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
- `MSISE00` — NRLMSISE-00 empirical atmosphere. **Enterprise**.

# Writing your own
Define a type that subtypes `AbstractDensityModel` and give it a `density(model, jd, r_eci, eop)`
method returning density in kg/m³. `AtmosphericDrag` then works with it unchanged.

!!! note "Enterprise"
    `MSISE00` comes from the `EpicycleEnterprise` package (commercial license). See the Force Models
    guide for details.
"""
abstract type AbstractDensityModel end

"""
    density(model::AbstractDensityModel, jd, r_eci, eop) -> ρ   [kg/m³]

Air density at the spacecraft: `jd` the UTC Julian date, `r_eci` the inertial position [km], `eop`
the Earth-orientation data. Each atmosphere provides this method; `AtmosphericDrag` calls it, and the
open-source version implements it for [`Exponential`](@ref).
"""
function density end

# ────────────────────────────── atmospheric drag ─────────────────────────────
"""
    AtmosphericDrag(body; model = Exponential())

Atmospheric drag on the spacecraft, from its velocity relative to a rigidly rotating atmosphere.

The cannonball drag acceleration is

```math
\\vec{a}_{\\mathrm{drag}} = -\\tfrac{1}{2}\\, \\rho\\, \\frac{C_d A}{m}\\,
                              |\\vec{v}_{\\mathrm{rel}}|\\, \\vec{v}_{\\mathrm{rel}},
\\qquad \\vec{v}_{\\mathrm{rel}} = \\vec{v} - \\vec{\\omega}_\\oplus \\times \\vec{r},
```

with ``\\rho`` from `model`, ``\\vec{\\omega}_\\oplus`` Earth's rotation rate (see
[`SatelliteToolboxBase.EARTH_ANGULAR_SPEED`]), and ``\\vec{r}, \\vec{v}`` the spacecraft's
inertial position and velocity.

# Arguments
- `body::CelestialBody`: the central body whose atmosphere acts on the spacecraft (positional,
  required — the atmosphere is that body's, so there is no default).
- `model::AbstractDensityModel`: which atmosphere gives the air density in kg/m³.
  Open: [`Exponential`](@ref). Enterprise: `MSISE00`.

# Notes
Reads the drag coefficient ``C_d`` and area ``A`` from the spacecraft's `SphericalDrag`, and
total mass from the spacecraft, at each integration step — set `sc.drag` before propagating.
Earth-orientation data is fetched once at construction via `fetch_iers_eop`.

# Examples
```julia
using AstroProp, AstroUniverse
drag = AtmosphericDrag(earth; model = Exponential())
```
"""
struct AtmosphericDrag{EoT, DM<:AbstractDensityModel} <: OrbitODE
    central_body::CelestialBody
    model::DM
    eop_data::EoT
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int
end

function AtmosphericDrag(body::CelestialBody;
                        model::AbstractDensityModel = Exponential())
    eop = fetch_iers_eop()
    return AtmosphericDrag(body, model, eop,
                           Type{<:AbstractVar}[PosVel], 6)
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
    jd  = t.utc.jd
    r   = SVector{3}(x̄[1], x̄[2], x̄[3])
    v   = SVector{3}(x̄[4], x̄[5], x̄[6])
    ρ   = density(force.model, jd, x̄, force.eop_data)
    BC  = geom.c_d * geom.drag_area / total_mass(sc)         # Cd·A/m [m²/kg]
    ω   = SVector{3}(0.0, 0.0, EARTH_ANGULAR_SPEED)
    v_app = v - cross(ω, r)                                   # transport theorem
    dfac  = -0.5 * BC * ρ * norm(v_app) * 1.0e3               # → km/s²
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = dfac * v_app[1]; x̄̇[5] = dfac * v_app[2]; x̄̇[6] = dfac * v_app[3]
    return x̄̇
end
