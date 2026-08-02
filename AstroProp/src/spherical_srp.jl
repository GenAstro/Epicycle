# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Solar radiation pressure — cannonball (spherical) SRP with a pluggable shadow model.
# Shadow-factor formulas adapted from Hammerhead-Space AstroForceModels.jl (MIT) and
# Montenbruck & Gill §3.4.2.
#   Copyright (c) Hammerhead Space, MIT License.

const _C_M_S = 2.99792458e8        # speed of light [m/s]

# ─────────────────────────────── shadow model seam ───────────────────────────
"""
    AbstractShadowModel

The eclipse model used by [`SolarRadiationPressure`](@ref), chosen with its `shadow` keyword — it
gives the fraction of the Sun's disk that is unobscured by the occulting body.

# Available models
- [`DualCone`](@ref) — umbra + penumbra dual-cone eclipse. Open.

# Writing your own
Define a type that subtypes `AbstractShadowModel` and give it a
`_shadow_factor(model, r_sat, r_sun, R_sun, R_occ)` method returning a lighting factor in ``[0, 1]``.
`SolarRadiationPressure` then works with it unchanged.
"""
abstract type AbstractShadowModel end

"""
    DualCone()

A dual-cone shadow with both umbra and penumbra — a realistic eclipse (Montenbruck & Gill,
*Satellite Orbits*, §3.4.2).

# Examples
```julia
SolarRadiationPressure(earth; shadow = DualCone())
```
"""
struct DualCone <: AbstractShadowModel end


# ─────────────────────────── solar radiation pressure ────────────────────────
"""
    SolarRadiationPressure(body; shadow = DualCone(),
                             solar_flux = 1367.0, nominal_sun = 149597870.691)

Solar radiation pressure on the spacecraft — the push of sunlight, with an eclipse model for the
central body's shadow.

The acceleration follows the cannonball model:

```math
\\vec{a}_{\\mathrm{SRP}} = F \\cdot \\frac{C_r A}{m} \\cdot \\frac{\\Phi}{c}
                          \\left(\\frac{r_{\\mathrm{AU}}}{|\\vec{r}_{s/c \\to \\odot}|}\\right)^{\\!2}
                          \\hat{r}_{\\odot \\to s/c}
```

where ``F \\in [0, 1]`` is the lighting factor from `shadow`, ``\\Phi`` is `solar_flux`, ``c`` is the
speed of light, and ``r_{\\mathrm{AU}}`` is `nominal_sun`.

# Arguments
- `body::CelestialBody`: positional, required — the central body whose shadow can eclipse the
  spacecraft.
- `shadow::AbstractShadowModel`: the eclipse model — [`DualCone`](@ref) (umbra + penumbra).
- `solar_flux::Real`: the solar irradiance at 1 AU, W/m². Must be positive. Default `1367.0`
  (GMAT `SRP.Flux`).
- `nominal_sun::Real`: the reference Sun distance used for 1/r² scaling, km. Must be positive.
  Default `149597870.691` (GMAT `SRP.Nominal_Sun`).

# Notes
Reads reflectivity ``C_r`` and area ``A`` from the spacecraft's `SphericalSRP`, and mass from the
spacecraft — set `sc.srp` before propagating. The Sun's position is looked up from the ephemeris
via `AstroUniverse.translate(body, sun, jd_tdb)` in km, inertial.

Cross-validated against GMAT (`SRPModel = Spherical`); see `test/force_srp_spherical.jl`.

# Examples
```julia
using AstroProp, AstroUniverse
srp = SolarRadiationPressure(earth; shadow = DualCone(), solar_flux = 1367.0)
```
"""
struct SolarRadiationPressure{S<:AbstractShadowModel} <: OrbitODE
    central_body::CelestialBody
    shadow::S
    solar_flux::Float64        # W/m² at 1 AU   (GMAT SRP.Flux)
    nominal_sun::Float64       # km             (GMAT SRP.Nominal_Sun; the AU used for 1/r² scaling)
    R_sun::Float64             # km
    R_occ::Float64             # km (occulting body radius)
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int
end

function SolarRadiationPressure(body::CelestialBody;
                                  shadow::AbstractShadowModel = DualCone(),
                                  solar_flux::Real = 1367.0,           # GMAT default
                                  nominal_sun::Real = 149597870.691)   # GMAT default [km]
    solar_flux > 0 || throw(ArgumentError(
        "solar_flux must be positive; got solar_flux = $solar_flux W/m²"))
    nominal_sun > 0 || throw(ArgumentError(
        "nominal_sun must be positive; got nominal_sun = $nominal_sun km"))
    return SolarRadiationPressure(body, shadow, Float64(solar_flux), Float64(nominal_sun),
                                  sun.equatorial_radius, body.equatorial_radius,
                                  Type{<:AbstractVar}[PosVel], 6)
end

function Base.show(io::IO, ::MIME"text/plain", f::SolarRadiationPressure)
    println(io, "SolarRadiationPressure:")
    println(io, "  central_body = ", f.central_body.name)
    println(io, "  shadow       = ", nameof(typeof(f.shadow)), "()")
    println(io, "  solar_flux   = ", f.solar_flux,  " W/m²")
    println(io, "  nominal_sun  = ", f.nominal_sun, " km")
end
Base.show(io::IO, f::SolarRadiationPressure) = show(io, MIME"text/plain"(), f)

# Numerically stable angle between two vectors (Montenbruck/AstroForceModels form).
@inline function _angle_between(v1, v2)
    u1 = v1 ./ norm(v1)
    u2 = v2 ./ norm(v2)
    return 2.0 * atan(norm(u1 .- u2), norm(u1 .+ u2))
end

function _shadow_factor(::DualCone, r_sat, r_sun, R_sun, R_occ)
    R_ss = r_sat .- r_sun
    a = asin(R_sun / norm(R_ss))          # apparent radius of the Sun
    b = asin(R_occ / norm(r_sat))         # apparent radius of the occulting body
    c = _angle_between(R_ss, r_sat)       # apparent separation
    if c ≥ (b + a)
        return 1.0                        # full sun
    elseif c < (b - a)
        return 0.0                        # umbra
    elseif c < (a - b)
        return 1.0 - b^2 / a^2            # occulting body fully within the Sun disk
    else                                  # penumbra
        x = (c^2 + a^2 - b^2) / (2.0 * c)
        y = sqrt(a^2 - x^2)
        area = a^2 * acos(x / a) + b^2 * acos((c - x) / b) - c * y
        return 1.0 - area / (π * a^2)
    end
end

function accel_eval!(force::SolarRadiationPressure, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    geom = sc.srp
    geom === nothing && throw(ArgumentError(
        "sc.srp must be a SphericalSRP for SolarRadiationPressure; got nothing. " *
        "Set sc.srp = SphericalSRP(; c_r, srp_area) before propagating."))
    jd_tdb = t.tdb.jd
    r_sat  = SVector{3}(x̄[1], x̄[2], x̄[3])
    r_sun  = SVector{3}(translate(force.central_body, sun, jd_tdb))   # Earth→Sun, km, inertial
    RC = geom.c_r * geom.srp_area / total_mass(sc)                    # Cr·A/m [m²/kg]
    Ψ  = force.solar_flux / _C_M_S                                    # N/m²
    F  = _shadow_factor(force.shadow, r_sat, r_sun, force.R_sun, force.R_occ)
    R_ss  = r_sat - r_sun
    d     = norm(R_ss)
    F_srp = F * RC * Ψ * (force.nominal_sun / d)^2 / 1.0e3            # → km/s², away from Sun
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = F_srp * R_ss[1] / d
    x̄̇[5] = F_srp * R_ss[2] / d
    x̄̇[6] = F_srp * R_ss[3] / d
    return x̄̇
end
