# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Spherical-harmonic gravity and atmospheric drag force models.
#
# The numeric kernels (ECI↔ECEF rotation, geopotential evaluation, NRLMSISE-00 density)
# are computed directly with the SatelliteToolbox packages. The call structure is adapted
# from Hammerhead-Space AstroForceModels.jl (MIT licensed) — see NOTICE / attribution.
#   Copyright (c) Hammerhead Space, MIT License.

using SatelliteToolboxTransformations: r_eci_to_ecef, ecef_to_geodetic, fetch_iers_eop,
                                       J2000, ITRF, DCM
using SatelliteToolboxGravityModels: GravityModels, IcgemFile, fetch_icgem_file
import SatelliteToolboxBase: EARTH_ANGULAR_SPEED
using StaticArrays: SVector
using LinearAlgebra: cross, norm, dot
using AstroUniverse: earth, sun, CelestialBody
using AstroModels: total_mass

const _JD_J2000 = 2451545.0
const _C_M_S    = 2.99792458e8        # speed of light [m/s]
const _AU_KM    = 1.495978707e8       # astronomical unit [km]

# ─────────────────────────── model-selection tags ────────────────────────────
"Supertype for published geopotential coefficient sets selected by `HarmonicGravity`."
abstract type AbstractGeopotential end
struct EGM96   <: AbstractGeopotential end
struct EGM2008 <: AbstractGeopotential end
_icgem_symbol(::EGM96)   = :EGM96
_icgem_symbol(::EGM2008) = :EGM2008
_load_geopotential(m::AbstractGeopotential) =
    GravityModels.load(IcgemFile, fetch_icgem_file(_icgem_symbol(m)))

"""
    AbstractDensityModel

Interface for atmospheric density models consumed by `AtmosphericDrag`. A density model
implements `density(model, jd, r_eci, eop) -> ρ` [kg/m³]. `Exponential` is the open
reference model; higher-fidelity models (e.g. NRLMSISE-00) are added by the enterprise
package, which subtypes this and adds a `density` method — no change to the open package.
"""
abstract type AbstractDensityModel end
struct Exponential <: AbstractDensityModel end

"Supertype for space-weather providers on the drag force."
abstract type AbstractSpaceWeather end

"""
    ConstantSpaceWeather(; f107 = 150.0, f107a = 150.0, magnetic_index = 3.0)

Constant space-weather indices (GMAT-style `Drag.F107` / `F107A` / `MagneticIndex`).
Held on `AtmosphericDrag`. (NRLMSISE-00 here reads the SpaceIndices tables; a constant
provider that overrides those tables is future work — see the force-model spec.)
"""
Base.@kwdef struct ConstantSpaceWeather <: AbstractSpaceWeather
    f107::Float64           = 150.0
    f107a::Float64          = 150.0
    magnetic_index::Float64 = 3.0
end

# ────────────────────────── spherical-harmonic gravity ───────────────────────
"""
    HarmonicGravity(body; degree, order, model = EGM96())

Central-body spherical-harmonic gravity to a user-specified `degree` and `order`, using
a published coefficient set (`model`). Loads the coefficients and Earth-orientation data
at construction.
"""
struct HarmonicGravity{GT, EoT, MT<:AbstractGeopotential} <: OrbitODE
    central_body::CelestialBody
    degree::Int
    order::Int
    model::MT
    gravity_model::GT
    eop_data::EoT
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int
end

function HarmonicGravity(body::CelestialBody; degree::Int, order::Int,
                         model::AbstractGeopotential = EGM96())
    coeffs = _load_geopotential(model)
    eop    = fetch_iers_eop()
    return HarmonicGravity(body, degree, order, model, coeffs, eop,
                           Type{<:AbstractVar}[PosVel], 6)
end

function accel_eval!(force::HarmonicGravity, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    jd = t.utc.jd
    R  = r_eci_to_ecef(J2000(), ITRF(), jd, force.eop_data)
    r_itrf = R * SVector{3}(x̄[1], x̄[2], x̄[3]) .* 1.0e3          # km → m
    tsec   = (jd - _JD_J2000) * 86400.0
    a_itrf = GravityModels.gravitational_acceleration(
                 force.gravity_model, r_itrf, tsec;
                 max_degree = force.degree, max_order = force.order) ./ 1.0e3   # m/s² → km/s²
    a_eci = R' * a_itrf
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = a_eci[1]; x̄̇[5] = a_eci[2]; x̄̇[6] = a_eci[3]
    return x̄̇
end

# ───────────────────────────── atmospheric drag ──────────────────────────────
"""
    AtmosphericDrag(; body = earth, model = Exponential(), space_weather = nothing)

Atmospheric drag force. Reads the spacecraft's `CannonballDrag` geometry and
`total_mass(sc)` at evaluation time; `model` selects the density model (any
`AbstractDensityModel`). The open package provides `Exponential`; enterprise density models
(e.g. `MSISE00`) plug in through the same field.
"""
struct AtmosphericDrag{EoT, DM<:AbstractDensityModel, SW} <: OrbitODE
    central_body::CelestialBody
    model::DM
    eop_data::EoT
    space_weather::SW
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int
end

function AtmosphericDrag(; body::CelestialBody = earth,
                           model::AbstractDensityModel = Exponential(),
                           space_weather = nothing)
    eop = fetch_iers_eop()
    return AtmosphericDrag(body, model, eop, space_weather,
                           Type{<:AbstractVar}[PosVel], 6)
end

function _geodetic(jd, x̄, eop)
    R    = r_eci_to_ecef(DCM, J2000(), ITRF(), jd, eop)
    ecef = R * SVector{3}(x̄[1], x̄[2], x̄[3])
    return ecef_to_geodetic(ecef .* 1.0e3)          # (lat, lon, alt) [rad, rad, m]
end

"""
    density(model::AbstractDensityModel, jd, r_eci, eop) -> ρ   [kg/m³]

Atmospheric mass density at the spacecraft: `jd` the UTC Julian date, `r_eci` the inertial
position [km], `eop` the Earth-orientation data. This is the extension point enterprise
density models add methods to; the open package implements `Exponential`.
"""
function density end

density(::Exponential, jd, x̄, eop) = _exponential_density(_geodetic(jd, x̄, eop)[3])

# Exponential atmosphere (Vallado, Fundamentals of Astrodynamics 4th ed., Table 8-4).
# Native reimplementation adapted from SatelliteToolboxAtmosphericModels.jl (MIT).
# _EXP_H0 base altitude [km], _EXP_RHO0 nominal density [kg/m³], _EXP_H scale height [km].
const _EXP_H0 = (0.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0,
                 130.0, 140.0, 150.0, 180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
                 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0)
const _EXP_RHO0 = (1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3, 3.206e-4, 8.770e-5,
                   1.905e-5, 3.396e-6, 5.297e-7, 9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9,
                   2.070e-9, 5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11, 9.518e-12,
                   3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.170e-14,
                   5.245e-15, 3.019e-15)
const _EXP_H  = (7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549, 5.799, 5.382, 5.877,
                 7.263, 9.473, 12.636, 16.149, 22.523, 29.740, 37.105, 45.546, 53.628,
                 53.298, 58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268.00)

function _exponential_density(h_m::Real)
    h_m < 0 && throw(ArgumentError("altitude must be ≥ 0, got $h_m m"))
    h  = h_m / 1000                          # → km
    id = 28
    @inbounds for i in 1:28
        if _EXP_H0[i] - h > 0
            id = i - 1
            break
        end
    end
    id = max(id, 1)
    @inbounds return _EXP_RHO0[id] * exp(-(h - _EXP_H0[id]) / _EXP_H[id])
end

function accel_eval!(force::AtmosphericDrag, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    geom = sc.drag
    geom === nothing && error("AtmosphericDrag requires spacecraft drag geometry: " *
                              "set sc.drag = CannonballDrag(; c_d, drag_area).")
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

# ─────────────────────────── solar radiation pressure ────────────────────────
"Supertype for eclipse/shadow models used by `SolarRadiationPressure`."
abstract type AbstractShadowModel end
struct Cylindrical <: AbstractShadowModel end
"Dual-cone (umbra + penumbra) conical shadow — the \"bi-conic\" model (Montenbruck & Gill §3.4.2)."
struct DualCone    <: AbstractShadowModel end
struct NoShadow    <: AbstractShadowModel end

"""
    SolarRadiationPressure(; body = earth, shadow = DualCone(), solar_flux = 1361.0)

Cannonball solar-radiation-pressure force. Reads the spacecraft's `CannonballSRP` geometry
and `total_mass(sc)`; gets the Sun position from the ephemeris via `translate(body, sun, …)`.
`shadow` selects the eclipse model; `solar_flux` is the irradiance at 1 AU [W/m²].
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

function SolarRadiationPressure(; body::CelestialBody = earth,
                                  shadow::AbstractShadowModel = DualCone(),
                                  solar_flux::Real = 1367.0,           # GMAT default
                                  nominal_sun::Real = 149597870.691)   # GMAT default [km]
    return SolarRadiationPressure(body, shadow, Float64(solar_flux), Float64(nominal_sun),
                                  sun.equatorial_radius, body.equatorial_radius,
                                  Type{<:AbstractVar}[PosVel], 6)
end

# Numerically stable angle between two vectors (Montenbruck/AstroForceModels form).
@inline function _angle_between(v1, v2)
    u1 = v1 ./ norm(v1)
    u2 = v2 ./ norm(v2)
    return 2.0 * atan(norm(u1 .- u2), norm(u1 .+ u2))
end

# Lighting factor F ∈ [0,1]: 1 = full sun, 0 = umbra.
_shadow_factor(::NoShadow, r_sat, r_sun, R_sun, R_occ) = 1.0

function _shadow_factor(::Cylindrical, r_sat, r_sun, R_sun, R_occ)
    ŝ  = r_sun ./ norm(r_sun)
    dp = dot(ŝ, r_sat)
    return (dp >= 0.0 || norm(r_sat .- dp .* ŝ) > R_occ) ? 1.0 : 0.0
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
    geom === nothing && error("SolarRadiationPressure requires spacecraft SRP geometry: " *
                              "set sc.srp = CannonballSRP(; c_r, srp_area).")
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
