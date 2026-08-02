# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Spherical-harmonic gravity — API type plus the extension seam pluggable geopotential fields
# implement. Concrete fields (e.g. Zonal, EGM96, EGM2008) live in their own files.
#
# The call structure is adapted from Hammerhead-Space AstroForceModels.jl (MIT licensed) —
# see NOTICE / attribution.
#   Copyright (c) Hammerhead Space, MIT License.

const _JD_J2000 = 2451545.0

# ─────────────────────────── geopotential model seam ─────────────────────────
"""
    AbstractGeopotential

The gravity field used by [`HarmonicGravity`](@ref), chosen with its `model` keyword.

# Available models
- [`Zonal`](@ref) — Earth's zonal harmonics J2–J5. Open.
- `EGM96`, `EGM2008` — full gravity fields. **Enterprise**.

# Writing your own
Define a type that subtypes `AbstractGeopotential` and give it `max_degree`, `max_order`,
`geopotential_data`, and `geopotential_accel`. `HarmonicGravity` then works with it unchanged.

!!! note "Enterprise"
    `EGM96` and `EGM2008` come from the `EpicycleEnterprise` package (commercial license). Load it to
    use them; without it, `EGM96()` is undefined. See the Force Models guide for details.
"""
abstract type AbstractGeopotential end

"Maximum spherical-harmonic degree a geopotential model supports."
function max_degree end
"Maximum spherical-harmonic order a geopotential model supports."
function max_order end
"Model-specific data cached on the force at construction (coefficients, loaded field, …)."
function geopotential_data end

"""
    geopotential_accel(model, data, r_itrf, tsec, degree, order) -> SVector{3}

Total gravitational acceleration from the central body in the Earth-fixed frame [m/s²] — the central
term plus the harmonics up to `degree` and `order`. Each gravity field provides this method, and
`HarmonicGravity` calls it.
"""
function geopotential_accel end

# ────────────────────────── spherical-harmonic gravity ───────────────────────
"""
    HarmonicGravity(body; degree, order, model = Zonal())

Gravity from a body's non-spherical field, evaluated to the degree and order you choose.

# Arguments
- `body::CelestialBody`: the central body (Earth is the tested case).
- `degree::Int`, `order::Int`: how far to evaluate the field; checked against what the chosen model
  supports.
- `model::AbstractGeopotential`: which gravity field. Open: [`Zonal`](@ref), the J2–J5 zonal field.
  Enterprise: `EGM96`, `EGM2008`, the full fields.

# Notes
The gravity field and Earth-orientation data are read once, when you construct the force, so set them
up first. Don't also add `PointMassGravity` for the same body — that counts the central gravity
twice, and `ForceModel` will stop you.

!!! note "Enterprise"
    `EGM96` and `EGM2008` come from the `EpicycleEnterprise` package (commercial license). The
    open-source version includes `Zonal` (J2–J5); switching to a full field changes only `model` —
    the rest of the call is the same.

# Examples
```julia
using AstroProp, AstroModels, AstroUniverse
grav = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())     # open, J2–J5

using EpicycleEnterprise
grav = HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())   # Enterprise, full field
```
"""
struct HarmonicGravity{MT<:AbstractGeopotential, GD, EoT} <: OrbitODE
    central_body::CelestialBody
    degree::Int
    order::Int
    model::MT
    data::GD
    eop_data::EoT
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int
end

function _validate_degree_order(model::AbstractGeopotential, degree::Int, order::Int)
    degree < 0 && throw(ArgumentError("degree must be ≥ 0, got $degree"))
    order  < 0 && throw(ArgumentError("order must be ≥ 0, got $order"))
    order > degree &&
        throw(ArgumentError("order ($order) cannot exceed degree ($degree)"))
    md, mo = max_degree(model), max_order(model)
    degree > md && throw(ArgumentError(
        "$(nameof(typeof(model))) supports degree ≤ $md, got $degree"))
    order > mo && throw(ArgumentError(
        "$(nameof(typeof(model))) supports order ≤ $mo, got $order"))
    return nothing
end

function HarmonicGravity(body::CelestialBody; degree::Int, order::Int,
                         model::AbstractGeopotential = Zonal())
    _validate_degree_order(model, degree, order)
    data = geopotential_data(model, body, degree, order)
    eop  = fetch_iers_eop()
    return HarmonicGravity(body, degree, order, model, data, eop,
                           Type{<:AbstractVar}[PosVel], 6)
end

function accel_eval!(force::HarmonicGravity, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    jd = t.utc.jd
    R  = r_eci_to_ecef(J2000(), ITRF(), jd, force.eop_data)
    r_itrf = R * SVector{3}(x̄[1], x̄[2], x̄[3]) .* 1.0e3          # km → m
    tsec   = (jd - _JD_J2000) * 86400.0
    a_itrf = geopotential_accel(force.model, force.data, r_itrf, tsec,
                                force.degree, force.order) ./ 1.0e3   # m/s² → km/s²
    a_eci  = R' * a_itrf
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = a_eci[1]; x̄̇[5] = a_eci[2]; x̄̇[6] = a_eci[3]
    return x̄̇
end
