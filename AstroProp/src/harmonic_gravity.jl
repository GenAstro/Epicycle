# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
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
- `EGM96`, `EGM2008` for the Earth, `GL0660B` for the Moon, `JGM85F01` for Mars, and
  `IcgemGravity` for a field in any ICGEM file — full gravity fields. **Enterprise**.

# Writing your own
Define a type that subtypes `AbstractGeopotential` and give it `max_degree`, `max_order`,
`geopotential_data`, and `geopotential_accel`. `HarmonicGravity` then works with it unchanged. If
the coefficients are defined in axes other than the body's default orientation, also give it
[`field_orientation`](@ref).

!!! note "Enterprise"
    The full fields come from the `EpicycleEnterprise` package (commercial license). Load it to
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
    geopotential_accel(model, data, r_fixed, tsec, degree, order) -> SVector{3}

Total gravitational acceleration from the central body in the field's body-fixed axes [m/s²] — the
central term plus the harmonics up to `degree` and `order`. `r_fixed` is the position in those
axes [m]. Each gravity field provides this method, and `HarmonicGravity` calls it.
"""
function geopotential_accel end

"""
    field_orientation(model, body) -> AbstractOrientationModel

The body-fixed axes `model`'s coefficients are defined in, which `HarmonicGravity` evaluates it in.

# Notes
A gravity field is estimated in particular axes, and its coefficients describe the body's mass
only in those axes; evaluated in others, the field is rotated away from the mass that produced it.
The default is the body's orientation model, `orientation_model(body)`: the frame theory for the
Earth (ITRF), `LunarPA()` for the Moon, `IAU2015()` for the planets. A field defined in other axes
overrides this, as `JGM85F01` does with the IAU 1991 Mars axes.

`HarmonicGravity`'s `orientation` keyword takes precedence over both.
"""
field_orientation(::AbstractGeopotential, body) = orientation_model(body)

# ─────────────────────── body-fixed axes, shared with drag ───────────────────

# The rotation `R` and its rate `Ṙ` from a 6×6 body-fixed rotation `[R 0; Ṙ R]`, as static
# 3×3 matrices.
@inline function _rotation_blocks(M::AbstractMatrix)
    R = SMatrix{3,3}(M[1,1], M[2,1], M[3,1], M[1,2], M[2,2], M[3,2], M[1,3], M[2,3], M[3,3])
    Ṙ = SMatrix{3,3}(M[4,1], M[5,1], M[6,1], M[4,2], M[5,2], M[6,2], M[4,3], M[5,3], M[6,3])
    return R, Ṙ
end

# A force's axes must be the body's, checked when the force is built rather than at the first
# evaluation. AstroUniverse knows which bodies each shipped model orients; a model a user writes
# is taken at its word, as `set_orientation!` does.
function _check_axes(axes::AbstractOrientationModel, body::CelestialBody)
    AstroUniverse._orients(axes, body.naifid) || throw(ArgumentError(
        "$(axes) does not give the axes of $(body.name) (NAIF $(body.naifid)); pass " *
        "`orientation` a model of that body, or leave it out for the body's own."))
    return nothing
end

# ────────────────────────── spherical-harmonic gravity ───────────────────────
"""
    HarmonicGravity(body; degree, order, model = Zonal(), orientation = nothing)

Gravity from a body's non-spherical field, evaluated to the degree and order you choose.

# Arguments
- `body::CelestialBody`: the central body.
- `degree::Int`, `order::Int`: how far to evaluate the field; checked against what the chosen model
  supports.
- `model::AbstractGeopotential`: which gravity field. Open: [`Zonal`](@ref), the Earth's J2–J5 zonal
  field. Enterprise: the full fields of the Earth, the Moon, Mars, and any ICGEM file.
- `orientation::AbstractOrientationModel`: the body-fixed axes to evaluate the field in. Leave it
  out to use the axes the field is defined in, [`field_orientation`](@ref)`(model, body)`.

# Notes
The field is evaluated in body-fixed axes and the acceleration rotated back to the propagation
axes, ICRF. The axes are chosen when you construct the force: the `orientation` keyword if given,
otherwise the field's own, otherwise the body's orientation model. For the Earth that is the frame
theory in force at construction, so a GMAT comparison sets `set_frame_theory!(FK5())` first.
Changing a body's orientation or the frame theory afterwards does not change a force already built.

Pass `orientation` only deliberately. A field's coefficients describe the body's mass in the axes
they were estimated in, and in other axes the field is rotated away from it.

The gravity field is read once, at construction. Don't also add `PointMassGravity` for the same body
— that counts the central gravity twice, and `ForceModel` will stop you.

!!! note "Enterprise"
    The full fields come from the `EpicycleEnterprise` package (commercial license). The
    open-source version includes `Zonal` (J2–J5); switching to a full field changes only `model` —
    the rest of the call is the same.

# Examples
```julia
using AstroProp, AstroModels, AstroUniverse
grav = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())     # open, J2–J5

using EpicycleEnterprise
grav = HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())   # Enterprise, full field
grav = HarmonicGravity(mars; degree = 50, order = 50, model = JGM85F01()) # in its IAU 1991 axes
```
"""
struct HarmonicGravity{MT<:AbstractGeopotential, GD, OT<:AbstractOrientationModel} <: OrbitODE
    central_body::CelestialBody
    degree::Int
    order::Int
    model::MT
    data::GD
    orientation::OT
    dependencies::Vector{Type{<:AbstractVarTag}}
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
                         model::AbstractGeopotential = Zonal(),
                         orientation::Union{Nothing,AbstractOrientationModel} = nothing)
    _validate_degree_order(model, degree, order)
    data = geopotential_data(model, body, degree, order)
    axes = orientation === nothing ? field_orientation(model, body) : orientation
    _check_axes(axes, body)
    return HarmonicGravity(body, degree, order, model, data, axes,
                           Type{<:AbstractVarTag}[PosVel], 6)
end

function accel_eval!(force::HarmonicGravity, t::Time, x̄::Vector, x̄̇::Vector,
                     sc::Spacecraft, params; jac::Dict = Dict())
    # ICRF to the field's body-fixed axes; for the Earth, the frame theory's ITRF chain.
    R, _ = _rotation_blocks(body_fixed_rotation(force.orientation, force.central_body.naifid, t))
    r_fixed = R * SVector{3}(x̄[1], x̄[2], x̄[3]) .* 1.0e3          # km → m
    tsec    = (t.utc.jd - _JD_J2000) * 86400.0
    a_fixed = geopotential_accel(force.model, force.data, r_fixed, tsec,
                                 force.degree, force.order) ./ 1.0e3   # m/s² → km/s²
    a_icrf  = R' * a_fixed
    x̄̇[1] = x̄[4]; x̄̇[2] = x̄[5]; x̄̇[3] = x̄[6]
    x̄̇[4] = a_icrf[1]; x̄̇[5] = a_icrf[2]; x̄̇[6] = a_icrf[3]
    return x̄̇
end
