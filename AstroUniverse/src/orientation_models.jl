# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Body orientation models.
#
# How a body is oriented in space is a property of the body, so it is set here
# and not on a coordinate system — the same rule frame theory and EOP follow.
# `AstroFrames` asks the universe how a body is oriented; it does not decide.
#
# For most planets and major moons the published models are good enough and
# ship as defaults. For a body whose orientation is poorly known — an asteroid,
# a comet, a newly visited moon — you supply your own model, and for a body
# whose orientation is being estimated you also need to read and write the
# parameters being solved for. Both go through the same interface:
#
#     set_orientation!(bennu, MySpinModel(...))
#     orientation_parameters(bennu)                       # what can be solved for
#     set_orientation_parameters!(bennu, (; pole_ra = x))  # what the estimator writes
#
# Writing a model means one method, `body_axes_rotation`. If the model has the
# usual pole-and-prime-meridian shape, `pole_axes_rotation` builds the matrix
# from six angles and you supply only those.
#
# Everything here is plain Julia and differentiable, so a model whose
# parameters carry derivative information works without special handling. The
# SPICE-backed model is the exception and says so.
# =============================================================================

"""
    AbstractOrientationModel

Supertype for models describing how a body is oriented in space.

# Interface
An abstract type, so it has no fields of its own. A model must provide:

- [`body_axes_rotation`](@ref)`(model, naifid, jd_tdb)` — **required.** The
  rotation from inertial to body-fixed axes at that epoch.
- [`orientation_parameters`](@ref)`(model)` — optional. Defaults to the model's
  fields, which is right for a model written as a plain struct.
- [`set_orientation_parameters`](@ref)`(model, nt)` — optional. Defaults to
  copying the model with those fields replaced.

The two optional ones are what an estimation reads and writes, and the defaults
mean a plain struct needs no extra code for them.

Shipped models: [`IauPolynomialOrientation`](@ref) for the Sun, the planets and
Pluto, and [`SpiceOrientation`](@ref) for a body whose orientation comes from a
kernel.

# Notes
Set a model on a body with [`set_orientation!`](@ref); read it back with
[`orientation_model`](@ref). Until a body has one, body-fixed axes for it raise
and say so.

A custom model covers a body Epicycle ships none for, such as an asteroid,
a comet or a newly visited moon. It also covers a body whose orientation is being
estimated, since a kernel cannot carry derivative information.

# Example
```julia
struct SimpleSpin{T<:Real} <: AbstractOrientationModel
    pole_ra::T          # deg
    pole_dec::T         # deg
    pm0::T              # deg at J2000
    spin_rate::T        # deg/day
end

function AstroUniverse.body_axes_rotation(m::SimpleSpin, naifid, jd_tdb)
    d = jd_tdb - 2451545.0
    Ẇ = deg2rad(m.spin_rate)
    return pole_axes_rotation(deg2rad(m.pole_ra), deg2rad(m.pole_dec),
                              deg2rad(m.pm0) + Ẇ * d,
                              zero(Ẇ), zero(Ẇ), Ẇ / 86_400)
end

set_orientation!(bennu, SimpleSpin(85.46, -60.36, 89.6, 2011.145))
```
"""
abstract type AbstractOrientationModel end

# --- Elementary rotations, local to this file -------------------------------
#
# Kept private rather than shared with AstroFrames' equivalents: this package
# sits below it and cannot depend on it, and these are four lines of
# trigonometry rather than a fact with a home.

@inline function _rot_z(θ::Real)
    s, c = sincos(θ)
    return @SMatrix [ c    s   zero(θ)
                     -s    c   zero(θ)
                     zero(θ) zero(θ) one(θ)]
end

@inline function _rot_x(θ::Real)
    s, c = sincos(θ)
    return @SMatrix [one(θ) zero(θ) zero(θ)
                     zero(θ)  c    s
                     zero(θ) -s    c]
end

@inline function _drot_z(θ::Real)
    s, c = sincos(θ)
    return @SMatrix [-s    c   zero(θ)
                     -c   -s   zero(θ)
                     zero(θ) zero(θ) zero(θ)]
end

@inline function _drot_x(θ::Real)
    s, c = sincos(θ)
    return @SMatrix [zero(θ) zero(θ) zero(θ)
                     zero(θ) -s    c
                     zero(θ) -c   -s]
end

"""
    pole_axes_rotation(α₀, δ₀, W, α̇₀, δ̇₀, Ẇ) -> SMatrix{6,6}

Build the rotation from inertial axes to body-fixed axes, given where the
body's north pole points and how far its prime meridian has turned.

# Arguments
- `α₀`, `δ₀` — right ascension and declination of the body's north pole in
  inertial axes, in radians.
- `W` — angle from the node of the body's equator on the inertial equator to
  the body's prime meridian, in radians.
- `α̇₀`, `δ̇₀`, `Ẇ` — their rates, in **radians per second**. A spin rate quoted
  in degrees per day needs `deg2rad(rate) / 86_400`; getting this wrong is the
  common mistake and shows up as a rotation period out by 86 400.

The result is a 6×6 that rotates position and velocity together, so the
`ω × r` term every rotating frame contributes is already in it. It suits a model
of the usual shape, a pole direction and a spin about it, which covers nearly
every planet, moon and asteroid:

```julia
function AstroUniverse.body_axes_rotation(m::MySpin, naifid, jd_tdb)
    d = jd_tdb - 2451545.0
    W = m.W₀ + m.spin_rate * d
    return pole_axes_rotation(m.α₀, m.δ₀, W, 0.0, 0.0, m.spin_rate / 86400)
end
```

The construction is `R_z(W) · R_x(π/2 − δ₀) · R_z(π/2 + α₀)`, the convention
of the IAU working group on cartographic coordinates and rotational elements
(Archinal et al., 2018 §2.1).

# Returns
A 6×6 matrix that rotates a state, position stacked on velocity, from inertial
axes into the body's rotating axes. Position in the top three rows,
velocity in the bottom three, with the `ω × r` term every rotating frame
contributes already included.

# Notes
The rotation is `R_z(W) · R_x(π/2 − δ₀) · R_z(π/2 + α₀)`, the convention of the
IAU working group on cartographic coordinates and rotational elements
(Archinal et al., *Celestial Mechanics and Dynamical Astronomy* 130:22, 2018,
§2.1). Row 3 of the rotation is the body's north pole in inertial axes.

Differentiable in all six arguments, so an estimated pole carries derivative
information through the rotation.

Angles are not reduced to a range; `W` may be many turns' worth and the
trigonometry handles it.

# Example
```jldoctest
julia> using AstroUniverse, LinearAlgebra

julia> M = pole_axes_rotation(deg2rad(317.68), deg2rad(52.89), 1.0, 0.0, 0.0, 7.088e-5);

julia> size(M)
(6, 6)

julia> round(det(M[1:3, 1:3]), digits = 12)   # a proper rotation
1.0
```
"""
function pole_axes_rotation(α₀::Real, δ₀::Real, W::Real,
                            α̇₀::Real, δ̇₀::Real, Ẇ::Real)
    α₀, δ₀, W, α̇₀, δ̇₀, Ẇ = promote(α₀, δ₀, W, α̇₀, δ̇₀, Ẇ)

    θ1, θ2, θ3 = π/2 + α₀, π/2 - δ₀, W

    Rz1, dRz1 = _rot_z(θ1), _drot_z(θ1)
    Rx2, dRx2 = _rot_x(θ2), _drot_x(θ2)
    Rz3, dRz3 = _rot_z(θ3), _drot_z(θ3)

    R = Rz3 * Rx2 * Rz1
    # Chain rule over the three factors: dθ1/dt = α̇₀, dθ2/dt = −δ̇₀, dθ3/dt = Ẇ.
    Ṙ = α̇₀ * Rz3 * Rx2 * dRz1 +
        (-δ̇₀) * Rz3 * dRx2 * Rz1 +
        Ẇ * dRz3 * Rx2 * Rz1

    Z = zero(eltype(R))
    return @SMatrix [R[1,1] R[1,2] R[1,3]  Z      Z      Z
                     R[2,1] R[2,2] R[2,3]  Z      Z      Z
                     R[3,1] R[3,2] R[3,3]  Z      Z      Z
                     Ṙ[1,1] Ṙ[1,2] Ṙ[1,3] R[1,1] R[1,2] R[1,3]
                     Ṙ[2,1] Ṙ[2,2] Ṙ[2,3] R[2,1] R[2,2] R[2,3]
                     Ṙ[3,1] Ṙ[3,2] Ṙ[3,3] R[3,1] R[3,2] R[3,3]]
end

"""
    body_axes_rotation(model, naifid, jd_tdb) -> SMatrix{6,6}

Rotation from inertial axes to the body-fixed axes of body `naifid` at TDB
Julian date `jd_tdb`.

This is the one method an orientation model has to implement.

# Arguments
- `model` — the orientation model, an [`AbstractOrientationModel`](@ref).
- `naifid::Integer` — NAIF ID of the body being oriented. Passed so one model
  type can serve several bodies. A model for one body may leave it unused.
- `jd_tdb::Real` — epoch as a Julian date in barycentric dynamical time (TDB).

# Returns
A 6×6 matrix rotating a state from inertial axes into the body's rotating
axes, as [`pole_axes_rotation`](@ref) returns.

# Notes
[`pole_axes_rotation`](@ref) builds the result from the pole direction and the
prime-meridian angle with their rates, and does the rest. Only a model of an
unusual shape needs to build the matrix itself.

Write it in plain Julia if the body's orientation is being estimated, or if a
trajectory's derivatives run through these axes. A model that calls out to
SPICE cannot carry derivative information; see [`SpiceOrientation`](@ref).

# Example
```julia
struct SimpleSpin <: AbstractOrientationModel
    pole_ra::Float64      # deg
    pole_dec::Float64     # deg
    pm0::Float64          # deg at J2000
    spin_rate::Float64    # deg/day
end

function AstroUniverse.body_axes_rotation(m::SimpleSpin, naifid, jd_tdb)
    d = jd_tdb - 2451545.0
    Ẇ = deg2rad(m.spin_rate)
    return pole_axes_rotation(deg2rad(m.pole_ra), deg2rad(m.pole_dec),
                              deg2rad(m.pm0) + Ẇ * d, 0.0, 0.0, Ẇ / 86_400)
end
```
"""
function body_axes_rotation end

# --- Reading and writing the parameters being solved for --------------------

"""
    orientation_parameters(model_or_body) -> NamedTuple

The orientation parameters that can be estimated, by name.

Defaults to the model's fields, so a model written as a plain struct needs no
extra code. A model with nothing to estimate returns an empty `NamedTuple`.

# Arguments
- `model_or_body` — an [`AbstractOrientationModel`](@ref), or a body, in which
  case its currently set model is read.

# Returns
A `NamedTuple` of parameter names and their current values. Units follow the
model's documented convention.

# Notes
Defaults to the model's fields, so a model written as a plain struct needs no
extra code for this to work. A model with nothing to estimate, such as the shipped
polynomials or a SPICE frame, returns an empty `NamedTuple`, which is the honest
answer rather than an omission.

Whatever manages an estimation reads the names here and writes back through
[`set_orientation_parameters!`](@ref). This package does not track which of
them are being solved for.

# Example
```jldoctest
julia> using AstroUniverse

julia> struct Spin <: AbstractOrientationModel
           pole_ra::Float64
           pole_dec::Float64
       end

julia> orientation_parameters(Spin(85.46, -60.36))
(pole_ra = 85.46, pole_dec = -60.36)

julia> orientation_parameters(IauPolynomialOrientation())
NamedTuple()
```
"""
function orientation_parameters(m::AbstractOrientationModel)
    names = fieldnames(typeof(m))
    return NamedTuple{names}(ntuple(i -> getfield(m, i), length(names)))
end

"""
    set_orientation_parameters(model, nt::NamedTuple) -> model

A copy of `model` with the named parameters replaced. Names not given keep
their current values.

# Arguments
- `model` — the model to copy.
- `nt::NamedTuple` — parameter names and their new values. Every name must be
  one the model has; an unknown name raises rather than being ignored.

# Returns
A new model of the same kind, with those parameters replaced.

# Notes
Returns a new model rather than mutating one, so a value carrying derivative
information can change the model's numeric type without invalidating anything
holding the old model.

A model whose fields are all one numeric type, which is the usual way to write
one, cannot be rebuilt from a mixture, and an estimated value arrives on its own
while the others stay as they were. The numeric fields are widened to their
common type when that happens, so writing a single solved-for value works on
the model most people write.

To write onto a body rather than a model, use
[`set_orientation_parameters!`](@ref).

# Example
```jldoctest
julia> using AstroUniverse

julia> struct Spin2 <: AbstractOrientationModel
           pole_ra::Float64
           pole_dec::Float64
       end

julia> orientation_parameters(set_orientation_parameters(Spin2(85.46, -60.36), (; pole_ra = 85.5)))
(pole_ra = 85.5, pole_dec = -60.36)
```
"""
function set_orientation_parameters(m::M, nt::NamedTuple) where {M<:AbstractOrientationModel}
    names = fieldnames(M)
    for k in keys(nt)
        k in names || throw(ArgumentError(
            "`$(k)` is not an orientation parameter of $(nameof(M)). " *
            "It has $(isempty(names) ? "none" : join(string.(names), ", ")). " *
            "Check the name, or use a model that carries the parameter you mean."))
    end
    isempty(names) && return m

    vals = ntuple(i -> get(nt, names[i], getfield(m, i)), length(names))
    W = Base.typename(M).wrapper
    applicable(W, vals...) && return W(vals...)

    # A model written the obvious way — `struct M{T<:Real}` with every field
    # typed `T` — cannot be built from a mixture, and an estimated parameter
    # arrives on its own while the rest stay as they were. Widen the numeric
    # fields to a common type and try again. Without this, writing a single
    # solved-for value into a model fails on the model most people write.
    return W(_promote_numeric(vals)...)
end

"""
    _promote_numeric(vals::Tuple) -> Tuple

`vals` with its numbers widened to their common type and everything else left
alone. A non-numeric field, such as a frame name or a flag, is not a parameter
and must not be touched.
"""
function _promote_numeric(vals::Tuple)
    nums = filter(v -> v isa Real, vals)
    isempty(nums) && return vals
    T = reduce(promote_type, map(typeof, nums))
    return map(v -> v isa Real ? convert(T, v) : v, vals)
end

# --- Shipped model: the published polynomials -------------------------------

"""
    IauPolynomialOrientation()

The published orientation for the Sun, the planets and Pluto: pole direction
and prime meridian as polynomials in time, with the periodic terms where the
IAU report gives them (Archinal et al., 2018).

This is the default for those bodies, and it applies until
[`set_orientation!`](@ref) is called. Its coefficients are fixed, so it has no
parameters to estimate; solving for a pole needs a model that carries the pole as
a field.

# Fields
None. The coefficients are the published ones and are not settable.

# Notes
Differentiable in time, so a trajectory's derivatives pass through body-fixed
axes for these bodies.

Covers NAIF IDs 10, 199, 299, 499, 599, 699, 799, 899 and 999. Earth and the
Moon are deliberately absent: neither is described by a pole and a prime
meridian, and both have their own frames.

# Example
```jldoctest
julia> using AstroUniverse

julia> orientation_model(jupiter)
IauPolynomialOrientation()

julia> size(body_axes_rotation(IauPolynomialOrientation(), 599, 2458849.5))
(6, 6)
```
"""
struct IauPolynomialOrientation <: AbstractOrientationModel end

function body_axes_rotation(::IauPolynomialOrientation, naifid::Integer, jd_tdb::Real)
    o = iau2015_orientation(naifid, jd_tdb)
    return pole_axes_rotation(o.ra_pole, o.dec_pole, o.prime_meridian,
                              o.ra_pole_rate, o.dec_pole_rate, o.prime_meridian_rate)
end

# --- Shipped model: a SPICE frame -------------------------------------------

"""
    SpiceOrientation(frame_name)

Orientation read from a loaded SPICE frame, such as `"IAU_MARS"` or
`"MOON_PA"`.

Use this for a body whose orientation is published as a kernel. The kernel
supplies periodic and libration terms without duplicating its coefficients.

# Fields
- `frame::String` — the SPICE frame name, such as `"IAU_MARS"` or `"MOON_PA"`.
  The kernels defining it must already be loaded; a text PCK defines the
  `IAU_<BODY>` frames.

# Example
```julia
load_spice_kernel("pck00011.tpc")
set_orientation!(mars, SpiceOrientation("IAU_MARS"))

cs = CoordinateSystem(mars, CelestialBodyFixed())
```

# Notes
**Not differentiable.** SPICE computes in double precision behind a C
interface, so derivative information cannot pass through it. A body whose
orientation is being estimated, or a trajectory whose gradient runs through
body-fixed axes, needs a model written in Julia; see
[`pole_axes_rotation`](@ref). This model has no parameters to estimate, because
nothing it returns can be solved for.
"""
struct SpiceOrientation <: AbstractOrientationModel
    frame::String
end

# The frame name is not a solve-for.
orientation_parameters(::SpiceOrientation) = NamedTuple()
set_orientation_parameters(m::SpiceOrientation, nt::NamedTuple) =
    isempty(nt) ? m : throw(ArgumentError(
        "A SPICE frame has no parameters to estimate, because its orientation comes " *
        "from the kernel. Register a model written in Julia to estimate orientation; " *
        "see `pole_axes_rotation`."))

function body_axes_rotation(m::SpiceOrientation, naifid::Integer, jd_tdb::Real)
    # SPICE ephemeris time is seconds past the J2000 epoch.
    et = (jd_tdb - _J2000_TDB_JD) * _SECONDS_PER_DAY
    return SMatrix{6,6,Float64,36}(sxform("J2000", m.frame, et))
end

# --- Which model a body uses ------------------------------------------------
#
# Held here rather than as a field on `CelestialBody` so that registering a
# model does not change the layout of a type other packages construct, and so
# that two references to the same body cannot disagree. Keyed on NAIF ID,
# which is what identifies a body.

const _ORIENTATION_MODELS = Dict{Int,AbstractOrientationModel}()
const _ORIENTATION_LOCK   = ReentrantLock()

const _IAU_POLYNOMIAL_BODIES = (10, 199, 299, 499, 599, 699, 799, 899, 999)

@inline _naifid(x::Integer) = Int(x)
@inline _naifid(b::CelestialBody) = b.naifid

"""
    set_orientation!(body, model)

Set how `body` is oriented in space.

```julia
set_orientation!(mars,  SpiceOrientation("IAU_MARS"))
set_orientation!(bennu, MySpinModel(85.46, -60.36, 89.6, 2011.145))
```

# Arguments
- `body` — a body, or a NAIF ID.
- `model` — an [`AbstractOrientationModel`](@ref).

# Returns
The model, so a registration can be written inline.

# Notes
Replaces whatever the body used before, including the shipped default. The
setting is keyed on the body's NAIF ID, so every reference to that body sees it.
Orientation is a property of the body rather than of one variable holding it.

Setting an orientation is what makes body-fixed axes available for a body
Epicycle does not ship a model for. Until one is set, asking for those axes
raises and says so.

# Example
```julia
bennu = CelestialBody("Bennu", 4.892e-9, 0.2825, 0.0, 2101955)
set_orientation!(bennu, MySpinModel(85.46, -60.36, 89.6, 2011.145))

cs = CoordinateSystem(bennu, CelestialBodyFixed())   # now works
```
"""
function set_orientation!(body, model::AbstractOrientationModel)
    n = _naifid(body)
    lock(_ORIENTATION_LOCK) do
        _ORIENTATION_MODELS[n] = model
    end
    return model
end

"""
    orientation_model(body) -> AbstractOrientationModel

# Arguments
- `body` — a body, or a NAIF ID.

# Returns
The body's [`AbstractOrientationModel`](@ref) — what was set with
[`set_orientation!`](@ref), or the shipped default for a body that has one.

# Notes
The Sun, the planets, and Pluto default to [`IauPolynomialOrientation`](@ref).
Earth and the Moon have no default and raise pointing at their own frames,
which carry Earth orientation parameters and lunar libration respectively —
neither is described by a pole and a prime meridian.

Any other body raises until a model is registered. Use
[`has_orientation_model`](@ref) to ask without raising.

# Example
```jldoctest
julia> using AstroUniverse

julia> orientation_model(mars)
IauPolynomialOrientation()

julia> has_orientation_model(earth)
false
```
"""
function orientation_model(body)
    n = _naifid(body)
    m = lock(_ORIENTATION_LOCK) do
        get(_ORIENTATION_MODELS, n, nothing)
    end
    m === nothing || return m

    n in _IAU_POLYNOMIAL_BODIES && return IauPolynomialOrientation()

    if n == 399
        throw(ArgumentError(
            "Earth's rotation is not described by a pole-and-meridian model. " *
            "Use the Earth frames instead, ITRF and the chain reaching it, which read " *
            "Earth orientation parameters."))
    elseif n == 301
        throw(ArgumentError(
            "The Moon has its own frames, which carry libration: MoonME for surface " *
            "and mapping work, MoonPA for gravity and dynamics. They differ by about " *
            "875 m on the surface, so pick deliberately. To use body-fixed axes " *
            "anyway, call `set_orientation!(moon, SpiceOrientation(\"MOON_PA\"))`."))
    end
    throw(ArgumentError(
        "No orientation model for NAIF ID $(n). Bodies outside the Sun, the planets, " *
        "and Pluto need one supplied: `set_orientation!(body, model)`. " *
        "Use `SpiceOrientation(\"IAU_<NAME>\")` if a kernel defines the frame, or write " *
        "your own model — see `pole_axes_rotation`."))
end

"""
    has_orientation_model(body) -> Bool

Whether `body` can be used with body-fixed axes.

# Arguments
- `body` — a body, or a NAIF ID.

# Returns
`true` if a model is set or a shipped default applies, `false` otherwise.

# Notes
The question [`orientation_model`](@ref) answers by raising. `false` for Earth
and the Moon, which have their own frames rather than a pole-and-meridian
model.

# Example
```jldoctest
julia> using AstroUniverse

julia> has_orientation_model(mars), has_orientation_model(earth)
(true, false)
```
"""
function has_orientation_model(body)
    n = _naifid(body)
    lock(_ORIENTATION_LOCK) do
        haskey(_ORIENTATION_MODELS, n)
    end && return true
    return n in _IAU_POLYNOMIAL_BODIES
end

"""
    orientation_parameters(body) -> NamedTuple

The orientation parameters of `body` that can be estimated, by name.
"""
orientation_parameters(body::CelestialBody) = orientation_parameters(orientation_model(body))

"""
    set_orientation_parameters!(body, nt::NamedTuple)

Write orientation parameters of `body` — the estimated values.

```julia
set_orientation_parameters!(bennu, (; pole_ra = 85.51, pole_dec = -60.30))
```

# Arguments
- `body` — a body, or a NAIF ID.
- `nt::NamedTuple` — parameter names and their new values. Names not given keep
  their current values.

# Returns
The updated model.

# Notes
**Mutates the body's registered orientation** — every reference to that body
sees the new values afterwards. This is how whatever manages an estimation
writes a solution back; this package does not track which parameters are being
solved for.

Raises if the body has no orientation model, or if a name is not one of its
parameters.

# Example
```julia
orientation_parameters(bennu)                          # what can be written
set_orientation_parameters!(bennu, (; pole_ra = 85.51))
```
"""
function set_orientation_parameters!(body, nt::NamedTuple)
    updated = set_orientation_parameters(orientation_model(body), nt)
    return set_orientation!(body, updated)
end
