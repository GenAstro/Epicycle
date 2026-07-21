# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

"""
    module EpicycleBase

Core abstract types shared across Epicycle (variables, states, controls, time, functions, points).
These form the public type hierarchy used by higher-level packages (AstroStates, AstroEpochs, AstroFrames,
AstroProp, AstroSolve, …).
"""
module EpicycleBase

export AbstractVar, AbstractState, AbstractControl, AbstractTime, AbstractParam
export AbstractFun, AlgebraicFun
export AbstractCalcVariable, AbstractOrbitVar, AbstractBodyVar, AbstractManeuverVar
export AbstractOrbitStateType
export AbstractPoint

export AbstractVarTag, AbstractStateTag, AbstractParamTag, AbstractControlTag, AbstractTimeTag
export ModelVariable, DirectVariable
export get_field, set_field!, differentiate_wrt, fd_differentiate_wrt
export state_jac!, param_jac!

"""
    AbstractVar

Base tag for all variable kinds (states, controls, time, parameters).

# Notes:
- Serves as the common supertype for variable categories.
"""
abstract type AbstractVar end

"""
    AbstractState <: AbstractVar

Base type for all state variables.

"""
abstract type AbstractState <: AbstractVar end

"""
    AbstractControl <: AbstractVar

Base type for all control variables.
"""
abstract type AbstractControl <: AbstractVar end

"""
    AbstractTime <: AbstractVar

Base type for time variables.
"""
abstract type AbstractTime <: AbstractVar end

"""
    AbstractParam <: AbstractVar

Base type for parameter variables.
"""
abstract type AbstractParam <: AbstractVar end

"""
    AbstractFun

Base type for function objects (e.g., dynamics, outputs).
"""
abstract type AbstractFun end

"""
    AlgebraicFun <: AbstractFun

Base type for algebraic (non-differential) function objects.
"""
abstract type AlgebraicFun <: AbstractFun end

"""
    AbstractCalcVariable

Base type for calculation variables (orbit, body, maneuver).
"""
abstract type AbstractCalcVariable end

"""
    AbstractOrbitVar <: AbstractCalcVariable

Base type for orbit calculation variables.
"""
abstract type AbstractOrbitVar    <: AbstractCalcVariable end

"""
    AbstractBodyVar <: AbstractCalcVariable 

Base type for body calculation variables.
"""
abstract type AbstractBodyVar     <: AbstractCalcVariable end

"""
    AbstractManeuverVar <: AbstractCalcVariable 

Base type for maneuver calculation variables.
"""
abstract type AbstractManeuverVar <: AbstractCalcVariable end

"""
    AbstractOrbitStateType <: AbstractOrbitVar

Base type for orbit state representation types (e.g., `Cartesian()`, `Keplerian()`, etc.).
"""
abstract type AbstractOrbitStateType <: AbstractOrbitVar end  

"""
    AbstractPoint

Base type for geometric points (e.g., Spacecraft, CelestialBody).
"""
abstract type AbstractPoint end

"""
    no_op()

A no-op function that returns nothing. Useful as a default callback placeholder.

# Notes:
- Unexported utility; reference as EpicycleBase.no_op.
"""
function no_op()
    return nothing
end

# =============================================================================
# Tag / Variable System
# =============================================================================

"""
    AbstractVarTag

Root supertype for all singleton tag types that identify scalar fields on model
objects.  Tags are zero-size structs used as keys for dispatch and Jacobian
block indexing.

# Hierarchy
    AbstractVarTag
    ├── AbstractStateTag   (state components, e.g. PosVel)
    ├── AbstractParamTag   (parameters, e.g. Mu, Cd, Cr, Mass)
    ├── AbstractControlTag (control variables, e.g. DvX, DvY, DvZ)
    └── AbstractTimeTag    (time variables — stub only)
"""
abstract type AbstractVarTag end

"""
    AbstractStateTag <: AbstractVarTag

Supertype for tags identifying components of the ODE state vector (e.g., `PosVel`).
"""
abstract type AbstractStateTag   <: AbstractVarTag end

"""
    AbstractParamTag <: AbstractVarTag

Supertype for tags identifying scalar parameters (e.g., `Mu`, `Cd`, `Cr`, `Mass`).
"""
abstract type AbstractParamTag   <: AbstractVarTag end

"""
    AbstractControlTag <: AbstractVarTag

Supertype for tags identifying scalar control variables (e.g., `DvX`, `DvY`, `DvZ`).
"""
abstract type AbstractControlTag <: AbstractVarTag end

"""
    AbstractTimeTag <: AbstractVarTag

Supertype for tags identifying time-related variables.  Stub — not yet used.
"""
abstract type AbstractTimeTag    <: AbstractVarTag end

"""
    ModelVariable{T<:AbstractVarTag}

Pairs a model object with a tag singleton to form a unique variable identity.

# Fields
- `model`: The owning object (e.g., `CelestialBody`, `Spacecraft`)
- `tag::T`: Tag singleton identifying which scalar field on `model`

# Identity key
    (objectid(model), typeof(tag))

Two `ModelVariable`s with the same model instance and tag type are equal.
Two different model instances (even of the same concrete type) produce distinct keys.
"""
struct ModelVariable{T<:AbstractVarTag}
    model
    tag::T
end

"""
    DirectVariable

A scalar decision variable that is not a field on any model object.  Used for
NLP variables such as ΔV components, time-of-flight, or scale factors.

# Fields
- `value::Float64`
- `lower_bounds::Float64`
- `upper_bounds::Float64`
- `name::String`
"""
struct DirectVariable
    value::Float64
    lower_bounds::Float64
    upper_bounds::Float64
    name::String

    function DirectVariable(; value::Real, lower_bounds::Real, upper_bounds::Real, name::AbstractString)
        return new(Float64(value), Float64(lower_bounds), Float64(upper_bounds), String(name))
    end
end

"""
    get_field(model, tag::AbstractVarTag) -> Float64

Return the scalar field on `model` identified by `tag`.
Concrete methods are defined in each package alongside their model struct.
"""
function get_field end

"""
    get_field(var::ModelVariable) -> Float64

Convenience: read the field value through a `ModelVariable`.
"""
get_field(v::ModelVariable) = get_field(v.model, v.tag)

"""
    set_field!(model, tag::AbstractVarTag, value::Real)

Set the scalar field on `model` identified by `tag` to `value`.
Concrete methods are defined in each package alongside their model struct.
"""
function set_field! end

"""
    set_field!(var::ModelVariable, value::Real)

Convenience: write the field value through a `ModelVariable`.
"""
set_field!(v::ModelVariable, value::Real) = set_field!(v.model, v.tag, value)

"""
    differentiate_wrt(rhs_fn, model, tag::AbstractVarTag) -> Vector{Float64}

Partial derivative of `rhs_fn()` with respect to the scalar field on `model`
identified by `tag`.

The **default** implementation uses 2-point central finite differences.  This
fallback makes any user-defined custom model immediately usable without analytic
work.  Packages that ship built-in models define a more specific analytic method
alongside the model struct — Julia dispatch selects the analytic override
automatically.

# Arguments
- `rhs_fn`: Zero-argument closure returning the current RHS as a `Vector{Float64}`
- `model`:  The owning model object
- `tag`:    Tag singleton identifying the scalar parameter

# Notes
- The model field is temporarily mutated and then restored.  The model must be
  mutable (i.e., `set_field!` must be defined) for this fallback to work.
- Step size: `h = cbrt(eps(Float64)) * max(1.0, |p₀|)` — optimal for 2nd-order
  central FD on well-scaled parameters.
"""
function differentiate_wrt(rhs_fn, model, tag::AbstractVarTag)
    p0 = get_field(model, tag)
    h  = cbrt(eps(Float64)) * max(1.0, abs(p0))
    set_field!(model, tag, p0 + h);  f_plus  = rhs_fn()
    set_field!(model, tag, p0 - h);  f_minus = rhs_fn()
    set_field!(model, tag, p0)
    return (f_plus .- f_minus) ./ (2h)
end

"""
    fd_differentiate_wrt(rhs_fn, model, tag::AbstractVarTag) -> Vector{Float64}

Finite-difference oracle for `differentiate_wrt`.  Always uses 2-point central
finite differences, even when an analytic method exists for `(model, tag)`.
Intended for use in tests to validate analytic implementations.

See also: [`differentiate_wrt`](@ref)
"""
function fd_differentiate_wrt(rhs_fn, model, tag::AbstractVarTag)
    p0 = get_field(model, tag)
    h  = cbrt(eps(Float64)) * max(1.0, abs(p0))
    set_field!(model, tag, p0 + h);  f_plus  = rhs_fn()
    set_field!(model, tag, p0 - h);  f_minus = rhs_fn()
    set_field!(model, tag, p0)
    return (f_plus .- f_minus) ./ (2h)
end

"""
    state_jac!(out, force, t, y, sc)

Analytic contribution of `force` to A = ∂f/∂y.  Accumulates into `out` (6×6).

Defining a method for a concrete force type IS the registration — no separate
call is needed.  `eval_jacobian!` detects the method via `hasmethod` at
`JacobianResult` construction and calls it on every ODE step.

If no method is defined for a force type, `eval_jacobian!` falls back to
`ForwardDiff.jacobian` on `accel_eval!` for that force in isolation.
"""
function state_jac! end

"""
    param_jac!(out, force, tag, t, y, sc)

Analytic contribution of `force` to B = ∂f/∂p for the parameter identified by
`tag`.  Accumulates into `out` (6×1).

Defining a method for `(force_type, tag_type)` IS the registration.
`eval_jacobian!` detects it via `hasmethod` and calls it on every ODE step.

If no method is defined for a `(force, tag)` pair, `eval_jacobian!` falls back
to a central-FD perturbation of `accel_eval!` for that force/tag pair only.
"""
function param_jac! end

end