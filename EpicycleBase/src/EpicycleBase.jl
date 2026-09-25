# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

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
export get_field, set_field!
export state_jac!, param_jac!

# Traits on quantity functions — see quantity_traits.jl.
export label, set_quantity!, is_settable, is_cyclic, cycle, tag
export output_partial, has_output_partial

"""
    AbstractVar

`AbstractVar` is the common type for states, controls, time variables, and parameters.

# Example
```jldoctest
isabstracttype(AbstractVar)

# output
true
```
"""
abstract type AbstractVar end

"""
    AbstractState <: AbstractVar

`AbstractState` identifies state variables.

# Example
```jldoctest
AbstractState <: AbstractVar

# output
true
```
"""
abstract type AbstractState <: AbstractVar end

"""
    AbstractControl <: AbstractVar

`AbstractControl` identifies control variables.

# Example
```jldoctest
AbstractControl <: AbstractVar

# output
true
```
"""
abstract type AbstractControl <: AbstractVar end

"""
    AbstractTime <: AbstractVar

`AbstractTime` identifies time variables.

# Example
```jldoctest
AbstractTime <: AbstractVar

# output
true
```
"""
abstract type AbstractTime <: AbstractVar end

"""
    AbstractParam <: AbstractVar

`AbstractParam` identifies model parameters.

# Example
```jldoctest
AbstractParam <: AbstractVar

# output
true
```
"""
abstract type AbstractParam <: AbstractVar end

"""
    AbstractFun

`AbstractFun` is the common type for dynamics, outputs, and other model functions.

# Example
```jldoctest
isabstracttype(AbstractFun)

# output
true
```
"""
abstract type AbstractFun end

"""
    AlgebraicFun <: AbstractFun

`AlgebraicFun` identifies model functions that do not define differential equations.

# Example
```jldoctest
AlgebraicFun <: AbstractFun

# output
true
```
"""
abstract type AlgebraicFun <: AbstractFun end

"""
    AbstractCalcVariable

`AbstractCalcVariable` is the common type for calculated orbit, body, and maneuver quantities.

# Example
```jldoctest
isabstracttype(AbstractCalcVariable)

# output
true
```
"""
abstract type AbstractCalcVariable end

"""
    AbstractOrbitVar <: AbstractCalcVariable

`AbstractOrbitVar` identifies calculated orbit quantities.

# Example
```jldoctest
AbstractOrbitVar <: AbstractCalcVariable

# output
true
```
"""
abstract type AbstractOrbitVar    <: AbstractCalcVariable end

"""
    AbstractBodyVar <: AbstractCalcVariable 

`AbstractBodyVar` identifies calculated celestial-body quantities.

# Example
```jldoctest
AbstractBodyVar <: AbstractCalcVariable

# output
true
```
"""
abstract type AbstractBodyVar     <: AbstractCalcVariable end

"""
    AbstractManeuverVar <: AbstractCalcVariable 

`AbstractManeuverVar` identifies calculated maneuver quantities.

# Example
```jldoctest
AbstractManeuverVar <: AbstractCalcVariable

# output
true
```
"""
abstract type AbstractManeuverVar <: AbstractCalcVariable end

"""
    AbstractOrbitStateType <: AbstractOrbitVar

`AbstractOrbitStateType` identifies orbital state representations such as Cartesian and Keplerian states.

# Example
```jldoctest
AbstractOrbitStateType <: AbstractOrbitVar

# output
true
```
"""
abstract type AbstractOrbitStateType <: AbstractOrbitVar end  

"""
    AbstractPoint

`AbstractPoint` is the common type for geometric points such as spacecraft and celestial bodies.

# Example
```jldoctest
isabstracttype(AbstractPoint)

# output
true
```
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

`AbstractVarTag` is the common type for singleton tags that identify scalar fields on model
objects. Tags are zero-size structs used for method selection and Jacobian block indexing.

# Hierarchy
    AbstractVarTag
    ├── AbstractStateTag   (state components, e.g. PosVel)
    ├── AbstractParamTag   (parameters, e.g. Mu, Mass)
    ├── AbstractControlTag (control variables; none defined yet)
    └── AbstractTimeTag    (time variables; none defined yet)

# Example
```jldoctest
isabstracttype(AbstractVarTag)

# output
true
```
"""
abstract type AbstractVarTag end

"""
    AbstractStateTag <: AbstractVarTag

`AbstractStateTag` identifies tags for components of the ODE state vector, such as `PosVel`.

# Example
```jldoctest
AbstractStateTag <: AbstractVarTag

# output
true
```
"""
abstract type AbstractStateTag   <: AbstractVarTag end

"""
    AbstractParamTag <: AbstractVarTag

`AbstractParamTag` identifies tags for scalar parameters such as `Mu` and `Mass`.

# Example
```jldoctest
AbstractParamTag <: AbstractVarTag

# output
true
```
"""
abstract type AbstractParamTag   <: AbstractVarTag end

"""
    AbstractControlTag <: AbstractVarTag

`AbstractControlTag` identifies tags for control variables and is not yet used by an Epicycle model.

# Example
```jldoctest
AbstractControlTag <: AbstractVarTag

# output
true
```
"""
abstract type AbstractControlTag <: AbstractVarTag end

"""
    AbstractTimeTag <: AbstractVarTag

`AbstractTimeTag` identifies tags for time variables and is not yet used by an Epicycle model.

# Example
```jldoctest
AbstractTimeTag <: AbstractVarTag

# output
true
```
"""
abstract type AbstractTimeTag    <: AbstractVarTag end

"""
    ModelVariable{T<:AbstractVarTag}

`ModelVariable` pairs a model object with a tag to identify one scalar variable.

# Fields
- `model`: The owning object (e.g., `CelestialBody`, `Spacecraft`)
- `tag::T`: Tag singleton identifying which scalar field on `model`

# Notes
Two `ModelVariable`s are equal, and are the same key in a `Dict`, when their tags have the same
type and their models are identical. For a mutable model such as a `Spacecraft` or a
`CelestialBody`, identical means the same instance, so two instances of the same type are
distinct variables. For an immutable model such as a `NamedTuple`, identical means equal field
values, so two separately built models with the same values are the same variable.

# Example
```julia
struct DryMassTag <: AbstractParamTag end
model = (dry_mass = 850.0,)
variable = ModelVariable(model, DryMassTag())
```
"""
struct ModelVariable{T<:AbstractVarTag}
    model
    tag::T
end

"""
    DirectVariable(; value::Real, lower_bound::Real, upper_bound::Real, name::AbstractString)

`DirectVariable` represents a scalar decision variable that is not a field on a model object,
such as a ΔV component, time of flight, or scale factor.

# Fields
- `value::Float64`
- `lower_bound::Float64`
- `upper_bound::Float64`
- `name::String`

# Example
```julia
tof = DirectVariable(value = 3.0, lower_bound = 1.0,
                     upper_bound = 10.0, name = "Time of flight")
```
"""
struct DirectVariable
    value::Float64
    lower_bound::Float64
    upper_bound::Float64
    name::String

    function DirectVariable(; value::Real, lower_bound::Real, upper_bound::Real, name::AbstractString)
        return new(Float64(value), Float64(lower_bound), Float64(upper_bound), String(name))
    end
end

"""
    get_field(model, tag::AbstractVarTag) -> value

Return the field on `model` identified by `tag`.
Concrete methods are defined in each package alongside their model struct.

# Returns
The field's current value: a scalar for a parameter tag such as `Mu`, and a vector for a state
tag such as `PosVel`, whose value is the 6-element position and velocity.

# Example
```jldoctest
struct RadiusTag <: AbstractParamTag end
EpicycleBase.get_field(model::NamedTuple, ::RadiusTag) = model.radius
get_field((radius = 6378.1363,), RadiusTag())

# output
6378.1363
```
"""
function get_field end

"""
    get_field(var::ModelVariable) -> value

Convenience: read the field value through a `ModelVariable`.
"""
get_field(v::ModelVariable) = get_field(v.model, v.tag)

"""
    set_field!(model, tag::AbstractVarTag, value)

Set the field on `model` identified by `tag` to `value`: a scalar for a parameter tag, a vector
for a state tag such as `PosVel`.
Concrete methods are defined in each package alongside their model struct.

# Returns
Nothing a caller relies on; `model` is mutated.

# Example
```jldoctest
struct MassTag <: AbstractParamTag end
mutable struct VehicleMass
    mass::Float64
end
EpicycleBase.set_field!(model::VehicleMass, ::MassTag, value::Real) =
    (model.mass = value; nothing)
vehicle = VehicleMass(1000.0)
set_field!(vehicle, MassTag(), 950.0)
vehicle.mass

# output
950.0
```
"""
function set_field! end

"""
    set_field!(var::ModelVariable, value)

Convenience: write the field value through a `ModelVariable`.
"""
set_field!(v::ModelVariable, value) = set_field!(v.model, v.tag, value)

"""
    state_jac!(out, force, t, y, sc)

Analytic contribution of `force` to A = ∂f/∂y.  Accumulates into `out` (6×6).

# Arguments
- `out`: The 6×6 Jacobian being assembled. Add this force's contribution to it.
- `force`: The force, a subtype of AstroProp's `OrbitODE`.
- `t`: The epoch, an AstroEpochs `Time`.
- `y`: The state, position and velocity in km and km/s.
- `sc`: The spacecraft the force acts on.

# Notes
Defining a method for a concrete force type is the registration; no separate
call is needed. A force with no method is differentiated automatically by
AstroProp's `eval_jacobian!`.

Only rows 4 to 6, the acceleration, are used. The caller writes rows 1 to 3,
∂ṙ/∂v = I, once for the whole force model, so anything a method adds there is
discarded.

# Returns
Nothing a caller relies on; the contribution is added into `out`.

# Example
```jldoctest
struct LinearForce end
EpicycleBase.state_jac!(out, ::LinearForce, t, y, sc) =
    (out[4, 1] += -2.0; nothing)
A = zeros(6, 6)
state_jac!(A, LinearForce(), 0.0, zeros(6), nothing)
A[4, 1]

# output
-2.0
```
"""
function state_jac! end

"""
    param_jac!(out, force, tag, t, y, sc)

Analytic contribution of `force` to B = ∂f/∂p for the parameter identified by
`tag`.  Accumulates into `out` (6-element vector).

# Arguments
- `out`: The 6-element column being assembled. Add this force's contribution to it.
- `force`: The force, a subtype of AstroProp's `OrbitODE`.
- `tag`: The parameter's tag, such as `Mu`.
- `t`: The epoch, an AstroEpochs `Time`.
- `y`: The state, position and velocity in km and km/s.
- `sc`: The spacecraft the force acts on.

# Notes
Defining a method for `(force_type, tag_type)` is the registration. A pair with
no method is differentiated by a central finite difference in AstroProp's
`eval_jacobian!`.

Only elements 4 to 6, the acceleration, are used. The caller sets elements 1 to
3 to zero, since the kinematics do not depend on a force parameter.

# Returns
Nothing a caller relies on; the contribution is added into `out`.

# Example
```jldoctest
struct GravityForce end
struct GravityParameterTag <: AbstractParamTag end
EpicycleBase.param_jac!(out, ::GravityForce, ::GravityParameterTag, t, y, sc) =
    (out[4, 1] += 0.25; nothing)
B = zeros(6, 1)
param_jac!(B, GravityForce(), GravityParameterTag(), 0.0, zeros(6), nothing)
B[4, 1]

# output
0.25
```
"""
function param_jac! end

include("quantity_traits.jl")

end
