# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

__precompile__()

"""
    module AstroCallbacks

Defines interfaces to define and evaluate solve-for variables and constraint functions.
"""
module AstroCallbacks

# TODO: REFACTOR THIS CODE TO USE GET/SET METHODS ON MODELS AND CALL THEM. 
#       CURRENTLY THIS CODE ACCESSES FIELDS DIRECTLY, WHICH IS BAD PRACTICE.

using LinearAlgebra

using AstroStates
using AstroEpochs
using EpicycleBase
using AstroUniverse
using AstroFrames
using AstroModels: Spacecraft, to_posvel, set_posvel!
using AstroManeuvers
using StaticArrays

import AstroStates: state_tag_to_type, state_type_to_tag

export AbstractCalcVariable, AbstractOrbitVar, AbstractBodyVar, AbstractManeuverVar
export AbstractCalc, OrbitCalc, BodyCalc, ManeuverCalc
export get_calc, set_calc!, calc_numvars, calc_is_settable, calc_input_statetag
export PositionVector, VelocityVector, PosMag, SMA, TA, RAAN
export position_vector, velocity_vector, raan
export state, state!

# Quantity readers, and the deferred form. See quantity_readers.jl and calc.jl.
export element, epoch, epoch!
export delta_v, delta_v!
export semi_major_axis, eccentricity, inclination, argument_of_periapsis, true_anomaly
export semi_major_axis!, eccentricity!, inclination!, raan!, argument_of_periapsis!,
       true_anomaly!
export position_x, position_y, position_z, position_magnitude, velocity_magnitude
export position_x!, position_y!, position_z!, position_vector!, velocity_vector!
export position_magnitude!, velocity_magnitude!
export position_dot_velocity
export mean_long_sma, mean_long_sma!, outgoing_rla, delta_v_magnitude, delta_v_magnitude!
export gravitational_parameter, gravitational_parameter!
export Calc, reapply
export history, quantities
export label, set_quantity!, is_settable, output_partial   # from EpicycleBase, so writing a quantity needs only AstroCallbacks
export OutGoingRLA, PosX, PosZ, VelMag, Ecc, Inc, PosDotVel, MeanSMA

# export Maneuver Variables
export DeltaVMag, DeltaVVector

# export Body Variables 
export GravParam

import EpicycleBase: AbstractCalcVariable, AbstractOrbitVar, AbstractBodyVar, AbstractManeuverVar

"""
    AbstractCalc

The common type of the older calculation structs, `OrbitCalc`, `BodyCalc` and `ManeuverCalc`,
each of which pairs a subject with a calculation tag. `get_calc` reads one and `set_calc!` writes
it. New scripts use quantity functions and [`Calc`](@ref) instead.

# Example
```julia
calc = OrbitCalc(sat, SMA())
calc isa AbstractCalc
```
"""
abstract type AbstractCalc end

# Trait to default number of variables for a calc variable to 1
calc_numvars(::AbstractCalcVariable) = 1     # COV_EXCL_LINE unreachable but safe check.

"""
    OrbitCalc(sc::Spacecraft, var::AbstractOrbitVar; dependency = nothing)

Calc struct for set/get orbit-derived variables for a spacecraft.

Fields
- sc::Spacecraft: Subject spacecraft.
- var::AbstractOrbitVar: Orbit variable tag to evaluate (e.g., PositionVector()).
- dep: Optional dependency (e.g., another Spacecraft) or `nothing`.

# Notes:
- For a list of all supported BodyCalc variables, see `?AbstractOrbitVar`.
- Conversions that require μ obtain it from the `CoordinateSystem` origin when available.
- Use `get_calc(::OrbitCalc)` to evaluate 
- Use `set_calc!(::OrbitCalc, values)` to assign when supported.

# Examples
```julia
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
    time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
)
oc = OrbitCalc(sc, PositionVector())
r  = get_calc(oc)
set_calc!(oc, [7000.0, 300.0, 0.0])
```
"""
struct OrbitCalc{S<:Spacecraft,V<:AbstractOrbitVar,D} <: AbstractCalc
    sc::S
    var::V
    dep::D
end

""" 
    OrbitCalc(sc::Spacecraft, var::V; dependency=nothing) where {V<:AbstractOrbitVar}

Outer constructor for OrbitCalc
"""
OrbitCalc(sc::Spacecraft, var::V; dependency=nothing) where {V<:AbstractOrbitVar} =
    OrbitCalc{typeof(sc),V,typeof(dependency)}(sc, var, dependency)

"""
    BodyCalc(body::CelestialBody, var::AbstractBodyVar)

Calc struct for set/get body-derived variables for a CelestialBody

Fields
- body::B where B is a `CelestialBody`
- var::V where V<:AbstractBodyVar

# Notes:
- For a list of all supported BodyCalc variables, see `?AbstractBodyVar`.

# Examples
```julia
mu_calc = BodyCalc(earth, GravParam())
muval = get_calc(mu_calc)
set_calc!(mu_calc, 3.986e5)
```
"""
struct BodyCalc{B,V<:AbstractBodyVar} <: AbstractCalc
    body::B
    var::V
end

"""
    ManeuverCalc(man, sc::Spacecraft, var::AbstractManeuverVar)

Calc struct for set/get maneuver-derived variables for a Maneuver and Spacecraft.

Fields
- man::M where M is a maneuver model (e.g., ImpulsiveManeuver)
- sc::S where S<:Spacecraft
- var::V where V<:AbstractManeuverVar

# Notes:
- See ?AbstractManeuverVar for supported maneuver variables (e.g., DeltaVVector).
- Use get_calc(::ManeuverCalc) to evaluate and set_calc!(::ManeuverCalc, values)
  to assign when supported.

# Examples
```julia
m = ImpulsiveManeuver(axes=Inertial(), Isp=300.0, element1=0.01, element2=0.02,
    element3=-0.03)
sc = Spacecraft(
    state     = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
    time      = Time("2020-01-01T00:00:00", TAI(), ISOT()),
)

mc_vec = ManeuverCalc(m, sc, DeltaVVector())
dv = get_calc(mc_vec)
set_calc!(mc_vec, [0.01, 0.02, -0.03])

mc_mag = ManeuverCalc(m, sc, DeltaVMag())
dvm = get_calc(mc_mag)
```
"""
struct ManeuverCalc{M,S,V<:AbstractManeuverVar} <: AbstractCalc
    man::M
    sc::S
    var::V
end

# Import calcs for each variable
for f in (
    "orbitcalc_orbitstates.jl",
    "orbitcalc_posmag.jl",
    "orbitcalc_positionvector.jl",
    "orbitcalc_velocityvector.jl",
    "orbitcalc_sma.jl",
    "orbitcalc_meansma.jl",
    "orbitcalc_inc.jl",
    "orbitcalc_ta.jl",
    "orbitcalc_raan.jl",
    "orbitcalc_outgoingrla.jl",
    "orbitcalc_posx.jl",
    "orbitcalc_posz.jl",
    "orbitcalc_velmag.jl",
    "orbitcalc_ecc.jl",
    "orbitcalc_pos_dot_vel.jl",
    "bodycalc_gravparam.jl",
    "maneuvercalc_deltavvector.jl",
    "maneuvercalc_deltavmag.jl",
    "frame_aware_quantities.jl",
    "quantity_readers.jl",
    "calc.jl",
    "history.jl",
)
    include(f)
end

""" 
    function _subjects_from_calc(c::AbstractCalc)

Extract the references used in a Calc instance. For example, an OrbitCalc
will return the spacecraft it references.
"""
_subjects_from_calc(c::AstroCallbacks.OrbitCalc)    = (c.sc,)
_subjects_from_calc(c::AstroCallbacks.ManeuverCalc) = (c.man, c.sc)
_subjects_from_calc(c::AstroCallbacks.BodyCalc)     = (c.body,)

# A `Calc` states its subject in the same slot every spec does, so there is
# nothing to dig for. Note what is *not* here: `ManeuverCalc` above also
# reports `c.sc`, because the old form carried a spacecraft it never used to
# evaluate a delta-V. A spacecraft that matters is named by some other spec in
# the sequence — a constraint on it, or the propagate that moves it.
_subjects_from_calc(c::Calc) = (first(c.args),)

"""
    function get_calc(c::AbstractCalc)

Compute the value of a Calc variable. Generic get_calc methods 
interface simple Calc types. More complicated types (i.e. OrbitCalc
have custom interfaces)  
"""
@inline get_calc(c::BodyCalc)     = _evaluate(c.var, c.body)
@inline get_calc(c::ManeuverCalc) = _evaluate(c.var, c.man, c.sc)

""" 
    calc_input_statetag(v::AbstractOrbitVar)

State representation required by an older orbit-variable tag.

# Returns
An orbit-state type tag. Throws an error when the variable has not declared one.

# Example
```jldoctest
calc_input_statetag(SMA())

# output

Keplerian()
```
"""
calc_input_statetag(v::AbstractOrbitVar) = begin
    tname = string(nameof(typeof(v)))
    error("There is no AbstractOrbitVar defined named $(tname)(). " *
          "See `subtypes(AbstractOrbitVar)` for defined quantities.")
end

"""
    function _extract_mu(cs::CoordinateSystem)

Extract μ from a coordinate system origin if it exists
"""
function _extract_mu(cs::CoordinateSystem)
    origin = cs.origin
    T = typeof(origin)
    if hasfield(T, :mu)
        return getfield(origin, :mu)
    else
        # If Nothing, and not OK in context, needs to be trapped in caller
        return nothing
    end
end

"""
    function to_concrete_state(os::OrbitState)::AbstractState 

Converts an OrbitState to its equivalent concrete AbstractState

# Arguments
- `os::OrbitState`: The OrbitState instance to convert.
"""
function to_concrete_state(os::OrbitState)::AbstractState 
    # TODO. Probably should be in AstroStates
    TS = state_tag_to_type(os.statetype) 
    v  = os.state                     
    return TS(copy(v))
end

"""
    function to_concrete_state(os::AbstractState)::AbstractState 

Pass-through for already-concrete state type conversion

# Arguments
- `os::AbstractState`: The AbstractState instance to convert.
"""
# TODO. Probably should be in AstroStates
to_concrete_state(s::AbstractState) = s

"""
    function convert_orbitcalc_state(st::AbstractOrbitState, cs::CoordinateSystem, 
                            target::AbstractOrbitStateType)::AbstractOrbitState

Convert AbstractOrbitState to state types and coordinates needed for OrbitCalc variable.
"""
function convert_orbitcalc_state(
    st::AbstractOrbitState,
    cs::CoordinateSystem,
    target::AbstractOrbitStateType,
)::AbstractOrbitState
    tostate = state_tag_to_type(target)

    # Fast path: already in target concrete type
    if st isa tostate
        return st  # COV_EXCL_LINE
    end

    # Try conversion without μ first; if method is missing, retry with μ
    try
        return tostate(st)
    catch e
        if e isa MethodError
            μ = _extract_mu(cs)
            if μ === nothing
                typein = typeof(st)
                msg = "AstroCallbacks.convert_orbitcalc_state: μ is required but missing; " *
                      "cannot convert orbit state from $(typein) to $(target) using " *
                      "CoordinateSystem(origin=$(typeof(cs.origin)), " *
                      "axes=$(typeof(cs.axes))). Use a celestial-body origin (with μ) " *
                      "or request a state type that does not require μ."
                error(msg)
            end
            return tostate(st, μ)
        else
            rethrow(e)
        end
    end
end

"""
    function convert_orbitcalc_state(os::OrbitState, cs::CoordinateSystem, target::AbstractOrbitStateType)::AbstractOrbitState

Convert OrbitState to state types and coordinates needed for OrbitCalc variable.
"""
function convert_orbitcalc_state(
    os::OrbitState,
    cs::CoordinateSystem,
    target::AbstractOrbitStateType,
)::AbstractOrbitState
    return convert_orbitcalc_state(to_concrete_state(os), cs, target)
end

"""
    function _state_for_calc(c::OrbitCalc)::AbstractOrbitState

Compute the state type required for the orbit calc variable.  Performs conversion if needed.
"""
function _state_for_calc(c::OrbitCalc)::AbstractOrbitState
    sc   = c.sc
    dep  = c.dep
    csys = sc.coord_sys
    tag  = calc_input_statetag(c.var)

    # Convert the spacecraft state to the type required by the Calc
    return convert_orbitcalc_state(sc.state, csys, tag)
end

"""
    get_calc(c::OrbitCalc)

Evaluate an orbit-derived variable.

Arguments
- c::OrbitCalc: Container holding the spacecraft, variable tag, and optional dependency.

Returns
- Number or Vector{<:Real}: Computed value of the requested orbit variable.

# Notes:
- Required input state is determined by `calc_input_statetag(c.var)` and converted as needed.
- If a conversion requires μ, it is taken from `c.sc.coord_sys.origin` when available.

# Examples
```julia
sc = Spacecraft(
    state     = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
    time      = Time("2020-09-21T12:23:12", TAI(), ISOT()),
)
oc = OrbitCalc(sc, PositionVector())
val = get_calc(oc)
```
"""
@inline function get_calc(c::OrbitCalc)
    # Compute the state type required for the variable
    st = _state_for_calc(c)
    # Evaluate the variable 
    return _evaluate(c.var, st) 
end

"""
    set_calc!(c::AbstractCalc, newval)

Write a new value through a settable legacy calculation.

# Returns
The value returned by the calculation-specific setter. Throws an error when
the variable is read-only.

# Example
```julia
calc = BodyCalc(earth, GravParam())
set_calc!(calc, 398600.0)
```
"""
function set_calc!(c::AbstractCalc, newval)
    # Check if calc is settable
    if !calc_is_settable(c.var)
        error("Variable $(typeof(c.var)) is not settable.")
    end
    
    # Dispatch to type-specific implementation
    return _set_calc_type!(c, newval)
end

"""
    set_calc!(c::OrbitCalc, newval::Vector{<:Real})

Assign a new value to a settable orbit-derived variable.

Arguments
- c::OrbitCalc: Container holding the spacecraft and variable tag.
- newval::Vector{<:Real}: New value(s) for the variable.

# Notes:
- Required input state is determined by `calc_input_statetag(c.var)` and converted as needed.
- After assignment, the spacecraft’s state is converted back to its original type and updated.
- Non-settable variables throw an error (see `calc_is_settable`).

# Returns
- OrbitState: The updated spacecraft state.

# Examples
```julia 
sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
    time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
)
oc = OrbitCalc(sc, PositionVector())
set_calc!(oc, [7000.0, 300.0, 0.0])
```
"""
@inline function _set_calc_type!(c::OrbitCalc, newval::Vector{<:Real})
    # Compute the state type required for the variable
    st = _state_for_calc(c)
    # Update the state with the new value
    newst = _set!(c.var, st, newval) 
    # Convert the state back to the original type of the spacecraft state
    back = convert_orbitcalc_state(newst, c.sc.coord_sys, c.sc.state.statetype)
    # Constuct an OrbitState and set it back to the spacecraft
    c.sc.state = OrbitState(to_vector(back), c.sc.state.statetype)
    return c.sc.state
end

""" 
    _set_calc_type!(c::OrbitCalc, newval::Real)

Type-specific implementation for OrbitCalc scalar assignments.
Scalar values are forwarded as one-element vectors.
"""
@inline function _set_calc_type!(c::OrbitCalc, newval::Real)
    return _set_calc_type!(c, [newval])
end

"""
    _set_calc_type!(c::ManeuverCalc, vals::AbstractVector{<:Real})

Assign new value(s) to a settable maneuver-derived variable.

Arguments
- c::ManeuverCalc: Container holding the maneuver model, spacecraft, and variable tag.
- vals::AbstractVector{<:Real}: New value(s) for the variable.

# Notes:
- Delegates to the variable-specific `_set!(c, vals)` implementation.
- Settability is checked by the generic `set_calc!` method.

# Returns
- Returns nothing. Performs inplace update of Maneuver.
- Spacecraft state is not modified.  Use maneuver!() to apply maneuver effects.

# Examples
```julia
m  = ImpulsiveManeuver(axes=Inertial(), Isp=300.0, element1=0.01, element2=0.02, element3=-0.03)
sc = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                time  = Time("2020-01-01T00:00:00", TAI(), ISOT()))
mc = ManeuverCalc(m, sc, DeltaVVector())
set_calc!(mc, [0.01, 0.02, -0.03])

# output
```
"""
@inline _set_calc_type!(c::ManeuverCalc, vals::AbstractVector{<:Real}) = _set!(c, vals)



"""
    _set_calc_type!(c::AbstractCalc, newval)

Fallback implementation for calc types that don't have specific `_set_calc_type!` methods.
This should only be reached if settability checking is bypassed.
"""
function _set_calc_type!(c::AbstractCalc, newval)
    error("No _set_calc_type! implementation found for $(typeof(c)). " *
          "This calc type may need a specific implementation.")
end

"""
    calc_is_settable(variable::AbstractCalcVariable) -> Bool

Whether a legacy calculation variable supports assignment.

# Returns
`false` unless the variable type provides a settable method.

# Example
```jldoctest
calc_is_settable(SMA())

# output

true
```
"""
calc_is_settable(::AbstractCalcVariable) = false   # COV_EXCL_LINE


end
