# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

__precompile__()

"""
    module AstroProp

Propagation interfaces, stopping conditions, and the interface to the ODE solvers.
"""
module AstroProp

# The solvers we name, and nothing else. OrdinaryDiffEq the metapackage pulls 176 packages
# against 91 for these, so eighty-five packages of stiff and specialist solvers used to load on
# every `using` and were never called. SciMLBase carries ODEProblem and the callbacks;
# CommonSolve carries `solve`.
using SciMLBase, CommonSolve
using OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner
using LinearAlgebra
using ForwardDiff

using EpicycleBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroFrames
using AstroModels: Spacecraft, to_posvel, set_posvel!, total_mass, SphericalDrag, AbstractDragGeometry, SphericalSRP, AbstractSRPGeometry
using AstroModels: HistorySegment, SpacecraftHistory, push_segment!
using AstroCallbacks: OrbitCalc, get_calc, Calc
import AstroCallbacks

using SatelliteToolboxTransformations: r_eci_to_ecef, ecef_to_geodetic, fetch_iers_eop,
                                       J2000, ITRF, DCM
import SatelliteToolboxBase: EARTH_ANGULAR_SPEED
using StaticArrays: SVector

import AstroCallbacks: AbstractFun, AbstractCalcVariable, AbstractOrbitVar
import AstroUniverse: translate
import AstroUniverse: Mu                   # re-export below
import AstroModels: Mass                   # re-export below
import EpicycleBase: AbstractStateTag, AbstractParamTag, AbstractVarTag
import EpicycleBase: ModelVariable, get_field, set_field!
# state_jac! and param_jac! are EpicycleBase's extension interface, declared there as
# bare generics. They have to be imported to be extended: defined at this module's top
# level after a plain `using`, each definition creates a NEW AstroProp function, and the
# eight force-model methods below hang off it where nobody dispatching on the interface
# can see them.
import EpicycleBase: state_jac!, param_jac!
import CommonSolve: solve                  # extended below for OrbitODEProblem

export IntegratorConfig
export Tsit5, Vern7, Vern9                 # the integrators AstroProp loads, so `using AstroProp` runs
export propagate!, ForceModel
export PosVel
export PointMassGravity, accel_eval!
export HarmonicGravity, AtmosphericDrag, SolarRadiationPressure
export Zonal, Exponential
export DualCone
export AbstractGeopotential, AbstractDensityModel, density
export geopotential_accel, geopotential_data, max_degree, max_order
export gravity_center, includes_central
export SphericalDrag, SphericalSRP, total_mass

export OrbitODE
export OrbitPropagator, StopAt
export PropDurationSeconds, PropDurationDays

export Mu, Mass                            # pass-through from AstroUniverse / AstroModels
export JacobianConfig, JacobianResult
export eval_jacobian!, eval_jacobian, state_jac!, param_jac!
export fd_differentiate_wrt
export STMConfig, OrbitODEProblem, PropagationResult, solve

"""
    PosVel <: AbstractStateTag

Tag singleton identifying the 6-element position–velocity state on a `Spacecraft`.
Use with `get_field` / `set_field!` and as the `tag` for state-vector solve-fors
and sensitivity blocks.

# Identity
`(sc, PosVel())` is the canonical `(obj, tag)` pair for the propagated
spacecraft state.  `get_field(sc, PosVel())` returns `[r..., v...]` (length 6);
`set_field!(sc, PosVel(), x)` writes back through `set_posvel!`.
"""
struct PosVel <: AbstractStateTag end

"""Return the position–velocity state of `sc` as a length-6 vector `[r..., v...]`."""
get_field(sc::Spacecraft, ::PosVel) = to_posvel(sc)

"""Write the position–velocity state of `sc` from a length-6 vector `[r..., v...]`."""
function set_field!(sc::Spacecraft, ::PosVel, x::AbstractVector)
    set_posvel!(sc, x)
    return nothing
end

"""
    OrbitODE

Abstract type for a force that acts on a spacecraft during propagation.

A force of your own subtypes `OrbitODE` and adds one method of [`accel_eval!`](@ref). A
`ForceModel` then sums it with the built-in forces for `propagate!`, for a collocation phase flown
by a `ForceModel`, and for the state transition matrix and sensitivity propagation.
"""
abstract type OrbitODE <: AbstractFun end

"""
    accel_eval!(force, t, y, dy, sc, params) -> dy

Write one force's contribution to the time derivative of a spacecraft's state.

# Arguments
- `force`: The force, a subtype of [`OrbitODE`](@ref).
- `t::Time`: The epoch of the evaluation.
- `y`: The state: position and velocity in km and km/s, followed by any further components the
  caller integrates, such as mass in a collocation phase.
- `dy`: The buffer this force writes into, the same length as `y`.
- `sc`: The spacecraft the force acts on, which carries its mass and drag and SRP geometry.
- `params`: Passed through by the caller. A collocation phase flown by a `ForceModel` passes
  `(control = u, params = p)`; `propagate!` passes internal bookkeeping a force does not read.

# Returns
`dy`.

# Notes
Every caller hands each force a buffer of zeros and sums the forces afterwards, so a method may
assign its contribution or add it, with the same result. Rows 4 to 6 are this force's acceleration
in km/s². Rows past 6 are summed across forces in the same way, so two thrusters that each write
a mass flow rate both count. Rows 1 to 3, ṙ = v, are written once by the caller, and whatever a
force writes there is ignored.

A method should accept `y` and `dy` of any element type, not only `Float64`: the Jacobian of a
force with no analytic `state_jac!` is taken by automatic differentiation through this method.
"""
function accel_eval! end


# Common supertype for gravity force models recognised by `_find_center`.
# Defined here so `PointMassGravity` can subtype it without pulling in the
# (currently blocked) external-force adapter.
abstract type AbstractGravityForce <: OrbitODE end

# external_force.jl is intentionally not included until the AstroForceModels,
# SatelliteToolboxGravityModels and ForwardDiff version conflict is resolved upstream.
# The file remains in the source tree.
include("point_mass_gravity.jl")
include("harmonic_gravity.jl")
include("zonal_gravity.jl")
include("atmospheric_drag.jl")
include("exponential_atmosphere.jl")
include("spherical_srp.jl")

"""
    IntegratorConfig

Configuration parameters for orbital integration using DifferentialEquations.jl solvers.

# Fields
- `integrator::Any`: ODE solver algorithm (e.g., `Vern7()`, `Tsit5()`)
- `dt::Union{Nothing, Float64}`: Fixed step size in seconds, or `Nothing` for adaptive stepping
- `reltol::Float64`: Relative tolerance for integration accuracy
- `abstol::Float64`: Absolute tolerance for integration accuracy

# Example
```julia
integ = IntegratorConfig(Vern7(); dt=3600.0, abstol = 1e-12, reltol=1e-12)
```
"""
struct IntegratorConfig
    integrator::Any                         
    dt::Union{Nothing, Float64}             
    reltol::Float64                       
    abstol::Float64                        

    function IntegratorConfig(
        integrator;
        dt = 5000.0,
        reltol = 1e-9,
        abstol = 1e-9
    )
        new(integrator, dt, reltol, abstol)
    end
end

"""
    ForceModel{N} <: OrbitODE

A collection of orbital dynamics forces and perturbations for spacecraft propagation.

# Fields
- `forces::NTuple{N, OrbitODE}`: Tuple of force models (e.g., gravity, drag, solar radiation pressure)
- `center::Union{CelestialBody, Nothing}`: Central gravitational body, determined automatically from forces

# Constructor
    ForceModel(forces...)
    ForceModel(force_tuple)

The central body is automatically determined from the primary gravitational force in the model.

# Example
```julia
gravity = PointMassGravity(earth, ())
model = ForceModel(gravity)
```
"""
struct ForceModel{N} <: OrbitODE
    forces::NTuple{N, OrbitODE}
    center::Union{CelestialBody, Nothing}
end

include("orbit_propagator.jl")
include("jacobian_config.jl")
include("orbit_ode_problem.jl")
include("variational.jl")

"""
    ForceModel(forces...)
    ForceModel(force_tuple)

Create a force model for orbital propagation from one or more force components.

# Arguments
- `forces...`: Variable number of force objects (e.g., `PointMassGravity`)
- `force_tuple`: Tuple of force objects

# Returns
- `ForceModel{N}`: Force model with N force components

The central gravitational body is automatically determined from the primary 
gravitational force in the collection.

# Examples
```julia
gravity = PointMassGravity(earth)
model = ForceModel(gravity)
```
"""
function ForceModel(forces::Tuple{Vararg{OrbitODE}})
    center = _find_center(forces)
    return ForceModel{length(forces)}(forces, center)
end

ForceModel(force::T) where {T<:OrbitODE} = ForceModel((force,))

# Varargs form so heterogeneous forces compose as `ForceModel(gravity, drag, ...)`.
ForceModel(forces::OrbitODE...) = ForceModel(forces)

"""
    gravity_center(force) -> Union{CelestialBody, Nothing}

The central body a gravity force is referenced to, or `nothing` if the force is not central-body
gravity (drag, SRP, …). Forces implement this so `ForceModel` can check that every gravity force
agrees on one central body. A new gravity model participates in that check by adding a method here.
"""
gravity_center(::OrbitODE)             = nothing
gravity_center(f::PointMassGravity)    = f.central_body
gravity_center(f::HarmonicGravity)     = f.central_body
gravity_center(f::AtmosphericDrag)     = f.central_body
gravity_center(f::SolarRadiationPressure) = f.central_body

"""
    includes_central(force) -> Bool

Whether a gravity force adds the central body's own gravity (the monopole term). `ForceModel` uses
this to ensure exactly one model provides the central term for a given body. Defaults to `false`.
"""
includes_central(::OrbitODE)          = false
includes_central(f::PointMassGravity) = f.include_center
includes_central(f::HarmonicGravity)  = true

function _find_center(forces::Tuple)
    centers      = CelestialBody[]   # every gravity force's reference body
    with_central = CelestialBody[]   # the forces that add the central monopole
    for f in forces
        c = gravity_center(f)
        c === nothing && continue
        push!(centers, c)
        includes_central(f) && push!(with_central, c)
    end
    # FR-FORCE-16: only one force may add the central term for a given body.
    for i in eachindex(with_central), j in (i + 1):lastindex(with_central)
        with_central[i] === with_central[j] &&
            error("ForceModel: two gravity forces both include the central term for " *
                  "'$(with_central[i].name)'. Set `include_center = false` on the point-mass " *
                  "force so only one model provides the central gravity for that body.")
    end
    isempty(centers) && return nothing
    all(c -> c === centers[1], centers) ||
        error("ForceModel: the gravity forces reference different central bodies.")
    return centers[1]
end

# The spacecraft and forces of one propagation, as the engine below consumes them. Internal:
# `propagate!(::OrbitPropagator, ...)` builds one per call.
struct DynSys
    spacecraft::Vector{<:Spacecraft}
    forces::OrbitODE

    function DynSys(; spacecraft::Vector{<:Spacecraft},
                      forces::OrbitODE,
                      )
        return new(spacecraft, forces)
    end
end

function _build_odereg(spacecraft_list::Vector{<:Spacecraft})
    reg = Dict{Spacecraft, Dict{Symbol, UnitRange{Int}}}()
    index = 1
    for sc in spacecraft_list
        reg[sc] = Dict(:posvel => index:(index+5))
        index += 6
    end
    return reg
end

function _build_odes!(model::ForceModel, start_epoch, du, u, p, t, spacecraft_list::Vector{<:Spacecraft})
    odereg = p[:odereg]
    for sc in spacecraft_list
        idxs = odereg[sc][:posvel]
        posvel = u[idxs[1:6]]

        current_time = start_epoch + t/86400.0
        acc = zeros(eltype(posvel), 6)
        a_sum = zeros(eltype(posvel), 3)
        for force in model.forces
            # A fresh buffer per force. Built-in forces overwrite what they write and a user force
            # may add to it; without the reset an additive force counted every earlier force again.
            fill!(acc, zero(eltype(acc)))
            accel_eval!(force, current_time, posvel, acc, sc, p)
            # Superimpose acceleration contributions only; kinematics set once below.
            a_sum[1] += acc[4]; a_sum[2] += acc[5]; a_sum[3] += acc[6]
        end
        du[idxs[1]] = posvel[4]; du[idxs[2]] = posvel[5]; du[idxs[3]] = posvel[6]
        du[idxs[4]] = a_sum[1];  du[idxs[5]] = a_sum[2];  du[idxs[6]] = a_sum[3]
    end
end

function _build_state(model::ForceModel, spacecraft_list::Vector{<:Spacecraft}, odereg::Dict)
    max_index = maximum([maximum(v[:posvel]) for v in values(odereg)])

    #first_state = CartesianState(spacecraft_list[1].state, model.center.mu).posvel
    first_state = to_posvel(spacecraft_list[1])
    state_vector = zeros(eltype(first_state), max_index)
    for sc in spacecraft_list
        idxs = odereg[sc][:posvel]
        #state_vector[idxs] .= CartesianState(sc.state, model.center.mu).posvel
        state_vector[idxs] .= to_posvel(sc)
    end
    return state_vector
end

function _update_structs!(forces::ForceModel, sol_u::Vector{<:Real}, odereg::Dict, sol_t::Real = 0.0, full_sol::Union{Nothing,ODESolution} = nothing)
    for (sc, idx_map) in odereg
        # Update dynamic state
        if :posvel in keys(idx_map)
            # TODO: Fix this to handle type via OrbitState util
            #sc.state = CartesianState(sol_u[idx_map[:posvel]])
            set_posvel!(sc, sol_u[idx_map[:posvel]])
        end

        # Append to history
        # Use TT for Earth-centered, TDB for others (matches propagation time scale)
        center_body = forces.center
        time_scale = (center_body === earth) ? TT() : TDB()
        initialtime = (center_body === earth) ? sc.time.tt.jd : sc.time.tdb.jd
        
        if full_sol !== nothing && sc.save_history
            # Extract times and states from solution
            times = [Time(initialtime + t / 86400.0, time_scale, JD()) for t in full_sol.t]
            states = [CartesianState(copy(u[idx_map[:posvel]])) for u in full_sol.u]
            
            # Create HistorySegment and add to spacecraft history
            segment = HistorySegment(times, states, sc.coord_sys, name="propagate")
            push_segment!(sc.history, segment)
        end

        # Update current time using appropriate scale
        finaltime = initialtime + sol_t/86400.0
        sc.time = Time(finaltime, time_scale, JD())
    end
end

# The integration engine behind `propagate!(::OrbitPropagator, ...)`: integrates `model` with
# `config`, stopping on the callbacks and at most one time-based condition in `stop_conditions`.
function _propagate_dynsys!(model::DynSys, config::IntegratorConfig,
    stop_conditions...;
    direction::Symbol = :forward,
    kwargs...)

    # Validate direction keyword
    direction in (:forward, :backward, :infer) ||
        error("Invalid direction: $direction. Must be :forward, :backward, or :infer.")

    # Separate time-based from state-based stopping conditions
    time_conds = filter(_is_time_condition, stop_conditions)
    state_conds = filter(!_is_time_condition, stop_conditions)

    # Validate at most one time-based condition
    if length(time_conds) > 1
        error("Multiple time-based stopping conditions not allowed. Found $(length(time_conds)) conditions.")
    end

    # Infer direction if requested
    actual_direction = if direction == :infer
        if isempty(time_conds)
            :forward  # Default for state-based or no conditions
        else
            _infer_direction(time_conds[1])
        end
    else
        direction
    end

    # Validate explicit direction matches duration sign
    # Duration sign is semantically meaningful: positive = forward, negative = backward
    if direction != :infer && !isempty(time_conds)
        target = time_conds[1].target
        # Check for contradictions between sign and explicit direction
        if target >= 0 && actual_direction == :backward
            error("Duration is positive (forward) but explicit direction is :backward. Use negative duration or direction=:infer.")
        elseif target < 0 && actual_direction == :forward
            error("Duration is negative (backward) but explicit direction is :forward. Use positive duration or direction=:infer.")
        end
    end

    # Compute tspan from time condition or use default
    tf = if isempty(time_conds)
        actual_direction == :forward ? 1.0e12 : -1.0e12
    else
        _compute_tf(time_conds[1], actual_direction)
    end
    tspan = (0.0, tf)

    # Build callbacks only from state-based conditions
    callbackset = isempty(state_conds) ? nothing :
       length(state_conds) == 1 ? state_conds[1] :
       CallbackSet(state_conds...)

    odereg = _build_odereg(model.spacecraft)
    params = (forces = model.forces, odereg = odereg)
    state0 = _build_state(model.forces, model.spacecraft, odereg)
    
    # Use TT for Earth-centered dynamics, TDB for others
    center_body = model.forces.center
    start_epoch = (center_body === earth) ? model.spacecraft[1].time.tt : model.spacecraft[1].time.tdb

    actual_direction == :forward || actual_direction == :backward ||
        error("Unknown direction: $actual_direction. Use :forward or :backward.")

    prob = ODEProblem((du, u, p, t) -> _build_odes!(model.forces, start_epoch,
         du, u, p, t, model.spacecraft), state0, tspan, params)

    sol = solve(prob, config.integrator;
    callback=callbackset,
    dt=config.dt,
    reltol=config.reltol,
    abstol=config.abstol,
    kwargs...)

    _update_structs!(model.forces, sol.u[end], odereg, sol.t[end], sol)

    return sol
end

# Helper: infer propagation direction from time-based stop condition
function _infer_direction(stop::StopAt{<:Any, <:IntegratorTimeCalc})
    # Infer from sign of duration
    target = stop.target
    return target >= 0 ? :forward : :backward
end

# Helper: compute final time from time-based stopping condition
function _compute_tf(stop::StopAt{<:Any, <:IntegratorTimeCalc}, direction::Symbol)
    var = stop.var
    target = stop.target
    stop_dir = stop.direction
    
    # Compute elapsed time in seconds (magnitude)
    tf_magnitude = if var isa PropDurationSeconds
        abs(Float64(target))
    elseif var isa PropDurationDays
        abs(Float64(target)) * 86400.0
    else
        error("Unknown IntegratorTimeCalc type: $(typeof(var))")
    end
    
    # Validate target is non-zero
    if tf_magnitude == 0.0
        error("Time-based stopping condition duration must be non-zero")
    end
    
    # Validate direction compatibility with stop_dir parameter (if specified)
    if stop_dir != 0
        if stop_dir > 0 && direction == :backward
            error("StopAt direction parameter is positive (increasing) but propagation direction is :backward")
        elseif stop_dir < 0 && direction == :forward
            error("StopAt direction parameter is negative (decreasing) but propagation direction is :forward")
        end
    end
    
    # Apply sign based on propagation direction
    return direction == :backward ? -tf_magnitude : tf_magnitude
end

include("precompile.jl")

end
