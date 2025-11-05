# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()

"""
    module AstroProp

Propagation interfaces, stopping conditions, and interface to OrdinaryDiffEq.
"""
module AstroProp

using OrdinaryDiffEq, LinearAlgebra

using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroCoords
using AstroModels: Spacecraft, to_posvel, set_posvel!
import AstroModels: push_history_segment!
using AstroFun: OrbitCalc, get_calc

import AstroFun: AbstractFun, AbstractCalcVariable, AbstractOrbitVar
import AstroUniverse: translate
export TwoBodyGravity, ExponentialAtmosphere, CartesianODE, IntegratorConfig
export DynamicsSystem, propagate, DynSys, ForceModel
export StopAtApoapsis, StopAtAscendingNode, StopAtPeriapsis, StopAtDays
export PosVel
export StopAtSeconds, StopAtRadius
export nbody_perts
export PointMassGravity, compute_point_mass_gravity!, evaluate, accel_eval!

export OrbitODE
export IntegratorConfig
export OrbitPropagator, StopAt

"""
    PosVel <: AbstractState

Position-velocity state representation for orbital propagation. This is a legacy 
state type used internally by the propagation system.

# Fields
- `numvars::Int`: Number of state variables (always 6 for position and velocity)

# Notes
This struct is marked for deprecation and should be replaced with new state 
representations from AstroStates.
"""
struct PosVel <: AbstractState 
    numvars::Int
    PosVel() = new(6)
end

abstract type OrbitODE <: AbstractFun end

include("point_mass_gravity.jl")
include("stop_conditions.jl")

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
function ForceModel(forces::Tuple{Vararg{T}}) where {T<:OrbitODE}
    center = _find_center(forces)
    return ForceModel{length(forces)}(forces, center)
end

ForceModel(force::T) where {T<:OrbitODE} = ForceModel((force,))

function _find_center(forces::Tuple)
    centers = CelestialBody[]
    for f in forces
        if f isa PointMassGravity
            push!(centers, f.central_body)
        end
    end
    if isempty(centers)
        return nothing
    elseif length(unique(centers)) == 1
        return first(centers)
    else
        error("Multiple conflicting central bodies found in ForceModel.")
    end
end

"""
    DynSys(; spacecraft, forces)

Container for a dynamic system segment.

# Keyword Arguments
- `spacecraft`: Vector of spacecraft structs (subtypes of `Spacecraft`)
- `forces`: OrbitODE model (e.g. CartesianForceModel)

# Fields
- `spacecraft::Vector{<:Spacecraft}`: Collection of spacecraft to propagate
- `forces::OrbitODE`: Force model defining the dynamics

# Constructor
    DynSys(; spacecraft, forces)

# Example
```julia
sc = Spacecraft()
gravity = PointMassGravity(earth, ())
sys = DynSys(spacecraft=[sc], forces=gravity)
```

# Notes
This interface will be deprecated in future releases
"""
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
        acc_sum = zeros(eltype(posvel), 6);
        for force in model.forces
            accel_eval!(force, current_time, posvel, acc, sc, p)
            acc_sum .+= acc
        end
        du[idxs[1:6]] .= acc
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
        initialtime = sc.time.jd  # This returns Float64 (days)
        if full_sol !== nothing
            new_segment = [
                (Time(initialtime+ t / 86400.0, TDB(), JD()), copy(u[idx_map[:posvel]]))
                for (t, u) in zip(full_sol.t, full_sol.u)
            ]
            push_history_segment!(sc, new_segment)
        end

        # Update current time (TAI)
        finaltime = initialtime + sol_t/86400.0
        sc.time = Time(finaltime, TDB(), JD())
    end
end

"""
    This interface will be deprecated in future releases.

    propagate(model::DynSys, config::IntegratorConfig, stop_conditions...; 
              direction=:forward, prop_stm=false, kwargs...)

Propagate spacecraft orbital states using numerical integration.
"""
function propagate(model::DynSys, config::IntegratorConfig,
    stop_conditions...;
    direction::Symbol = :forward, prop_stm::Bool = false,
    kwargs...)

    callbackset = isempty(stop_conditions) ? nothing :
       length(stop_conditions) == 1 ? stop_conditions[1] :
       CallbackSet(stop_conditions...)

    odereg = _build_odereg(model.spacecraft)
    params = (forces = model.forces, odereg = odereg)
    state0 = _build_state(model.forces, model.spacecraft, odereg)
    start_epoch = model.spacecraft[1].time.tai

    tspan = direction == :forward ? (0.0, 1.0e12) :
    direction == :backward ? (0.0, -1.0e12) :
    error("Unknown direction: $direction. Use :forward or :backward.")

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

end