# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

__precompile__()
module AstroProp

# TODO, trap error when user propagates a spacecraft in propagate() that is not in dynsys

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
#export STMStore, add_stm!, get_stm
#export StateTransitionMatrix, STMConfig,set_stm!, get_stm,  reset_stm!, has_stm, composable, assemble_stm
export IntegratorConfig
export OrbitPropagator, StopAt

#TODO: Old tags used in ForceModels.  Should replace with new tags like Cartesian()
struct PosVel <: AbstractState 
    numvars::Int
    PosVel() = new(6)
end



#export get_state

abstract type OrbitODE <: AbstractFun end

include("point_mass_gravity.jl")
include("stop_conditions.jl")

#export StopAt

struct IntegratorConfig
    integrator::Any                         # e.g., Vern7(), Tsit5()
    dt::Union{Nothing, Float64}             # Optional step size (Nothing = adaptive)
    reltol::Float64                         # Relative tolerance
    abstol::Float64                         # Absolute tolerance

    function IntegratorConfig(
        integrator;
        dt = 5000.0,
        reltol = 1e-9,
        abstol = 1e-9
    )
        new(integrator, dt, reltol, abstol)
    end
end

# ==========================
# State/ODE Assembly & Update & Registry (Dispatched on CartesianForceModel)
# ==========================

struct ForceModel{N} <: OrbitODE
    forces::NTuple{N, OrbitODE}
    center::Union{CelestialBody, Nothing}
end

include("orbit_propagator.jl")

function ForceModel(forces::Tuple{Vararg{T}}) where {T<:OrbitODE}
    center = find_center(forces)
    return ForceModel{length(forces)}(forces, center)
end

ForceModel(force::T) where {T<:OrbitODE} = ForceModel((force,))

function find_center(forces::Tuple)
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

# =======================
# DynSys Definition
# =======================
"""
    DynSys(; spacecraft, forces)

Container for a dynamic system segment.

# Keyword Arguments
- `spacecraft`: Vector of spacecraft structs (subtypes of `Spacecraft`)
- `forces`: OrbitODE model (e.g. CartesianForceModel)
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

function build_odereg(spacecraft_list::Vector{<:Spacecraft})
    reg = Dict{Spacecraft, Dict{Symbol, UnitRange{Int}}}()
    index = 1
    for sc in spacecraft_list
        reg[sc] = Dict(:posvel => index:(index+5))
        index += 6
    end
    return reg
end

function build_odes!(model::ForceModel, start_epoch, du, u, p, t, spacecraft_list::Vector{<:Spacecraft})
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

function build_state(model::ForceModel, spacecraft_list::Vector{<:Spacecraft}, odereg::Dict)
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

function update_structs!(forces::ForceModel, sol_u::Vector{<:Real}, odereg::Dict, sol_t::Real = 0.0, full_sol::Union{Nothing,ODESolution} = nothing)
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

function propagate(model::DynSys, config::IntegratorConfig,
    stop_conditions...;
    direction::Symbol = :forward, prop_stm::Bool = false,
    kwargs...)

    callbackset = isempty(stop_conditions) ? nothing :
       length(stop_conditions) == 1 ? stop_conditions[1] :
       CallbackSet(stop_conditions...)

    odereg = build_odereg(model.spacecraft)
    params = (forces = model.forces, odereg = odereg)
    state0 = build_state(model.forces, model.spacecraft, odereg)
    start_epoch = model.spacecraft[1].time.tai

    tspan = direction == :forward ? (0.0, 1.0e12) :
    direction == :backward ? (0.0, -1.0e12) :
    error("Unknown direction: $direction. Use :forward or :backward.")

    prob = ODEProblem((du, u, p, t) -> build_odes!(model.forces, start_epoch, du, u, p, t, model.spacecraft),
           state0, tspan, params)
  

    sol = solve(prob, config.integrator;
    callback=callbackset,
    dt=config.dt,
    reltol=config.reltol,
    abstol=config.abstol,
    kwargs...)

    update_structs!(model.forces, sol.u[end], odereg, sol.t[end], sol)

    return sol
end

end