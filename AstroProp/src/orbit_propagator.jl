# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    OrbitPropagator(forces::ForceModel, integ::IntegratorConfig)

Lightweight propagator bundling a force model and integrator configuration.
Use with `propagate(op, sc, stops...)` where `sc` is a `Spacecraft` or a
vector of spacecraft, and `stops` are generic `StopAt` conditions.
"""
struct OrbitPropagator
    forces::ForceModel
    integ::IntegratorConfig
end

#=
"""
    StopAt(calc, target; direction=0)

Generic, calc-based stopping condition. `calc` is typically an `AstroFun.OrbitCalc`,
`target` is the scalar crossing value, and `direction ∈ {-1, 0, +1}` controls
root direction (DifferentialEquations convention).
"""
struct StopAt{C,T}
    calc::C
    target::T
    direction::Int
    function StopAt(calc::C, target::T; direction::Int=0) where {C,T}
        return new{C,T}(calc, target, direction)
    end
end
=#

struct StopAt{S,V,T}
    subject::S                  # Spacecraft | Maneuver | CelestialBody | …
    var::V                      # e.g., VelMag(), PosMag(), GravParam()
    target::T                   # numeric: Real or Vector matching calc output
    direction::Int              # -1, 0, +1
end

# Positional target (required)
StopAt(subject, var, target; direction::Int=0) = StopAt(subject, var, target, direction)

# Extensibility point: how to build a Calc from (subject, var)
# Users/packages extend this for new subject/var kinds.
make_calc(subject, var) = error("make_calc not implemented for $(typeof(subject)), $(typeof(var))")

# Common case: orbit variables on spacecraft
make_calc(sc::Spacecraft, v::AbstractOrbitVar) = OrbitCalc(sc, v)

# Rebind swaps the subject to the promoted owner; var/target/direction unchanged
rebind(x::StopAt, owner_map::Dict{Any,Any}) =
    StopAt(get(owner_map, x.subject, x.subject), x.var, x.target, x.direction)

# Update a subject from the integrator state u (specialize per subject type)
_subject_update_from_u!(subject, dynsys, u) = error("No _subject_update_from_u! for $(typeof(subject))")

# Spacecraft specialization: load 6x1 pos/vel slice into the Spacecraft
_subject_update_from_u!(sc::Spacecraft, dynsys, u) = begin
    pv = _posvel_from_u(u, dynsys, sc)
    set_posvel!(sc, pv)
    nothing
end
# --- Utilities ---

# Normalize to a Vector{Spacecraft}
_as_scvec(sc::Spacecraft) = Spacecraft[sc]
_as_scvec(v::Vector{<:Spacecraft}) = v

#=
function _calc_subject(oc::OrbitCalc)
    # Try common property names first (fast path)
    for name in (:spacecraft, :obj, :subject)
        if Base.hasproperty(oc, name)
            val = getproperty(oc, name)
            if val isa Spacecraft
                return val
            end
        end
    end
    # Fallback: scan fields and return the first Spacecraft-typed field
    T = typeof(oc)
    for fname in fieldnames(T)
        val = getfield(oc, fname)
        if val isa Spacecraft
            return val
        end
    end
    error("StopAt: could not locate a Spacecraft inside OrbitCalc($(T)). Provide a subject accessor or standard field.")
end
=# 

# Find spacecraft index in a DynSys
_sc_index(dynsys, sc::Spacecraft) = findfirst(x -> x === sc, getfield(dynsys, :spacecraft))

# Extract a 6x1 Cartesian pos/vel slice for a spacecraft from integrator state u
function _posvel_from_u(u, dynsys, sc::Spacecraft)
    idx = _sc_index(dynsys, sc)
    idx === nothing && error("StopAt: spacecraft not found in DynamicsSystem")
    i0 = 6*(idx-1) + 1
    return collect(@view u[i0:i0+5])
end

# --- StopAt callback builders ---
# Build a ContinuousCallback for StopAt over an OrbitCalc
function _build_callback(cond::StopAt, dynsys)
    subject = cond.subject
    var     = cond.var
    target  = cond.target
    dir     = cond.direction

    # Build a calc from (subject, var) after we know the dynsys/subject we’ll mutate
    calc = make_calc(subject, var)

    # Event function: zero when calc hits target
    function g(u, t, integ)
        # Bring subject up-to-date from integrator state
        _subject_update_from_u!(subject, dynsys, u)
        # Evaluate current value via AstroFun
        val = get_calc(calc)
        return val - target
    end

    term!(integ) = terminate!(integ)

    if dir == 0
        return ContinuousCallback(g, term!; rootfind=true)
    elseif dir > 0
        return ContinuousCallback(g, term!; affect_neg! = (_integ)->nothing, rootfind=true)
    else
        return ContinuousCallback(g, (_integ)->nothing; affect_neg! = term!, rootfind=true)
    end
end

# Fallback for unsupported calc types
#function _build_callback(cond::StopAt, dynsys)
#    calcT = typeof(cond.calc)
#    error("StopAt: no callback builder for calc type $calcT. Add a method for _build_callback(::StopAt{$calcT}, dynsys).")
#end

# --- New propagate API ---

"""
    propagate(op::OrbitPropagator, sc_or_scs, stops...; direction=:forward, kwargs...)

Propagate one or more spacecraft with the force/integrator in `op` until any
of the provided `StopAt` conditions triggers. `stops` may be one or many
`StopAt` values; both a single `Spacecraft` or a vector of `Spacecraft` are accepted.

Returns the ODE solution from DifferentialEquations.jl (same as legacy API).
"""
function propagate(op::OrbitPropagator, sc_or_scs, stops...;
                   direction::Symbol = :forward, kwargs...)
    scv = _as_scvec(sc_or_scs)
    dyn = DynSys(spacecraft=scv, forces=op.forces)

    # Build callbacks (CallbackSet if multiple)
    callbacks = map(s -> _build_callback(s, dyn), collect(stops))
    cbset = isempty(callbacks) ? nothing :
            length(callbacks) == 1 ? callbacks[1] : CallbackSet(callbacks...)

    # Delegate to existing DynSys-based propagate (reuses your ODE assembly and updates)
    return propagate(dyn, op.integ, cbset; direction=direction, kwargs...)
end