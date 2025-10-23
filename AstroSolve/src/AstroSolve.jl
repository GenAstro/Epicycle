# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

### Function Registry and Evaluation System
__precompile__()
module AstroSolve

using NLsolve
using LinearAlgebra

using AstroBase
#using AstroProp
using AstroMan
#using AstroStates
using AstroFun

using AstroModels: Spacecraft, to_posvel, set_posvel!

# Public structs 
export SolverVariable, Event, Sequence, SequenceManager 

# Public functions
export set_sol_var, get_sol_var
export add_events!, topo_sort
export order_unique_vars, apply_event
export get_var_values, get_var_shifts, get_var_scales
export set_var_values, get_var_lower_bounds, get_var_upper_bounds 
export solver_fun!, is_astrosolve_stateful, reset_stateful_structs!, trajectory_solve!
export get_fun_values, get_fun_upper_bounds, get_fun_lower_bounds

# trait for structs that can be used in solving iterations.  
# stateful means the structs must be reset after each solve iteration.
#is_astrosolve_stateful(::Type)
is_astrosolve_stateful(::Type) = false  # default
is_astrosolve_stateful(::Type{T}) where {T<:Spacecraft} = true
is_astrosolve_stateful(::Type{T}) where {T<:ImpulsiveManeuver} = true

"""
    Event

Defines an event node for the trajectory sequence.
- `name`: A human-readable name or identifier for the event.
- `event`: The function closure to execute when this node runs.
- `vars`: Vector of variable objects (`SolverVariable`) or similar.
- `funcs`: Vector of functions (constraints/objectives) associated with this event.
"""
struct Event
    name::String
    event::Function
    vars::Vector{Any}
    funcs::Vector{Any}
end

# Outer constructor for Event with kwargs for user
Event(; name::String = "", event::Function = () -> nothing, vars = [], funcs = []) =
    Event(name, event, vars, funcs)

"""
    SolverVariable(; calc, lower_bound, upper_bound, shift, scale, name)

Represents a solver-controlled variable defined by a calc container (e.g., OrbitCalc, ManeuverCalc, BodyCalc).

Fields
- calc::AbstractCalcVariable     # the calc container (object + variable tag [+ context])
- numvars::Int                   # number of scalar variables for this calc (from calc variable tag)
- lower_bound::Vector{Float64}   # length numvars (defaults to -Inf)
- upper_bound::Vector{Float64}   # length numvars (defaults to +Inf)
- shift::Vector{Float64}         # length numvars (defaults to 0.0)
- scale::Vector{Float64}         # length numvars (defaults to 1.0)
- name::String                   # optional label

Notes
- numvars is derived from the calc’s variable tag (assumes a `var` field on the calc).
- Bounds/shift/scale are stored as Float64 for optimizer compatibility.
"""

mutable struct SolverVariable{C<:AbstractCalc}
    calc::C
    numvars::Int
    lower_bound::Vector{Float64}
    upper_bound::Vector{Float64}
    shift::Vector{Float64}
    scale::Vector{Float64}
    name::String

    function SolverVariable{C}(; 
        calc::C,
        lower_bound = nothing, 
        upper_bound = nothing, 
        shift = nothing, 
        scale = nothing, 
        name::String = ""
    ) where {C<:AbstractCalc}
        var = getfield(calc, :var)
        n = calc_numvars(var)

        to_vec(x, n, default) = x === nothing ? fill(default, n) :
                                (x isa AbstractVector ? Float64.(x) : fill(Float64(x), n))

        lb = to_vec(lower_bound, n, -Inf)
        ub = to_vec(upper_bound, n, +Inf)
        sh = to_vec(shift,       n, 0.0)
        sc = to_vec(scale,       n, 1.0)

        (length(lb) == n && length(ub) == n && length(sh) == n && length(sc) == n) ||
            throw(ArgumentError("SolverVariable: all bound/shift/scale vectors must have length $n."))

        new{C}(calc, n, lb, ub, sh, sc, name)
    end
end

# Convenience ctor infers C
SolverVariable(; calc::C, lower_bound=nothing, upper_bound=nothing, shift=nothing, scale=nothing, name::String="") where {C<:AbstractCalc} =
    SolverVariable{C}(; calc, lower_bound, upper_bound, shift, scale, name)


function Base.show(io::IO, sv::SolverVariable)
    calc_ty = nameof(typeof(sv.calc))
    var_str = hasfield(typeof(sv.calc), :var) ? string(getfield(sv.calc, :var)) : "?"
    println(io, "SolverVariable($calc_ty; var=$var_str)")
    println(io, "  numvars:     ", sv.numvars)
    println(io, "  lower_bound: ", sv.lower_bound)
    println(io, "  upper_bound: ", sv.upper_bound)
    println(io, "  shift:       ", sv.shift)
    println(io, "  scale:       ", sv.scale)
    println(io, "  name:        ", sv.name)
end

function set_sol_var(var::SolverVariable,val::Vector)
    # Test this calc is settable
    if !calc_is_settable(var.calc.var)
        throw(ArgumentError("set_sol_var: calc type $(typeof(var.calc)) does not support setting values."))
    end

    # Delegate to AstroFun; accept vectors for both scalar and vector calcs
    n = var.numvars
    length(val) == n || throw(ArgumentError("set_sol_var: expected length $n (got $(length(val)))."))
    if n == 1
        set_calc!(var.calc, val[1])
    else
        set_calc!(var.calc, val)
    end
    return var
end

function get_sol_var(var::SolverVariable)
    # Always return a Vector for SequenceManager
    vals = get_calc(var.calc)
    return vals isa AbstractVector ? vals : [vals]
end


"""
    apply_event(event::Event)

Execute the event's function closure.
"""
function apply_event(event::Event)
    event.event()
    return nothing
end

include("sequence.jl")

end