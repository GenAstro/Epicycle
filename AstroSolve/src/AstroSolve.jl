# Copyright (C) 2025 Gen Astro LLC
SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()
"""
    module AstroSolve

Solver interfaces, event sequencing, solver variables and constraints 
and trajectory optimization utilities.  The module implements a Directed
Acyclic Graph for trajectory event sequencing.
"""
module AstroSolve

using NLsolve
using SNOW
using LinearAlgebra
using Printf

using AstroBase
using AstroMan
using AstroFun
using AstroModels: Spacecraft, to_posvel, set_posvel!

export SolverVariable, Event, Sequence, SequenceManager 

export set_sol_var, get_sol_var
export add_events!, topo_sort
export order_unique_vars, apply_event
export get_var_values, get_var_shifts, get_var_scales
export set_var_values, get_var_lower_bounds, get_var_upper_bounds 
export solver_fun!, is_astrosolve_stateful, reset_stateful_structs!, trajectory_solve!
export get_fun_values, get_fun_upper_bounds, get_fun_lower_bounds
export trajectory_solve, sequence_report, solution_report

# if stateful, struct must be reset after each solve iteration.
is_astrosolve_stateful(::Type) = false 
is_astrosolve_stateful(::Type{T}) where {T<:Spacecraft} = true
is_astrosolve_stateful(::Type{T}) where {T<:ImpulsiveManeuver} = true

"""
    Event

Defines an event node for trajectory sequence execution in mission design workflows.

# Fields
- `name::String`: Human-readable identifier for the event
- `event::Function`: Function closure to execute when this node runs
- `vars::Vector{Any}`: Vector of solver variables (`SolverVariable`) associated with this event
- `funcs::Vector{Any}`: Vector of constraint/objective functions for this event

# Notes
Events are the fundamental building blocks of trajectory sequences, encapsulating
both the action to be performed and the optimization variables/constraints
associated with that action.

# Examples
```julia
# Simple propagation event
prop_event = Event(
    name = "Propagate to Apoapsis",
    event = () -> propagate(dynsys, integ, StopAtApoapsis(sat)),
)

# Event with solver variables and constraints
sc = Spacecraft(); toi = ImpulsiveManeuver();
var_dv = SolverVariable(
    calc = ManeuverCalc(toi, sc, DeltaVVector()),
)
fun_dv = Constraint(
    calc = ManeuverCalc(toi, sc, DeltaVMag()),
    lower_bounds = [0.1],
    upper_bounds = [0.1],
    scale = [1.0],
)
maneuver_event = Event(
    event = () -> maneuver!(sc, toi),
    vars = [var_dv],
    funcs = [fun_dv],
    name = "Departure Maneuver",
)
```
"""
struct Event
    name::String
    event::Function
    vars::Vector{Any}
    funcs::Vector{Any}
end

"""
    Base.show(io::IO, event::Event)

Display method for `Event` showing name, actual variable count, and function count.
"""
function Base.show(io::IO, event::Event)
    name_str = isempty(event.name) ? "<unnamed>" : "\"$(event.name)\""
    
    # Count actual scalar variables (not just SolverVariable structs)
    total_vars = 0
    for var in event.vars
        if hasfield(typeof(var), :numvars)
            total_vars += var.numvars
        else
            total_vars += 1  # Fallback for unknown types
        end
    end
    
    # Count actual constraint functions (not just Constraint structs)  
    total_funcs = 0
    for func in event.funcs
        if hasfield(typeof(func), :numfuncs)
            total_funcs += func.numfuncs
        elseif hasmethod(length, (typeof(func),))
            # Try to get length if it's a collection
            try
                total_funcs += length(func)
            catch
                total_funcs += 1
            end
        else
            total_funcs += 1  # Fallback for unknown types
        end
    end
    
    print(io, "Event($name_str; $total_vars vars, $total_funcs funcs)")
end 

"""
    Event(; name::String = "", event::Function = () -> nothing, vars = [], funcs = [])

Outer constructor for Event with kwargs.
"""
Event(; name::String = "", event::Function = () -> nothing, vars = [], funcs = []) =
    Event(name, event, vars, funcs)

"""
    SolverVariable{C<:AbstractCalc}

Represents a solver-controlled variable defined by a calculation container for 
trajectory optimization.

# Fields
- `calc::C`: Calculation container (e.g., `OrbitCalc`, `ManeuverCalc`, `BodyCalc`)
- `numvars::Int`: Number of scalar variables for this calculation
- `lower_bound::Vector{Float64}`: Lower bounds for optimization (length `numvars`, 
  defaults to `-Inf`)
- `upper_bound::Vector{Float64}`: Upper bounds for optimization (length `numvars`, 
  defaults to `+Inf`)
- `shift::Vector{Float64}`: Variable shifting for numerical conditioning 
  (length `numvars`, defaults to `0.0`)
- `scale::Vector{Float64}`: Variable scaling for numerical conditioning 
  (length `numvars`, defaults to `1.0`)
- `name::String`: Optional human-readable identifier

# Constructor
    SolverVariable(; calc, lower_bound=nothing, upper_bound=nothing, shift=nothing,
    scale=nothing, name="")

# Arguments
- `calc`: Calculation container implementing `AbstractCalc` interface
- `lower_bound`: Scalar or vector of lower bounds (broadcast to `numvars` if scalar)
- `upper_bound`: Scalar or vector of upper bounds (broadcast to `numvars` if scalar)  
- `shift`: Scalar or vector of shift values for numerical conditioning
- `scale`: Scalar or vector of scale factors for numerical conditioning
- `name`: Optional descriptive name for the variable

# Notes
- See `AbstractCalc` for supported calculation types
- `numvars` is automatically determined from the calculation variable type
- Bounds, shift, and scale vectors are stored as `Float64` for optimizer compatibility
- The calculation must support setting values (`calc_is_settable(calc.var) == true`) 

# Examples
```julia
using Epicycle

# Spacecraft and maneuver objects
sc = Spacecraft(); toi = ImpulsiveManeuver()

# 3D delta-V vector variable with individual component bounds
var_toi = SolverVariable(
    calc = ManeuverCalc(toi, sc, DeltaVVector()),
    lower_bound = [-2.0, -2.0, -2.0],
    upper_bound = [2.0, 2.0, 2.0],
    scale = [1.0, 1.0, 1.0],  
    name = "Departure DV Vector"
)

# Orbital element variable
var_sma = SolverVariable(
    calc = OrbitCalc(sc, SMA()),
    lower_bound = 6700.0,   
    upper_bound = 42000.0,  
    name = "Target Semi-Major Axis"
)
```
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

"""
    SolverVariable(; calc, lower_bound=nothing, upper_bound=nothing, shift=nothing, scale=nothing, name="")

Convenience constructor that automatically infers the calculation type parameter `C` from the `calc` argument.

This eliminates the need to explicitly specify the type parameter when creating `SolverVariable` instances.

# Usage Comparison
```julia
# Without convenience constructor (explicit type parameter):
var1 = SolverVariable{OrbitCalc}(calc=OrbitCalc(sat, SemiMajorAxis()))

# With convenience constructor (type inferred):
var2 = SolverVariable(calc=OrbitCalc(sat, SemiMajorAxis()))  # Preferred
```

This constructor delegates to `SolverVariable{C}(...)` with the inferred type.
"""
SolverVariable(; calc::C, lower_bound=nothing, upper_bound=nothing, shift=nothing, scale=nothing, name::String="") where {C<:AbstractCalc} =
    SolverVariable{C}(; calc, lower_bound, upper_bound, shift, scale, name)


""" 
   Base.show(io::IO, sv::SolverVariable)

Show method for `SolverVariable` that displays key information about the variable.
"""
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

"""
    set_sol_var(var::SolverVariable, val::Vector)

Set the value(s) of the solver variable struct
"""
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

"""
    get_sol_var(var::SolverVariable)

Get the solver variable values from the struct.
"""
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
include("sequence_report.jl")

end