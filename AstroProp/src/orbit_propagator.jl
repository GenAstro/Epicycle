# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    OrbitPropagator

Orbital propagator configuration combining force models and integration settings.

# Fields
- `forces::ForceModel`: Force model defining the orbital dynamics
- `integ::IntegratorConfig`: Integration configuration (solver, tolerances, step size)

# Examples
```julia
gravity = PointMassGravity(earth,(moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)
```
"""
struct OrbitPropagator
    forces::ForceModel
    integ::IntegratorConfig
end

"""
    StopAt{S,V<:AbstractOrbitVar,T}

Generic stopping condition for orbital propagation based on calculated quantities.

# Fields
- `subject::S`: The object to evaluate (e.g., `Spacecraft`, `Maneuver`, `CelestialBody`)
- `var::V`: The calculation variable to monitor (must be <: AbstractOrbitVar, e.g., `PosX()`, `VelMag()`)
- `target::T`: Target value to stop at (numeric value or vector matching calc output)
- `direction::Int`: Crossing direction (-1: decreasing, 0: any, +1: increasing)

# Constructor
    StopAt(subject, var, target; direction::Int=0)

# Examples
```julia
# Stop when spacecraft reaches ascending node (z-position = 0, increasing)
sc = Spacecraft(state=CartesianState([7000.0, 0, 0, 0, 7.5, 0]))
stop_cond = StopAt(sc, PosZ(), 0.0; direction=+1)

# Stop at apoapsis (position dot velocity = 0, decreasing)  
stop_apo = StopAt(sc, PosDotVel(), 0.0; direction=-1)
```
"""
struct StopAt{S,V<:AbstractOrbitVar,T}
    subject::S                 
    var::V                     
    target::T                  
    direction::Int              
end

# Positional target (required) with validation
function StopAt(subject, var, target; direction::Int=0)
    var isa AbstractOrbitVar || error("var must be <: AbstractOrbitVar, got $(typeof(var))")
    StopAt(subject, var, target, direction)
end

"""
    make_calc(subject, var)

Create a calculation object from a subject and variable for use in stopping conditions.

# Arguments
- `subject`: The object to evaluate (e.g., `Spacecraft`, `Maneuver`, `CelestialBody`)
- `var`: The calculation variable (must be <: AbstractOrbitVar)

# Returns
A calculation object that can be used with `get_calc()` to evaluate the variable on the subject.

# Extensibility
This is an extensibility point for the stopping condition framework. Users and packages
should extend this function for new subject/variable combinations:

```julia
# Example extension for a custom subject type
make_calc(my_object::MyType, v::AbstractOrbitVar) = CustomCalc(my_object, v)
```

# Examples
```julia
# Built-in case: spacecraft orbital variables
calc = make_calc(spacecraft, PosX())
current_x = get_calc(calc)
```
"""
make_calc(subject, var) = error("make_calc not implemented for $(typeof(subject)), $(typeof(var))")

# Common case: orbit variables on spacecraft
make_calc(sc::Spacecraft, v::AbstractOrbitVar) = OrbitCalc(sc, v)

"""
    rebind(stop_condition::StopAt, owner_map::Dict{Any,Any})

Create a new StopAt condition with updated subject references based on an owner mapping.

# Arguments
- `stop_condition::StopAt`: The original stopping condition
- `owner_map::Dict{Any,Any}`: Mapping from old subjects to new subjects

# Returns
A new `StopAt` with the subject updated according to the mapping, while preserving
the variable, target, and direction unchanged.

# Usage
This function is used internally when transferring stopping conditions between
different propagation contexts where the subject objects may need to be remapped.

# Example
```julia
old_stop = StopAt(old_spacecraft, PosX(), 7000.0)
owner_map = Dict(old_spacecraft => new_spacecraft)
new_stop = rebind(old_stop, owner_map)
# new_stop.subject == new_spacecraft, other fields unchanged
```
"""
rebind(x::StopAt, owner_map::Dict{Any,Any}) =
    StopAt(get(owner_map, x.subject, x.subject), x.var, x.target, x.direction)

"""
   _subject_update_from_u!(subject, dynsys, u)

Update a subject from the integrator state u (specialize per subject type)
"""
_subject_update_from_u!(subject, dynsys, u) = error("No _subject_update_from_u! for $(typeof(subject))")

"""
   _subject_update_from_u!(sc::Spacecraft, dynsys, u)

Map Cartesian state slice to spacecraft struct state
"""
_subject_update_from_u!(sc::Spacecraft, dynsys, u) = begin
    pv = _posvel_from_u(u, dynsys, sc)
    set_posvel!(sc, pv)
    nothing
end

# Normalize to a Vector{Spacecraft}
_as_scvec(sc::Spacecraft) = Spacecraft[sc]
_as_scvec(v::Vector{<:Spacecraft}) = v

"""
   _sc_index(dynsys, sc::Spacecraft)

Find spacecraft index in a DynSys
"""
_sc_index(dynsys, sc::Spacecraft) = findfirst(x -> x === sc, getfield(dynsys, :spacecraft))

"""
   function _posvel_from_u(u, dynsys, sc::Spacecraft)
    
Extract a 6x1 Cartesian pos/vel slice for a spacecraft from integrator state u
"""
function _posvel_from_u(u, dynsys, sc::Spacecraft)
    idx = _sc_index(dynsys, sc)
    idx === nothing && error("StopAt: spacecraft not found in DynamicsSystem")
    i0 = 6*(idx-1) + 1
    return collect(@view u[i0:i0+5])
end

""" 
    _build_callback(cond::StopAt, dynsys)

Build a DifferentialEquations.jl ContinuousCallback for a StopAt condition.
"""
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

"""
    propagate(op::OrbitPropagator, sc_or_scs, stops...; direction=:forward, kwargs...)

Numerically integrate spacecraft equations of motion using propagator 
until stopping conditions are met.

# Arguments
- `op::OrbitPropagator`: Propagator configuration containing force model and integrator settings
- `sc_or_scs`: Single `Spacecraft` or `Vector{Spacecraft}` to propagate
- `stops...`: One or more `StopAt` stopping conditions

# Keyword Arguments
- `direction::Symbol=:forward`: Integration direction (`:forward` or `:backward`)
- `kwargs...`: Additional arguments passed to the underlying ODE solver

# Returns
`ODESolution` from DifferentialEquations.jl containing the complete trajectory solution.
Access final states via `sol.u[end]`, times via `sol.t`, or interpolate at any time.

# Examples
```julia
using Epicycle

# Spacecraft
sat = Spacecraft(
    state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time=Time("2015-09-21T12:23:12", TAI(), ISOT()),
    #name="SC-StopAt",
    coord_sys=CoordinateSystem(earth, ICRFAxes()),
)

# Forces + integrator
gravity = PointMassGravity(earth, (moon,sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Propagate to periapsis
propagate(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))

# Propagate backwards to node
propagate(prop, sat, StopAt(sat, PosX(), 0.0); direction=:backward)

# Propagate multiple spacecraft with multiple stopping conditions
sc1 = Spacecraft(); sc2 = Spacecraft() 
stop_sc1_node = StopAt(sc1, PosZ(), 0.0)
stop_sc2_periapsis = StopAt(sc2, PosDotVel(), 0.0; direction=+1)
propagate(prop, [sc1,sc2], stop_sc1_node, stop_sc2_periapsis)

```
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