# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

# ============================================================================
# Time-based stopping condition types
# ============================================================================

"""
    IntegratorTimeCalc <: AbstractCalcVariable

Base type for stopping conditions based on integrator time rather than spacecraft state.
These are specific to propagation contexts and cannot be used as general calculation variables.

# Subtypes
- `PropDurationSeconds`: Stop after a specified number of seconds
- `PropDurationDays`: Stop after a specified number of days

# Notes
These types are handled specially in `propagate!()` by setting the integration time span
directly rather than using callbacks. They derive from `AbstractCalcVariable` for type
safety in `StopAt` but do not have corresponding `make_calc()` implementations.
"""
abstract type IntegratorTimeCalc <: AbstractCalcVariable end

"""
    PropDurationSeconds <: IntegratorTimeCalc

Stopping condition based on elapsed time in seconds.

# Usage
```julia
# Stop after 1 hour (3600 seconds)
StopAt(sc, PropDurationSeconds(), 3600.0)
```
"""
struct PropDurationSeconds <: IntegratorTimeCalc end

"""
    PropDurationDays <: IntegratorTimeCalc

Stopping condition based on elapsed time in days.

# Usage
```julia
# Stop after 2.5 days
StopAt(sc, PropDurationDays(), 2.5)
```
"""
struct PropDurationDays <: IntegratorTimeCalc end

# ============================================================================
# Propagator and Stopping Condition Definitions
# ============================================================================

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
    StopAt{S,V<:AbstractCalcVariable,T}

Generic stopping condition for orbital propagation based on calculated quantities or time.

# Fields
- `subject::S`: The object to evaluate (e.g., `Spacecraft`, `Maneuver`, `CelestialBody`)
- `var::V`: The calculation variable to monitor (e.g., `PosX()`, `VelMag()`, `PropDurationSeconds()`)
- `target::T`: Target value to stop at (numeric value or vector matching calc output)
- `direction::Int`: Event crossing direction for state-based stops (-1: decreasing, 0: any, +1: increasing)
  * For time-based stops (PropDuration*), must be 0 (event crossing not applicable)
  * For state-based stops, controls which direction of zero-crossing triggers the event
- `detection::Symbol`: how the crossing is polled during integration.
  * `:discrete` (default): fast per-step polling on `integrator.u` with a `DiscreteCallback`;
    when a sign change is detected between two consecutive accepted steps, bisect on the
    Vern9 dense-output interpolant to locate the exact root, rewind the integrator to that
    root, and terminate. Cheapest option — ~5× less per-step overhead than `:continuous`.
  * `:continuous`: full `ContinuousCallback` — evaluates `g` on the interpolant at both
    step endpoints every accepted step, then root-finds when a sign change appears. Use
    only for signals that can double-cross within a single Vern9 step (fast oscillators
    with loose tolerances). Rare in astrodynamics.
- `rootfind_tol::Float64`: bisection convergence tolerance on the `g = calc - target`
  value at the root. Default `1e-9`. Applies to `:discrete`; `:continuous` inherits
  DiffEqBase's own root solver tolerance.

# Constructor
    StopAt(subject, var, target; direction::Int=0,
           detection::Symbol=:discrete, rootfind_tol::Real=1e-9)

# Notes
The `direction` field is for **event crossing direction** (state-based stops only).
For **time integration direction** (forward/backward), use the `direction` keyword in `propagate!()`.

Both `:discrete` and `:continuous` deliver the crossing to the same precision — they only
differ in how they poll for it. `sol.u[end]` and `sc.state` land at the interpolated root
under either mode; `sc.history`'s terminal entry is the root state.

# Examples
```julia
# Fast path (default): DiscreteCallback + interpolant bisection
stop_cond = StopAt(sc, PosZ(), 0.0; direction=+1)

# Escape hatch for a signal that oscillates faster than one Vern9 step
stop_fast = StopAt(sc, SomeFastOsc(), 0.0; detection=:continuous)

# Loosen the bisection tolerance if 1e-9 is tighter than needed
stop_loose = StopAt(sc, PosDotVel(), 0.0; direction=-1, rootfind_tol=1e-6)

# Time-based stops don't use root-finding at all; both fields are inert.
stop_time = StopAt(sc, PropDurationSeconds(), 3600.0)
```
"""
struct StopAt{S,V<:AbstractCalcVariable,T}
    subject::S
    var::V
    target::T
    direction::Int
    detection::Symbol
    rootfind_tol::Float64
end

# Positional target (required) with validation
function StopAt(subject, var, target;
                direction::Int = 0,
                detection::Symbol = :discrete,
                rootfind_tol::Real = 1e-9)
    var isa AbstractCalcVariable || error("var must be <: AbstractCalcVariable, got $(typeof(var))")

    detection in (:discrete, :continuous) ||
        throw(ArgumentError("detection must be :discrete or :continuous; got detection = $(repr(detection))"))
    rootfind_tol > 0 ||
        throw(ArgumentError("rootfind_tol must be > 0; got rootfind_tol = $rootfind_tol"))

    # Time-based stops don't use event crossing direction
    if var isa IntegratorTimeCalc && direction != 0
        error("Time-based stopping conditions must use direction=0 (event crossing direction not applicable for time-based stops)")
    end

    StopAt(subject, var, target, direction, detection, Float64(rootfind_tol))
end

# Convenience constructor for absolute time stopping
"""
    StopAt(subject::Spacecraft, target_time::Time; direction::Int=0)

Convenience constructor to stop at an absolute time by converting to elapsed seconds.
Supports both forward propagation (target after current) and backward propagation (target before current).

# Arguments
- `subject::Spacecraft`: The spacecraft being propagated
- `target_time::Time`: The absolute time to stop at (can be past or future)
- `direction::Int=0`: Event crossing direction (must be 0 for time-based stops)

# Notes
This `direction` parameter is for event crossing (not applicable to time stops, always use 0).
For **time integration direction** (forward/backward in time), use the `direction` keyword
in `propagate!()` with values `:forward`, `:backward`, or `:infer`.

# Example
```julia
using AstroEpochs

sat = Spacecraft(
    time=Time("2025-12-25T11:00:00", UTC(), ISOT()),
)
# Forward to future time (default direction=:forward in propagate! works)
stop_future = StopAt(sat, Time("2025-12-26T12:00:00", UTC(), ISOT()))
propagate!(prop, sat, stop_future)  # Uses default direction=:forward

# Backward to past time (use direction=:infer in propagate! to auto-detect)
stop_past = StopAt(sat, Time("2025-12-24T12:00:00", UTC(), ISOT()))
propagate!(prop, sat, stop_past; direction=:infer)  # Infers :backward from negative elapsed time
# OR explicitly:
propagate!(prop, sat, stop_past; direction=:backward)
```
"""
function StopAt(subject::Spacecraft, target_time::Time; direction::Int=0)
    # Use TT for Earth-centered, TDB for others (matches propagation)
    center_body = subject.coord_sys.origin
    
    # Convert both times to the appropriate dynamical time scale
    target_dyn = (center_body === earth) ? target_time.tt : target_time.tdb
    current_dyn = (center_body === earth) ? subject.time.tt : subject.time.tdb
    
    # Compute elapsed time in dynamical time seconds (can be negative for past times)
    elapsed_sec = (target_dyn.jd - current_dyn.jd) * 86400.0
    
    # No error for negative - supports backward propagation with direction=:infer
    return StopAt(subject, PropDurationSeconds(), elapsed_sec, direction, :discrete, 1e-9)
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

Build a callback for a `StopAt` condition. Dispatches on `cond.detection`:
`:discrete` (default) → `_build_hybrid_callback`; `:continuous` → `_build_continuous_callback`.
"""
function _build_callback(cond::StopAt, dynsys)
    return cond.detection === :discrete ? _build_hybrid_callback(cond, dynsys) :
                                          _build_continuous_callback(cond, dynsys)
end

"""
    _build_hybrid_callback(cond::StopAt, dynsys)

Fast path (default). One `g(u,t)` evaluation per accepted step on `integrator.u` (no dense
interpolant polling). When a sign change appears between two consecutive step endpoints,
bisect on the Vern9 dense output to locate the exact root, rewind the integrator to that
root, and terminate. `sol.u[end]`/`sol.t[end]` land at the root state and time; the existing
`_update_structs!` path then delivers the root to `sc.state` and `sc.history` unchanged.
"""
function _build_hybrid_callback(cond::StopAt, dynsys)
    subject = cond.subject
    var     = cond.var
    target  = cond.target
    dir     = cond.direction
    tol     = cond.rootfind_tol
    calc    = make_calc(subject, var)

    g_prev = Ref(NaN)     # (calc - target) at the previous accepted-step endpoint
    t_prev = Ref(NaN)     # elapsed time at that endpoint

    # g(u) = calc(u) - target. Mutates subject as a side effect via _subject_update_from_u!
    # (unavoidable given the get_calc(calc) contract; only subject.state is touched).
    function g_at(u)
        _subject_update_from_u!(subject, dynsys, u)
        return get_calc(calc) - target
    end

    function cond_fn(u, t, _integ)
        g_now = g_at(u)
        if isnan(g_prev[])
            g_prev[] = g_now; t_prev[] = t
            return false
        end
        crossed = dir < 0 ? (g_prev[] > 0 && g_now ≤ 0) :
                  dir > 0 ? (g_prev[] < 0 && g_now ≥ 0) :
                            (sign(g_prev[]) != sign(g_now))
        if !crossed
            g_prev[] = g_now; t_prev[] = t
        end
        return crossed
    end

    function affect!(integ)
        tl, th = t_prev[], integ.t
        gl     = g_prev[]                       # opposite sign of g_at(u_th) at entry
        t_mid  = 0.5 * (tl + th)
        u_mid  = integ(t_mid)                   # Vern9 dense output
        for _ in 1:60
            gm = g_at(u_mid)
            (abs(gm) < tol || (th - tl) < 1e-6) && break
            if (gl > 0) == (gm > 0)
                tl = t_mid; gl = gm
            else
                th = t_mid
            end
            t_mid = 0.5 * (tl + th)
            u_mid = integ(t_mid)
        end
        # Rewind integrator to the root. save_everystep already appended (u_overshoot, t_overshoot)
        # for this step, so overwrite that last sol entry in place — this keeps sol.t monotone
        # and avoids the alternatives of pop! (breaks solver internal indexing) or leaving the
        # overshoot in place (violates the contract that sol/history end at the root).
        integ.u .= u_mid
        integ.t  = t_mid
        if length(integ.sol.t) > 0
            integ.sol.t[end] = t_mid
            integ.sol.u[end] = collect(u_mid)
        end
        terminate!(integ)
    end

    return DiscreteCallback(cond_fn, affect!; save_positions=(false, false))
end

"""
    _build_continuous_callback(cond::StopAt, dynsys)

Escape hatch for signals that can double-cross within a single Vern9 step (fast oscillators
with loose tolerances). Full `ContinuousCallback` with interp_points=2 and save_positions
disabled — same root-find precision as the hybrid path, but pays for dense-output
evaluation at every accepted step whether a sign change exists or not.
"""
function _build_continuous_callback(cond::StopAt, dynsys)
    subject = cond.subject
    var     = cond.var
    target  = cond.target
    dir     = cond.direction

    calc = make_calc(subject, var)

    function g(u, t, _integ)
        _subject_update_from_u!(subject, dynsys, u)
        val = get_calc(calc)
        return val - target
    end
    term!(integ) = terminate!(integ)

    # interp_points=2: only evaluate g at the two accepted-step endpoints via the
    # interpolant. save_positions=(false,false): don't re-save the ODE state per eval.
    if dir == 0
        return ContinuousCallback(g, term!;
            interp_points=2, save_positions=(false, false), rootfind=true)
    elseif dir > 0
        return ContinuousCallback(g, term!; affect_neg! = (_integ)->nothing,
            interp_points=2, save_positions=(false, false), rootfind=true)
    else
        return ContinuousCallback(g, (_integ)->nothing; affect_neg! = term!,
            interp_points=2, save_positions=(false, false), rootfind=true)
    end
end

"""
    propagate!(op::OrbitPropagator, sc_or_scs, stops...; direction=:forward, kwargs...)

Numerically integrate spacecraft equations of motion using propagator 
until stopping conditions are met.

# Arguments
- `op::OrbitPropagator`: Propagator configuration containing force model and integrator settings
- `sc_or_scs`: Single `Spacecraft` or `Vector{Spacecraft}` to propagate
- `stops...`: One or more `StopAt` stopping conditions

# Keyword Arguments
- `direction::Symbol=:forward`: Time integration direction
  * `:forward` - Integrate forward in time (default)
  * `:backward` - Integrate backward in time
  * `:infer` - Automatically determine from time-based stop conditions (duration sign or time comparison)
- `kwargs...`: Additional arguments passed to the underlying ODE solver

# Notes
The `direction` keyword controls **time integration direction** (which way time moves).
This is different from `StopAt`'s `direction` field, which controls **event crossing direction**
for state-based stops (increasing/decreasing zero-crossing detection).

# Returns
`ODESolution` from DifferentialEquations.jl containing the complete trajectory solution.
Access final states via `sol.u[end]`, times via `sol.t`, or interpolate at any time.

# Examples
```julia
using AstroEpochs, AstroStates, AstroFrames, AstroUniverse 
using AstroModels, AstroCallbacks, AstroProp, OrdinaryDiffEq

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
propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=+1))

# Propagate backwards to node
propagate!(prop, sat, StopAt(sat, PosX(), 0.0); direction=:backward)

# Propagate multiple spacecraft with multiple stopping conditions
sc1 = Spacecraft(); sc2 = Spacecraft() 
stop_sc1_node = StopAt(sc1, PosZ(), 0.0)
stop_sc2_periapsis = StopAt(sc2, PosDotVel(), 0.0; direction=+1)
propagate!(prop, [sc1,sc2], stop_sc1_node, stop_sc2_periapsis)

```
"""
function propagate!(op::OrbitPropagator, sc_or_scs, stops...;
                   direction::Symbol = :forward, kwargs...)
    scv = _as_scvec(sc_or_scs)
    dyn = DynSys(spacecraft=scv, forces=op.forces)

    # Separate time-based from state-based stopping conditions
    # Time-based conditions will be handled by setting tspan in the main propagate!()
    time_conds = filter(_is_time_condition, stops)
    state_conds = filter(!_is_time_condition, stops)

    # Build callbacks only from state-based conditions
    callbacks = map(s -> _build_callback(s, dyn), collect(state_conds))
    cbset = isempty(callbacks) ? nothing :
            length(callbacks) == 1 ? callbacks[1] : CallbackSet(callbacks...)

    # Delegate to existing DynSys-based propagate! with both callbacks and time conditions
    # Note: We pass time_conds separately, not as callbacks
    return propagate!(dyn, op.integ, cbset, time_conds...; direction=direction, kwargs...)
end

# Helper: detect time-based stopping conditions
_is_time_condition(::StopAt{<:Any, <:IntegratorTimeCalc}) = true
_is_time_condition(::Any) = false