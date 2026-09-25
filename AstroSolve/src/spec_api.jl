# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The do-block layer: readers, partial, and the verbs both phase kinds share.
#

# ─────────────────────────────────────────────────────────────────────────────
# Internal registry — per-phase typed-API metadata
# Keyed by objectid(phase); stores state/control type constructors and
# maps registered constraint/objective handles to their internal objects.
# ─────────────────────────────────────────────────────────────────────────────

struct _PhaseRegistry
    state_type    ::Any    # e.g. PlanarState (the UnionAll, not an instance)
    control_type  ::Any    # e.g. PlanarControl
    target        ::Any    # opaque instance — framework calls set!(target, y_s); nothing if unused
    n_states      ::Int
    n_controls    ::Int
end

const _phase_registry = Dict{UInt64, _PhaseRegistry}()

"""
    n_components(::Type{T}) -> Int

How many scalars a state or control type holds, so a phase can size itself
from the type rather than being told twice. Defaults to the field count;
override it for a type whose fields are not all scalars.
"""
n_components(::Type{T}) where {T} = fieldcount(T)
n_components(T::UnionAll)         = fieldcount(Base.unwrap_unionall(T))

# ─────────────────────────────────────────────────────────────────────────────
# The readers a spec names — a derivative is a fact about a function
#
# Keyed by the function, not by the phase or a handle, because the derivative
# of `speed` does not change with where `speed` is used. Declare it once and
# every constraint, phase and objective that names the function gets it. Same
# reasoning that makes `setter` a trait on a quantity.
#
# A partial takes exactly the arguments its function takes. Where the framework
# hands over fewer (dynamics Jacobians get no `model`), the wrapper supplies
# them, so there is one rule to remember rather than two.
# ─────────────────────────────────────────────────────────────────────────────

"""
    PropulsionModel(; mu, Isp = 0.0, Tmax = 0.0, g0 = 9.80665)

The physical constants a phase flies under.

Collocation already puts these behind `model = …` — a user struct its dynamics
reads. MGAnDSMs and Sims-Flanagan took them as constructor keywords instead, so
the same fact lived in two places depending on the transcription. This is the
convenience type for the common case; any struct with the right property names
works, exactly as it does for collocation.

# Fields
- `mu`: the gravitational parameter of the central body, which sets the length and time units every
  other quantity in the phase follows.
- `Isp`: specific impulse, in the phase's time units.
- `Tmax`: maximum thrust per thruster, in units consistent with `mu` and the masses.
- `g0`: the reference gravity that turns `Isp` into an exhaust velocity, and which has to match the
  units the phase flies in, so a phase in kilometres and seconds takes `9.80665e-3`.

# Example
```julia
PropulsionModel(mu = 1.32712440018e11, Isp = 3000.0, Tmax = 1.0e-3, g0 = 9.80665e-3)
```
"""
struct PropulsionModel{T<:Real}
    mu   ::T
    Isp  ::T
    Tmax ::T
    g0   ::T
end
PropulsionModel(; mu, Isp = 0.0, Tmax = 0.0, g0 = 9.80665) =
    PropulsionModel(promote(float(mu), float(Isp), float(Tmax), float(g0))...)

"""Read a constant from a model object, falling back to an explicit keyword."""
_model_field(model, name::Symbol, explicit, default) =
    explicit !== nothing ? explicit :
    (model !== nothing && hasproperty(model, name)) ? getproperty(model, name) : default

# What a spec can name. These are the second argument to `@partial` and the
# first to `Vary`, so one name serves both: `Vary(state, phase)` and
# `@partial(f, state)`. The tag structs above are internal from here on.
#
# They are readers, not markers. `state(phase)` returns the phase's state, the
# way `delta_v(maneuver)` returns a delta-V — declared here with no methods so
# each transcription supplies its own, and so asking a phase for a quantity it
# does not have is an error rather than a silent echo of the argument.

# `state` is AstroCallbacks'. A spacecraft has one and so does a phase, and it
# means the same thing about both, so the phase methods extend that generic
# rather than declaring a rival of the same name.
"""
    control(phase_or_context)

The control block of a phase, which is what the transcription steers the dynamics with.

Varied with `Vary(control, phase; guess = ..., lower_bound = ..., upper_bound = ...)` and
differentiated against with `@partial(f, control)`.

# Returns
On a phase, the control history as a matrix with one row per control component and one column per
node. Inside a path function, the control at that node. Units are the caller's.

# Example
```julia
control(phase)
```
"""
function control end

"""
    parameter(phase_or_context)

A phase's static parameters, which are constant over the arc rather than functions of time.

Varied with `Vary(parameter, phase; guess = ...)`, which is how a constant in the equations of
motion becomes something the solver finds.

# Returns
The parameter vector. Units are the caller's, and are whatever the dynamics read them as.

# Example
```julia
parameter(phase)
```
"""
function parameter end

"""
    segment_durations(phase)

The durations of a shooting phase's segments, which together span the phase.

Varied with `Vary(segment_durations, phase; guess = ...)` when the segment lengths are themselves
to be found rather than split evenly.

# Returns
One duration per segment, in the time units of the dynamics.

# Example
```julia
segment_durations(phase)
```
"""
function segment_durations end

"""
    initial_time(phase_or_context)

The time at the start of a phase.

Varied with `Vary(initial_time, phase; guess = ..., lower_bound = ..., upper_bound = ...)`, which
is how a departure date becomes an optimization variable.

# Returns
The start time as a `Float64`, in the time units of the dynamics. Before `solve!` it is the start of
`tspan` or the guess; after it, the solution.

# Example
```julia
initial_time(phase)
```
"""
function initial_time end

"""
    final_time(phase_or_context)

The time at the end of a phase.

Varied with `Vary(final_time, phase; guess = ..., lower_bound = ..., upper_bound = ...)`, which is
how a minimum-time problem lets the arrival float. Minimizing it is the objective of a
minimum-time problem.

# Returns
The end time as a `Float64`, in the time units of the dynamics. Before `solve!` it is the end of
`tspan` or the guess; after it, the solution.

# Example
```julia
final_time(phase)
```
"""
function final_time end

"""
    initial_state(phase_or_context)

The state at a phase's left end, which is where its boundary conditions at departure are written.

Constrained with `Constraint(f, phase; at = Initial())`, where the function reads this quantity.

# Returns
The state at the left end, as the phase's state type when it declares one and as a plain vector
otherwise. Units are the caller's.

# Example
```julia
initial_state(phase)
```
"""
function initial_state end

"""
    final_state(phase_or_context)

The state at a phase's right end, which is where its boundary conditions at arrival are written.

Constrained with `Constraint(f, phase; at = Final())`, where the function reads this quantity.

# Returns
The state at the right end, as the phase's state type when it declares one and as a plain vector
otherwise. Units are the caller's.

# Example
```julia
final_state(phase)
```
"""
function final_state end

"""
    departure_vinf(phase_or_context)

The excess velocity with which a phase leaves the body at its left control point, as a vector.

Varied with `Vary(departure_vinf, phase; guess = ...)`. Its magnitude is what a launch vehicle's
capability bounds, and what a flyby matches across a linkage.

# Returns
A three-component excess velocity, in the velocity units of the dynamics.

# Example
```julia
departure_vinf(phase)
```
"""
function departure_vinf end

"""
    arrival_vinf(phase_or_context)

The excess velocity with which a phase arrives at the body at its right control point, as a vector.

Varied with `Vary(arrival_vinf, phase; guess = ...)`. Its magnitude is what an arrival capture
budget bounds, and what a flyby matches across a linkage.

# Returns
A three-component excess velocity, in the velocity units of the dynamics.

# Example
```julia
arrival_vinf(phase)
```
"""
function arrival_vinf end

"""
    initial_mass(phase_or_context)

The spacecraft mass at the start of a phase.

Varied with `Vary(initial_mass, phase; lower_bound = ..., upper_bound = ...)`, and pinned by giving
the two bounds the same value when the launch mass is known.

# Returns
The mass as a `Float64`, in the mass units of the dynamics.

# Example
```julia
initial_mass(phase)
```
"""
function initial_mass end

"""
    final_mass(phase_or_context)

The spacecraft mass at the end of a phase, which is the mass delivered.

Varied with `Vary(final_mass, phase; guess = ..., lower_bound = ..., upper_bound = ...)`. A phase
propagates mass from both ends to its match point, so this is a variable the match point constrains
rather than a value read off the forward propagation.

# Returns
The mass as a `Float64`, in the mass units of the dynamics.

# Example
```julia
final_mass(phase)
```
"""
function final_mass end

"""
    deep_space_dv(phase_or_context)

The impulses a gravity-assist phase applies between its Keplerian arcs.

Varied with `Vary(deep_space_dv, phase; guess = ...)`. There is one impulse per maneuver, and a
phase with `n` maneuvers has `n + 1` arcs.

# Returns
A matrix with three rows and one column per maneuver, in the velocity units of the dynamics.

# Example
```julia
deep_space_dv(phase)
```
"""
function deep_space_dv end

"""
    arc_fractions(phase_or_context)

The fractions of a gravity-assist phase's duration taken by each of its Keplerian arcs.

Varied with `Vary(arc_fractions, phase; guess = ...)`. The phase carries a built-in equality
holding the fractions to a sum of one, so they are a division of the duration rather than
independent lengths.

# Returns
One fraction per arc, dimensionless, summing to one at a feasible point.

# Example
```julia
arc_fractions(phase)
```
"""
function arc_fractions end

"""
    forward_control(phase_or_context)

The throttle vectors of the segments a Sims-Flanagan phase propagates forward from its left control
point.

Varied with `Vary(forward_control, phase; guess = ..., lower_bound = ..., upper_bound = ...)` and
differentiated against with `@partial(f, forward_control)`. A throttle is normalized, so the thrust
it commands is its magnitude times the phase's `Tmax`, and holding that magnitude inside the unit
ball is a `Constraint` the caller writes.

# Returns
A matrix with three rows and one column per forward segment, dimensionless.

# Example
```julia
forward_control(phase)
```
"""
function forward_control end

"""
    backward_control(phase_or_context)

The throttle vectors of the segments a Sims-Flanagan phase propagates backward from its right
control point.

Varied with `Vary(backward_control, phase; guess = ..., lower_bound = ..., upper_bound = ...)` and
differentiated against with `@partial(f, backward_control)`. Backward segments are stored in
forward time order, so the first column is the one nearest the match point.

# Returns
A matrix with three rows and one column per backward segment, dimensionless.

# Example
```julia
backward_control(phase)
```
"""
function backward_control end

# Collocation reads its own cached blocks.
state(p::CollocationPhase)        = p._Y
control(p::CollocationPhase)      = p._U
initial_time(p::CollocationPhase) = p._t0
final_time(p::CollocationPhase)   = p._tf

"""Map a quantity, plus where it is evaluated, onto the framework's own tag."""
_path_tag(q) = q === state     ? State()     :
               q === control   ? Control()   :
               q === parameter ? Parameter() : nothing

_bnd_tag(q, initial::Bool) =
    q === state        ? (initial ? InitialState() : FinalState()) :
    q === initial_time ? InitialTime() :
    q === final_time   ? FinalTime()   : nothing

_pname(f) = string(nameof(typeof(f)))

# ─────────────────────────────────────────────────────────────────────────────
# @partial — a derivative is a fact about a function
#
# Declared against the function rather than against a phase or a handle, because
# the derivative of `speed` does not change with where `speed` is used. uc5
# constrains `speed` twice and declares its partial once. Same reasoning that
# makes a label a trait on a quantity.
#
# A declaration is a method on `partial`. Nothing is registered, so re-running a
# script redefines it the way redefining any function does, and nothing outlives
# the session that the session did not put there. It replaced `partial!`, which
# wrote into a dictionary keyed by `typeof(f)` and made the ordinary REPL loop
# fail: editing a bound and running the file again raised on a partial the user
# had not touched.
# ─────────────────────────────────────────────────────────────────────────────

"""
    partial(f, wrt, args...)

The derivative of `f` with respect to `wrt`, evaluated at `args`, which are the
arguments `f` itself takes. Declared with [`@partial`](@ref).
"""
function partial end

"""
    @partial(f, wrt) do <same arguments as f> ... end

Declare the derivative of `f` with respect to `wrt`.

`wrt` names a variable the same way `Vary` does: `state`, `control`,
`initial_time`, `final_time`. Returns ∂f/∂wrt with a row per component of `f`
and a column per component of the variable, except for a dynamics partial,
which writes into the matrix it is given.

A partial not declared is taken by automatic differentiation, and
`check_partials` reports which is which.

# Examples
```julia
speed(c) = [state(c).v]

@partial(speed, state) do c
    [0.0 0.0 1.0]
end
```
"""
macro partial(body, f, wrt)
    body isa Expr && body.head === :-> || throw(ArgumentError(
        "usage: @partial(f, wrt) do <the arguments f takes> ... end"))
    ps   = body.args[1]
    args = ps isa Expr && ps.head === :tuple ? ps.args : [ps]
    # AstroSolve.partial by name, not `partial`. Escaped, the name resolves in
    # the caller's scope, so a declaration written at the top level of a script
    # defines a new function there and the framework never sees it — the
    # declaration is accepted and the derivative silently falls back to AD.
    return quote
        function $(GlobalRef(@__MODULE__, :partial))(
                ::typeof($(esc(f))), ::typeof($(esc(wrt))),
                $(map(esc, args)...))
            $(esc(body.args[2]))
        end
    end
end

"The `wrt` readers that `f` has a `partial` method for."
function _declared_wrts(f)
    out = Any[]
    for m in methods(partial)
        p = m.sig.parameters
        length(p) >= 3 && p[2] === typeof(f) || continue
        isdefined(p[3], :instance) || continue
        push!(out, p[3].instance)
    end
    return out
end

"""Every partial declared for `f`, as (reader, function) pairs."""
partials_of(f) = Tuple{Any,Function}[(w, (a...) -> partial(f, w, a...))
                                     for w in _declared_wrts(f)]

"""
    partial(f, wrt, args...)

The partial of `f` with respect to `wrt`, evaluated at `args`, which are the arguments `f` itself
takes. Declared with [`@partial`](@ref).

# Returns
The derivative, shaped as the declaration returns it, or `nothing` when no partial was declared for
that pair, which is how the framework knows to differentiate `f` itself instead.

# Example
<!-- doc-fragment -->
```julia
@partial(thrust_ball, state) do c
    [2.0 * state(c).x]
end

partial(thrust_ball, state, context)
```
"""
function partial(f, wrt)
    for (tag, fn) in partials_of(f)
        typeof(tag) === typeof(wrt) && return fn
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# _prepare_model! — single place to build typed structs and sync the target.
#
# Two-arg form  (dynamics, path constraints, Lagrange per-node):
#   y_s, u_s = _prepare_model!(reg, Ty, y_vec, u_vec)
#   - Builds reg.state_type{Ty}(y_vec...)  and reg.control_type{Ty}(u_vec...)
#   - Calls set!(reg.target, y_s)          (always, state is required)
#   - Calls set!(reg.target, u_s)          (only when the user has defined a
#                                            matching set! overload for the control)
#
# One-arg form  (boundary constraints, Mayer — state-only context):
#   y_s = _prepare_model!(reg, Ty, y_vec)
#   - Builds and syncs state only.
#   - For boundary contexts call once for y0, once for yf; target ends at yf.
#
# If reg.state_type / reg.control_type is nothing, the raw vector is returned
# unchanged and no set! is called.
# ─────────────────────────────────────────────────────────────────────────────

function _prepare_model!(reg::_PhaseRegistry, ::Type{Ty}, y_vec, u_vec) where {Ty}
    # If the value is already a typed struct (e.g. ctx.y0 in a BoundaryContext that
    # carries pre-built structs), skip construction to avoid splat-on-struct errors.
    y_s = (reg.state_type !== nothing && !(y_vec isa reg.state_type)) ?
              reg.state_type{Ty}(y_vec...) : y_vec
    u_s = (reg.control_type !== nothing && !(u_vec isa reg.control_type)) ?
              reg.control_type{Ty}(u_vec...) : u_vec
    if reg.target !== nothing
        set!(reg.target, y_s)
        applicable(set!, reg.target, u_s) && set!(reg.target, u_s)
    end
    return y_s, u_s
end

function _prepare_model!(reg::_PhaseRegistry, ::Type{Ty}, y_vec) where {Ty}
    y_s = (reg.state_type !== nothing && !(y_vec isa reg.state_type)) ?
              reg.state_type{Ty}(y_vec...) : y_vec
    if reg.target !== nothing
        set!(reg.target, y_s)
    end
    return y_s
end

# Handle types returned to the user — thin wrappers around internal objects
struct BoundaryConstraintHandle
    bc   ::BoundaryConstraint   # internal object
    phase::CollocationPhase
end

struct MayerHandle
    phase::CollocationPhase
end

struct LagrangeHandle
    phase   ::CollocationPhase
    fn      ::Function             # integrand (y, u, p, t) -> scalar
end

struct ObjectiveHandle
    phase   ::CollocationPhase
    sense   ::Union{Min,Max}
    mayer   ::Union{Nothing, MayerHandle}
    lagrange::Union{Nothing, LagrangeHandle}
    # mutable fields via Ref so we can attach terms after creation
    _mayer_ref   ::Ref{Union{Nothing,MayerHandle}}
    _lagrange_ref::Ref{Union{Nothing,LagrangeHandle}}
end

ObjectiveHandle(phase, sense) =
    ObjectiveHandle(phase, sense, nothing, nothing,
                    Ref{Union{Nothing,MayerHandle}}(nothing),
                    Ref{Union{Nothing,LagrangeHandle}}(nothing))

struct PathConstraintHandle
    pc   ::Any                    # internal path constraint record
    phase::CollocationPhase
end

# ─────────────────────────────────────────────────────────────────────────────
# Sequence constructor sugar:  Sequence(phase)
# ─────────────────────────────────────────────────────────────────────────────

function Sequence(phase::CollocationPhase)
    seq = Sequence()
    add_sequence!(seq, phase)
    return seq
end

# ─────────────────────────────────────────────────────────────────────────────
# solve! — alias for solve_trajectory!
# ─────────────────────────────────────────────────────────────────────────────

# solve!(seq) lives in the spec layer now, taking `method`. Leaving an alias
# here would be the same shadowing pattern that cost us AstroEpochs.Time and
# SolverVariable: two definitions, one silently winning.

# ─────────────────────────────────────────────────────────────────────────────
# set_state! — new API
#   set_state!(phase, StateType; guess, lower_bounds, upper_bounds, scale, shift)
#   StateType is the UnionAll (e.g. BrachState), not an instance.
#   Stores the type in the registry; builds a DirectSolverVariable internally.
# ─────────────────────────────────────────────────────────────────────────────

# Sentinel array types for the new API (not user-visible)
# S/C are the user's concrete state/control UnionAll types so that
# eltype_state / eltype_control return them correctly.
struct _NewStateArray{S}   <: AbstractStateArray{S};   data::Matrix{Float64} end
struct _NewControlArray{C} <: AbstractControlArray{C}; data::Matrix{Float64} end

function set_state!(phase::CollocationPhase, StateType::Type;
                    guess        = nothing,
                    lower_bounds,
                    upper_bounds,
                    scale        = nothing,
                    shift        = nothing,
                    target       = nothing)
    ns  = length(lower_bounds)
    sc  = scale === nothing ? ones(ns)  : Float64.(scale)
    sh  = shift === nothing ? zeros(ns) : Float64.(shift)
    g   = guess === nothing ? zeros(ns, 2) : Float64.(guess)

    v = DirectSolverVariable(
        _NewStateArray{StateType}(g),
        Float64.(lower_bounds), Float64.(upper_bounds),
        sc, sh, "state",
        vec(g[:, 1]),   # initial physical value = first column of guess
    )
    set_state!(phase, v)

    # Store the type constructor in the registry
    reg = get(_phase_registry, objectid(phase), nothing)
    nc  = reg === nothing ? 0 : reg.n_controls
    ct  = reg === nothing ? nothing : reg.control_type
    _phase_registry[objectid(phase)] = _PhaseRegistry(StateType, ct, target, ns, nc)
    nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# set_control! — new API
# ─────────────────────────────────────────────────────────────────────────────

function set_control!(phase::CollocationPhase, ControlType::Type;
                      guess        = nothing,
                      lower_bounds,
                      upper_bounds,
                      scale        = nothing,
                      shift        = nothing)
    nc  = length(lower_bounds)
    sc  = scale === nothing ? ones(nc)  : Float64.(scale)
    sh  = shift === nothing ? zeros(nc) : Float64.(shift)
    g   = guess === nothing ? zeros(nc, 2) : Float64.(guess)

    v = DirectSolverVariable(
        _NewControlArray{ControlType}(g),
        Float64.(lower_bounds), Float64.(upper_bounds),
        sc, sh, "control",
        vec(g[:, 1]),
    )
    set_control!(phase, v)

    reg = get(_phase_registry, objectid(phase), nothing)
    ns  = reg === nothing ? 0 : reg.n_states
    st  = reg === nothing ? nothing : reg.state_type
    tt  = reg === nothing ? nothing : reg.target
    _phase_registry[objectid(phase)] = _PhaseRegistry(st, ControlType, tt, ns, nc)
    nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# set_parameter! — new API (no user type required)
#   set_parameter!(phase; guess, lower_bounds, upper_bounds, scale, shift)
#   Parameters are accessed inside dynamics/Jacobians via get_param(ctx).
# ─────────────────────────────────────────────────────────────────────────────

struct _NewParam <: AbstractParameter; data::Vector{Float64} end

function set_parameter!(phase::CollocationPhase;
                        guess        = nothing,
                        lower_bounds,
                        upper_bounds,
                        scale        = nothing,
                        shift        = nothing)
    np = length(lower_bounds)
    sc = scale === nothing ? ones(np)  : Float64.(scale)
    sh = shift === nothing ? zeros(np) : Float64.(shift)
    g  = guess === nothing ? copy(sh)  : Float64.(guess)
    v  = DirectSolverVariable(
        _NewParam(g),
        Float64.(lower_bounds), Float64.(upper_bounds),
        sc, sh, "param", g,
    )
    set_parameter!(phase, v)
    nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# set_initial_time! / set_final_time! — new API
#   Fixed:  set_initial_time!(phase; equality = value)
#   Free:   set_final_time!(phase, guess; lower_bounds, upper_bounds, scale, shift)
# ─────────────────────────────────────────────────────────────────────────────

# Sentinel time type for the new API
struct _NewTime <: AbstractTime; value::Float64 end

function set_initial_time!(phase::CollocationPhase; equality::Float64)
    v = DirectSolverVariable(
        _NewTime(equality),
        [equality], [equality],
        [1.0], [0.0], "t0", [equality],
    )
    set_initial_time!(phase, v)
    nothing
end

function set_initial_time!(phase::CollocationPhase, guess::Float64;
                            lower_bounds::Float64,
                            upper_bounds::Float64,
                            scale ::Float64 = 1.0,
                            shift ::Float64 = 0.0)
    v = DirectSolverVariable(
        _NewTime(guess),
        [lower_bounds], [upper_bounds],
        [scale], [shift], "t0", [guess],
    )
    set_initial_time!(phase, v)
    nothing
end

function set_final_time!(phase::CollocationPhase; equality::Float64)
    v = DirectSolverVariable(
        _NewTime(equality),
        [equality], [equality],
        [1.0], [0.0], "tf", [equality],
    )
    set_final_time!(phase, v)
    nothing
end

function set_final_time!(phase::CollocationPhase, guess::Float64;
                         lower_bounds::Float64,
                         upper_bounds::Float64,
                         scale ::Float64 = 1.0,
                         shift ::Float64 = 0.0)
    v = DirectSolverVariable(
        _NewTime(guess),
        [lower_bounds], [upper_bounds],
        [scale], [shift], "tf", [guess],
    )
    set_final_time!(phase, v)
    nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# set_dynamics_jacobian! — new API, tag dispatch
# ─────────────────────────────────────────────────────────────────────────────

function set_dynamics_jacobian!(phase::CollocationPhase, ::State, fn::Function)
    ph  = phase
    var = phase.state_var
    # A partial takes the same arguments as the function it differentiates.
    # The dynamics is called f(dy, y, u, params, t, model) whether or not the
    # phase has a control, with nothing in the u slot when it does not, so a
    # dynamics partial is called the same way. This used to drop the u slot
    # on a control-free phase and hand over the whole context in place of the
    # parameters.
    wrapped = function _state_jac_wrapper!(dF, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            fn(dF, y_s, reg.control_type === nothing ? nothing : u_s, ctx.params, t)
        else
            fn(dF, y_vec, ctx, t)
        end
    end
    set_dynamics_jacobian!(phase, var, wrapped)
end

function set_dynamics_jacobian!(phase::CollocationPhase, ::Control, fn::Function)
    ph  = phase
    var = phase.control_var
    wrapped = function _control_jac_wrapper!(dF, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            fn(dF, y_s, u_s, ctx.params, t)
        else
            fn(dF, y_vec, ctx, t)
        end
    end
    set_dynamics_jacobian!(phase, var, wrapped)
end

function set_dynamics_jacobian!(phase::CollocationPhase, ::Parameter, fn::Function)
    ph  = phase
    var = phase.param_var
    var === nothing && throw(ArgumentError(
        "a parameter must be declared before its dynamics Jacobian can be registered; " *
        "call Vary(parameter, phase; ...) first"))
    wrapped = function _param_jac_wrapper!(dF, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            fn(dF, y_s, reg.control_type === nothing ? nothing : u_s, ctx.params, t)
        else
            fn(dF, y_vec, ctx, t)
        end
    end
    set_dynamics_jacobian!(phase, var, wrapped)
end

# ─────────────────────────────────────────────────────────────────────────────
# _typed_struct — instantiate a user struct from a raw vector using the
#   registered type constructor.  Calls StateType{Float64}(vals...) or
#   ControlType{Float64}(vals...).
# ─────────────────────────────────────────────────────────────────────────────

function _typed_state(phase::CollocationPhase, v::Vector{Float64})
    reg = _phase_registry[objectid(phase)]
    reg.state_type{Float64}(v...)
end

function _typed_control(phase::CollocationPhase, v::Vector{Float64})
    reg = _phase_registry[objectid(phase)]
    reg.control_type{Float64}(v...)
end

# ─────────────────────────────────────────────────────────────────────────────
# add_boundary_constraint! — new API
#   fn(y0, yf, p, t0, tf) -> Vector
#   Returns a BoundaryConstraintHandle for add_jacobian!
# ─────────────────────────────────────────────────────────────────────────────

function add_boundary_constraint!(fn::Function, phase::CollocationPhase;
                                  name::Symbol,
                                  equality    = nothing,
                                  lower_bounds = nothing,
                                  upper_bounds = nothing)
    # Resolve bounds
    if equality !== nothing
        lb = Float64.(equality)
        ub = Float64.(equality)
    else
        lb = Float64.(lower_bounds)
        ub = Float64.(upper_bounds)
    end
    n = length(lb)

    # Wrap user fn(y0, yf, p, t0, tf, model) into the internal ctx-based BoundaryFunction form.
    # IMPORTANT: use ctx.y0/ctx.yf/ctx.t0/ctx.tf — not phase fields —
    # so ForwardDiff can perturb through them in the AD Jacobian path.
    ph = phase
    wrapped = BoundaryFunction(phase; name = string(name)) do ctx
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            T    = eltype(ctx.y0)
            y0_s = _prepare_model!(reg, T, ctx.y0)   # target synced to y0
            yf_s = _prepare_model!(reg, T, ctx.yf)   # target synced to yf (overwrites; terminal is most relevant)
            fn(y0_s, yf_s, nothing, ctx.t0, ctx.tf, ph.model)
        else
            fn(ctx.y0, ctx.yf, nothing, ctx.t0, ctx.tf, ph.model)
        end
    end

    bc = BoundaryConstraint(wrapped, lb, ub)
    add_constraint!(phase, bc)
    return BoundaryConstraintHandle(bc, phase)
end

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! — new API, tag dispatch for BoundaryConstraintHandle
#   Jacobian closure receives only the relevant typed arg (tag selects it).
#   Framework receives a zero-arg closure returning the Jacobian matrix.
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::BoundaryConstraintHandle, ::InitialState, fn::Function)
    # fn(y0, yf, p, t0, tf, model) -> (n_bc × ns) compact matrix.
    # Framework expects full (n_bc × ns*N) block; embed at columns 1:ns.
    ph = h.phase
    zero_arg = () -> begin
        reg     = get(_phase_registry, objectid(ph), nothing)
        y0_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._y0) : ph._y0
        yf_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._yf) : ph._yf
        compact = fn(y0_s, yf_s, nothing, ph._t0, ph._tf, ph.model)
        ns = ph._n_states
        N  = n_control_nodes(ph)
        J  = zeros(size(compact, 1), ns * N)
        J[:, 1:ns] .= compact
        J
    end
    add_jacobian!(h.bc.calc, ph.state_var, zero_arg)
end

function add_jacobian!(h::BoundaryConstraintHandle, ::FinalState, fn::Function)
    # fn(y0, yf, p, t0, tf, model) -> (n_bc × ns) compact matrix; embed at columns (N-1)*ns+1:N*ns.
    ph = h.phase
    zero_arg = () -> begin
        reg     = get(_phase_registry, objectid(ph), nothing)
        y0_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._y0) : ph._y0
        yf_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._yf) : ph._yf
        compact = fn(y0_s, yf_s, nothing, ph._t0, ph._tf, ph.model)
        ns = ph._n_states
        N  = n_control_nodes(ph)
        J  = zeros(size(compact, 1), ns * N)
        J[:, (N-1)*ns+1 : N*ns] .= compact
        J
    end
    add_jacobian!(h.bc.calc, ph.state_var, zero_arg)
end

function add_jacobian!(h::BoundaryConstraintHandle, ::FinalTime, fn::Function)
    ph = h.phase
    zero_arg = () -> fn(ph._tf)   # tf is scalar → returns (n_bc × 1)
    add_jacobian!(h.bc.calc, ph.tf_var, zero_arg)
end

function add_jacobian!(h::BoundaryConstraintHandle, ::InitialTime, fn::Function)
    ph = h.phase
    zero_arg = () -> fn(ph._t0)
    add_jacobian!(h.bc.calc, ph.t0_var, zero_arg)
end

# do-block sugar for all add_jacobian! overloads
add_jacobian!(fn::Function, h::BoundaryConstraintHandle, tag) =
    add_jacobian!(h, tag, fn)

# ─────────────────────────────────────────────────────────────────────────────
# set_objective! — new API
#   Returns an ObjectiveHandle; sense is Min() or Max()
# ─────────────────────────────────────────────────────────────────────────────

function set_objective!(phase::CollocationPhase, sense::Union{Min,Max})
    return ObjectiveHandle(phase, sense)
end

# ─────────────────────────────────────────────────────────────────────────────
# add_mayer! — new API
#   fn(y0, yf, p, t0, tf) -> scalar
#   Returns a MayerHandle for add_jacobian!
# ─────────────────────────────────────────────────────────────────────────────

# The value of a phase's Mayer term at its current decision vector, as a BolzaObjective holds it.
_mayer_value(ph, m::MayerObjective) = () -> begin
    ctx = BoundaryContext(_state_named(ph, ph._Y[:, 1]), _state_named(ph, ph._Y[:, end]),
                          ph._t0, ph._tf, ph._params)
    Float64(m.fn(ctx))
end

function add_mayer!(fn::Function, obj::ObjectiveHandle)
    h = MayerHandle(obj.phase)
    obj._mayer_ref[] = h

    # IMPORTANT: use ctx.y0/ctx.yf/ctx.t0/ctx.tf — not phase fields —
    # so ForwardDiff can perturb through them in the AD gradient path.
    ph = obj.phase
    sense_sym = obj.sense isa Min ? :Min : :Max
    internal_obj = MayerObjective(obj.phase; sense = sense_sym) do ctx
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            T    = eltype(ctx.y0)
            y0_s = _prepare_model!(reg, T, ctx.y0)   # target synced to y0
            yf_s = _prepare_model!(reg, T, ctx.yf)   # target synced to yf
            fn(y0_s, yf_s, nothing, ctx.t0, ctx.tf, ph.model)
        else
            fn(ctx.y0, ctx.yf, nothing, ctx.t0, ctx.tf, ph.model)
        end
    end
    # A Lagrange term declared first is kept: the phase becomes Bolza with this as its Mayer
    # half, carrying the running cost and its partials over unchanged.
    existing = ph.objective
    if existing isa BolzaObjective
        set_objective!(ph, BolzaObjective(existing.phases, existing.sense,
                                          _mayer_value(ph, internal_obj),
                                          existing.lagrange_fn, existing.jac_fns))
    else
        set_objective!(ph, internal_obj)
    end
    obj.phase._mayer_internal = internal_obj   # store for add_jacobian!
    return h
end

add_mayer!(obj::ObjectiveHandle, fn::Function) = add_mayer!(fn, obj)  # do-block

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! for MayerHandle — tag dispatch
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::MayerHandle, ::FinalState, fn::Function)
    # fn(y0, yf, p, t0, tf, model) -> (1 × ns) row matrix for a scalar Mayer term.
    # Embed into full (ns*N,) gradient vector: only yf slot is non-zero.
    ph = h.phase
    internal_obj = ph._mayer_internal
    zero_arg = () -> begin
        reg  = get(_phase_registry, objectid(ph), nothing)
        y0_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._y0) : ph._y0
        yf_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._yf) : ph._yf
        compact = vec(fn(y0_s, yf_s, nothing, ph._t0, ph._tf, ph.model))   # (ns,) vector
        ns = ph._n_states
        N  = n_control_nodes(ph)
        g  = zeros(ns * N)
        g[(N-1)*ns+1 : N*ns] .= compact
        g
    end
    add_objective_jacobian!(internal_obj, ph.state_var, zero_arg)
end

function add_jacobian!(h::MayerHandle, ::InitialState, fn::Function)
    ph = h.phase
    internal_obj = ph._mayer_internal
    zero_arg = () -> begin
        reg  = get(_phase_registry, objectid(ph), nothing)
        y0_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._y0) : ph._y0
        yf_s = (reg !== nothing && reg.state_type !== nothing) ?
                   _prepare_model!(reg, Float64, ph._yf) : ph._yf
        compact = vec(fn(y0_s, yf_s, nothing, ph._t0, ph._tf, ph.model))
        ns = ph._n_states
        N  = n_control_nodes(ph)
        g  = zeros(ns * N)
        g[1:ns] .= compact
        g
    end
    add_objective_jacobian!(internal_obj, ph.state_var, zero_arg)
end

function add_jacobian!(h::MayerHandle, ::FinalTime, fn::Function)
    internal_obj = h.phase._mayer_internal
    ph = h.phase
    zero_arg = () -> fn(ph._tf)   # scalar tf → [1.0], length-1 vector
    add_objective_jacobian!(internal_obj, ph.tf_var, zero_arg)
end

function add_jacobian!(h::MayerHandle, ::InitialTime, fn::Function)
    internal_obj = h.phase._mayer_internal
    ph = h.phase
    zero_arg = () -> fn(ph._t0)
    add_objective_jacobian!(internal_obj, ph.t0_var, zero_arg)
end

add_jacobian!(fn::Function, h::MayerHandle, tag) = add_jacobian!(h, tag, fn)

# ─────────────────────────────────────────────────────────────────────────────
# add_lagrange! — new API
#   lag = add_lagrange!(obj) do y, u, p, t ... end
#   fn(y::StateType, u::ControlType, p, t) -> scalar
#   If add_mayer! was already called, upgrades the phase objective to BolzaObjective.
#   If not, installs a pure-Lagrange BolzaObjective (zero Mayer term).
#   Returns a LagrangeHandle for add_jacobian!.
# ─────────────────────────────────────────────────────────────────────────────

function add_lagrange!(fn::Function, obj::ObjectiveHandle)
    ph        = obj.phase
    sense_sym = obj.sense isa Min ? :Min : :Max
    reg       = get(_phase_registry, objectid(ph), nothing)

    # Wrap user fn(y_struct, u_struct, p, t) -> scalar
    # into BolzaObjective's expected lagrange_fn(y_vec, u_vec, t) -> scalar.
    # Must NOT cast to Float64 — ForwardDiff calls this with Dual arrays.
    lagrange_wrapped = function _lag_integrand(y_vec, u_vec, t)
        T = promote_type(eltype(y_vec), eltype(u_vec))
        if reg !== nothing && reg.state_type !== nothing
            y_s, u_s = _prepare_model!(reg, T, y_vec, u_vec)
            if reg.control_type !== nothing
                return fn(y_s, u_s, ph.model, t)
            else
                return fn(y_s, u_vec, ph.model, t)   # control-free: pass raw u_vec (empty)
            end
        else
            return fn(y_vec, u_vec, ph.model, t)
        end
    end

    # Mayer term: reuse existing internal MayerObjective fn if present, else zero
    existing = ph._mayer_internal
    mayer_wrapped = existing === nothing ? () -> 0.0 : _mayer_value(ph, existing)

    # The Mayer half's partials stay on the MayerObjective, where objective_gradient_chunk
    # reads them; copying them into this table would count them twice.
    bolza = BolzaObjective(mayer_wrapped, lagrange_wrapped, ph; sense = sense_sym)

    set_objective!(ph, bolza)

    h = LagrangeHandle(ph, fn)
    obj._lagrange_ref[] = h
    return h
end

add_lagrange!(obj::ObjectiveHandle, fn::Function) = add_lagrange!(fn, obj)   # do-block

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! for LagrangeHandle — tag dispatch
#   fn(y_struct) -> Vector  OR  fn(u_struct) -> Vector  (∂L/∂var at one node)
#   Wraps into BolzaObjective jac_fns as a zero-arg closure returning (nc*N,) vector.
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::LagrangeHandle, ::Control, fn::Function)
    ph  = h.phase
    reg = get(_phase_registry, objectid(ph), nothing)
    # BolzaObjective gradient_chunk calls what was registered, returning a Vector
    # The Bolza AD chunk assembles node-by-node; here we register a zero-arg closure
    # that returns the full (nc*N,) gradient assembled with quadrature weights.
    # A partial takes the same arguments as the function it comes from. The
    # integrand is called fn(y, u, model, t), so its partial is too. This used
    # to pass only the block being differentiated against, which no
    # context-shaped partial could accept.
    zero_arg = () -> begin
        mesh  = build_mesh(ph.transcription)
        N     = mesh.N
        nc    = ph._n_controls
        qs    = quadrature_weights(mesh, ph._t0, ph._tf)
        ts    = node_times(mesh, ph._t0, ph._tf)
        g     = zeros(nc * N)
        for k in 1:N
            if reg !== nothing && reg.state_type !== nothing
                y_s, u_s = _prepare_model!(reg, Float64, ph._Y[:, k], ph._U[:, k])
                gk  = vec(fn(y_s, u_s, ph.model, ts[k]))
            else
                gk  = vec(fn(ph._Y[:, k], ph._U[:, k], ph.model, ts[k]))
            end
            g[(k-1)*nc+1 : k*nc] .= qs[k] .* gk
        end
        g
    end
    obj = ph.objective
    obj isa BolzaObjective && add_objective_jacobian!(obj, ph.control_var, zero_arg)
    nothing
end

function add_jacobian!(h::LagrangeHandle, ::State, fn::Function)
    ph  = h.phase
    reg = get(_phase_registry, objectid(ph), nothing)
    # A partial takes the same arguments as the function it comes from. The
    # integrand is called fn(y, u, model, t), so its partial is too. This used
    # to pass only the block being differentiated against, which no
    # context-shaped partial could accept.
    zero_arg = () -> begin
        mesh  = build_mesh(ph.transcription)
        N     = mesh.N
        ns    = ph._n_states
        qs    = quadrature_weights(mesh, ph._t0, ph._tf)
        ts    = node_times(mesh, ph._t0, ph._tf)
        g     = zeros(ns * N)
        for k in 1:N
            if reg !== nothing && reg.state_type !== nothing
                y_s, u_s = _prepare_model!(reg, Float64, ph._Y[:, k], ph._U[:, k])
                gk  = vec(fn(y_s, u_s, ph.model, ts[k]))
            else
                gk  = vec(fn(ph._Y[:, k], ph._U[:, k], ph.model, ts[k]))
            end
            g[(k-1)*ns+1 : k*ns] .= qs[k] .* gk
        end
        g
    end
    obj = ph.objective
    obj isa BolzaObjective && add_objective_jacobian!(obj, ph.state_var, zero_arg)
    nothing
end

add_jacobian!(fn::Function, h::LagrangeHandle, tag) = add_jacobian!(h, tag, fn)


#   fn(y::StateType, u::ControlType, p, t) -> Vector   (functional, not mutating)
#   Returns a PathConstraintHandle for add_jacobian!
# ─────────────────────────────────────────────────────────────────────────────

function add_path_constraint!(fn::Function, phase::CollocationPhase;
                               name::Symbol,
                               equality     = nothing,
                               lower_bounds = nothing,
                               upper_bounds = nothing)
    if equality !== nothing
        lb = Float64.(equality isa Number ? [equality] : equality)
        ub = copy(lb)
    else
        lb = Float64.(lower_bounds isa Number ? [lower_bounds] : lower_bounds)
        ub = Float64.(upper_bounds isa Number ? [upper_bounds] : upper_bounds)
    end
    n_pc = length(lb)
    ph   = phase

    # Wrap functional fn(y_struct, u_struct, p, t, model) into internal mutating form
    # fn!(g, y_vec, ctx::EvalContext, t)
    # Uses eltype(y_vec)/eltype(ctx.u) so ForwardDiff AD passes through correctly.
    wrapped! = function _pc_wrapper!(g, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            g .= fn(y_s, u_s, nothing, t, ctx.model)
        else
            fn(g, y_vec, ctx, t)
        end
    end

    pc = PathConstraint(wrapped!, n_pc, lb, ub, Dict{UInt64,Function}(), string(name))
    add_path_constraint!(phase, pc)
    return PathConstraintHandle(pc, phase)
end

# do-block sugar
add_path_constraint!(phase::CollocationPhase, args...; kwargs...) =
    add_path_constraint!(identity, phase, args...; kwargs...)   # placeholder — do-block form is primary

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! for PathConstraintHandle — tag dispatch
#   fn(u_struct) -> (n_pc × nc) matrix  OR  fn(y_struct) -> (n_pc × ns) matrix
#   Wraps into internal mutating fn!(dg, y_vec, ctx, t) signature.
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::PathConstraintHandle, ::Control, fn::Function)
    ph  = h.phase
    var = ph.control_var
    wrapped! = function _pc_ctrl_jac!(dg, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing && reg.control_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            dg .= fn(y_s, u_s, nothing, t, ctx.model)
        else
            fn(dg, y_vec, ctx, t)
        end
    end
    set_path_jacobian!(h.pc, var, wrapped!)
end

function add_jacobian!(h::PathConstraintHandle, ::State, fn::Function)
    ph  = h.phase
    var = ph.state_var
    wrapped! = function _pc_state_jac!(dg, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            dg .= fn(y_s, u_s, nothing, t, ctx.model)
        else
            fn(dg, y_vec, ctx, t)
        end
    end
    set_path_jacobian!(h.pc, var, wrapped!)
end

add_jacobian!(fn::Function, h::PathConstraintHandle, tag) =
    add_jacobian!(h, tag, fn)

# get_final_state is defined once, with get_initial_state, in collocation.jl. It
# used to be defined a second time here, and being the later definition it was
# the one that ran. The two agree whenever the phase has a registered state type,
# which is every phase a spec builds; they differ only for a phase built without
# one, where this version threw a KeyError and the other returns the raw vector.

# ─────────────────────────────────────────────────────────────────────────────
# CollocationPhase extension — storage slot for _mayer_internal
# Added here rather than in the struct definition to avoid breaking existing code.
# ─────────────────────────────────────────────────────────────────────────────

# Extend phase with a side-channel dict for new-API internal objects
const _phase_internals = Dict{UInt64, Dict{Symbol,Any}}()

function _phase_set!(phase::CollocationPhase, key::Symbol, val)
    d = get!(() -> Dict{Symbol,Any}(), _phase_internals, objectid(phase))
    d[key] = val
end

function _phase_get(phase::CollocationPhase, key::Symbol)
    d = get(_phase_internals, objectid(phase), nothing)
    d === nothing && return nothing
    get(d, key, nothing)
end

# Redirect _mayer_internal field accesses to the side-channel
function Base.getproperty(phase::CollocationPhase, name::Symbol)
    if name === :_mayer_internal
        return _phase_get(phase, :_mayer_internal)
    end
    return getfield(phase, name)
end

function Base.setproperty!(phase::CollocationPhase, name::Symbol, val)
    if name === :_mayer_internal
        _phase_set!(phase, :_mayer_internal, val)
        return val
    end
    return setfield!(phase, name, val)
end
