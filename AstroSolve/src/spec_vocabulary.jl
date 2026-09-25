# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The vocabulary a problem is written in.
#
# Vary says what a solver may change, Constraint what must hold, Objective what
# to minimise, Link how two phases join, and solve! runs whichever machinery the
# problem turns out to need. One set of words for an event graph, a transcribed
# phase and an estimation arc, because a mission analyst should not have to know
# which subsystem is answering.


# `solve!` takes a method per problem kind. The sequence and estimation ones are
# below; a phase type in another package adds its own the same way.

# Methods on the existing generics, not rivals to them. This is what makes the
# claim literal: a script calls the same `Vary` and the same `Constraint`
# whether its subject is a maneuver or a collocation phase.


# --- where a constraint is evaluated -----------------------------------------

"""
    Initial()

Places a phase constraint at the start of the arc: `Constraint(f, phase; equals = x0, at =
Initial())` evaluates `f` at the initial state and time. See also [`Final`](@ref),
[`Path`](@ref) and [`Boundary`](@ref).

# Example
<!-- doc-fragment -->
```julia
Constraint(position, phase; equals = [0.0], at = Initial())
```
"""
struct Initial   end

"""
    Final()

Places a phase constraint or objective at the end of the arc, evaluated at the final state and
time. An `Objective` at `Final()`, the default, is a Mayer term. See also [`Initial`](@ref) and
[`Path`](@ref).

# Example
<!-- doc-fragment -->
```julia
Constraint(altitude, phase; equals = [185.0], at = Final())
```
"""
struct Final     end

"""
    Path()

Places a phase constraint at every mesh node, so it holds along the whole arc. An `Objective` at
`Path()` is a Lagrange term, integrated over the arc. See also [`Initial`](@ref) and
[`Final`](@ref).

# Example
<!-- doc-fragment -->
```julia
Constraint(thrust_ball, phase; upper_bound = [1.0], at = Path())
```
"""
struct Path end

"""
    Boundary()

The default `at` for a phase `Constraint`. A function taking one argument is evaluated at the
final point, as with [`Final`](@ref). A function taking six, `(y0, yf, p, t0, tf, model)`,
receives both ends of the arc at once, which is how a constraint relating the start to the end
is written. Partials declared with `@partial` are not used at `Boundary()`; the framework
differentiates the function itself. Give `at = Final()` for a declared partial to be used.

# Example
<!-- doc-fragment -->
```julia
# a function of both ends, which Boundary() passes the whole arc
period_closes(y0, yf, p, t0, tf, model) = [yf[1] - y0[1], yf[2] - y0[2]]

Constraint(period_closes, phase; equals = [0.0, 0.0])
```
"""
struct Boundary  end

# The state and control types a phase was built with. When the phase took them
# in its constructor (proposal C) they are already in the framework's registry;
# the two-step form records them here as it goes.
const _CONTROL_TYPE = IdDict{Any,Any}()

function _state_type(p)
    r = get(_phase_registry, objectid(p), nothing)
    r === nothing || r.state_type === nothing ?
        throw(ArgumentError(
            "phase :$(p.name) must be given a state type; pass `state = YourState` to " *
            "CollocationPhase")) : r.state_type
end

function _control_type(p)
    r = get(_phase_registry, objectid(p), nothing)
    r !== nothing && r.control_type !== nothing && return r.control_type
    haskey(_CONTROL_TYPE, p) && return _CONTROL_TYPE[p]
    throw(ArgumentError(
        "phase :$(p.name) must be given a control type to vary a control; pass " *
        "`control = YourControl` to CollocationPhase"))
end

_vec(x) = collect(float.(x))
_vec(x::Number) = [float(x)]

# =============================================================================
# Vary — collocation
# =============================================================================

# Every `Vary` accepts the same six keywords: guess, lower_bound, upper_bound,
# scale, shift, name. A subsystem that cannot use one takes it and says so,
# rather than leaving it out of the signature — a missing keyword raises a
# MethodError naming the wrong problem, which is how the `scale` gap
# survived for months. `test_vary_contract.jl` checks this by reflection.

"Name the variable a setter just built and stored on the phase."
_named!(v, name) = (isempty(name) || (v.name = String(name)); v)

"Reject a keyword this subject genuinely has no use for."
_no_kw(kw, subject, why) = throw(ArgumentError(
    "Vary: `$kw` has no meaning for $subject — $why"))

# `scale` is not decoration. In dimensional units a state vector mixes 1e4 km
# with 1e3 kg with 7 km/s, and an unscaled NLP spends its iterations on the
# conditioning rather than the problem.
function Vary(::typeof(state), p::CollocationPhase, T::Type;
              guess, lower_bound, upper_bound, scale = nothing, shift = nothing,
              name = "")
    set_state!(p, T; guess = guess, lower_bounds = _vec(lower_bound),
               upper_bounds = _vec(upper_bound), scale = scale, shift = shift)
    return _named!(p.state_var, name)
end

function Vary(::typeof(control), p::CollocationPhase, T::Type;
              guess, lower_bound, upper_bound, scale = nothing, shift = nothing,
              name = "")
    _CONTROL_TYPE[p] = T
    set_control!(p, T; guess = guess, lower_bounds = _vec(lower_bound),
                 upper_bounds = _vec(upper_bound), scale = scale, shift = shift)
    return _named!(p.control_var, name)
end

# Free time. The constructor's tspan was the guess; this makes it a variable,
# which is what Vary means everywhere else too.
"""Static parameters: constant over the arc, one column each, no node index.

Estimating a drag coefficient or identifying a plant constant is this, and it is
the same line as any other Vary.
"""
function Vary(::typeof(parameter), p::CollocationPhase; lower_bound, upper_bound,
              guess = nothing, scale = nothing, shift = nothing, name = "")
    set_parameter!(p; guess        = guess,
                      lower_bounds = _vec(lower_bound),
                      upper_bounds = _vec(upper_bound),
                      scale        = scale, shift = shift)
    return _named!(p.param_var, name)
end

function Vary(::typeof(initial_time), p::CollocationPhase; lower_bound, upper_bound,
              guess = nothing, scale = nothing, shift = nothing, name = "")
    set_initial_time!(p, guess === nothing ? p._t0 : float(guess);
                      lower_bounds = float(lower_bound),
                      upper_bounds = float(upper_bound),
                      scale = scale === nothing ? 1.0 : float(scale),
                      shift = shift === nothing ? 0.0 : float(shift))
    return _named!(p.t0_var, name)
end

function Vary(::typeof(final_time), p::CollocationPhase; lower_bound, upper_bound,
              guess = nothing, scale = nothing, shift = nothing, name = "")
    set_final_time!(p, guess === nothing ? p._tf : float(guess);
                    lower_bounds = float(lower_bound),
                    upper_bounds = float(upper_bound),
                    scale = scale === nothing ? 1.0 : float(scale),
                    shift = shift === nothing ? 0.0 : float(shift))
    return _named!(p.tf_var, name)
end

# =============================================================================
# Vary — shooting
#
#= gap: a shooting phase takes a DirectSolverVariable and nothing else, so the
   per-component bounds a user writes have to be repeated here, the guess
   flattened here, and the phase's working arrays filled here. On the
   collocation side `set_state!` does all three from the state type. =#
# =============================================================================

# `Vary` forwards `scale` to every phase type. Some methods did not accept
# it, so any two-argument `Vary(state, zoh_phase; ...)` — which forwards through
# the fallback below and always passes `scale` — failed with a MethodError.
# The variable already carries an explicit scale vector, so honour it.
_scale_vec(scale, n) = scale === nothing ? ones(n) :
                       scale isa Number  ? fill(float(scale), n) : _vec(scale)
_shift_vec(shift, n) = shift === nothing ? zeros(n) :
                       shift isa Number  ? fill(float(shift), n) : _vec(shift)

# `Vary(initial_time, p; equals = t)` is gone. A fixed interval is the
# constructor's `tspan`; `Vary` means the solver may change it. The form had
# been dead anyway: for CollocationPhase the same positional signature taking
# bounds was defined below it, and Julia does not dispatch on keywords, so the
# later definition had silently replaced it. That phase type now has no
# free-time method at all, which no use case needs; add one taking bounds when
# one does.

# --- the phase already knows its types ---------------------------------------
#
# When a phase was constructed complete, `Vary` says only what the solver may
# change. The type is not repeated because the phase has it.

# The fallback below takes any subject, so a spacecraft would match it and be
# asked for bounds a phase needs and an estimated state does not have. A
# spacecraft goes to AstroSolve's own Vary, named through invoke because this
# method is the more specific one and would otherwise call itself.
Vary(q::typeof(state), sc::Spacecraft; kwargs...) =
    invoke(AstroSolve.Vary, Tuple{Function, Any}, q, sc; kwargs...)

Vary(q::typeof(state), p; guess, lower_bound, upper_bound, scale = nothing,
     shift = nothing, name = "") =
    Vary(q, p, _state_type(p); guess = guess, lower_bound = lower_bound,
         upper_bound = upper_bound, scale = scale, shift = shift, name = name)

Vary(q::typeof(control), p; guess, lower_bound, upper_bound, scale = nothing,
     shift = nothing, name = "") =
    Vary(q, p, _control_type(p); guess = guess, lower_bound = lower_bound,
         upper_bound = upper_bound, scale = scale, shift = shift, name = name)

# =============================================================================
# Constraint
# =============================================================================

_bounds(x) = x isa Number ? [float(x)] : _vec(x)

# `equals` is the common case and sets both bounds; give lower_bound and
# upper_bound instead for an inequality. Same words as every other spec.
function Constraint(f::Function,
                    p::Union{CollocationPhase,AbstractShootingPhase};
                    equals = nothing, lower_bound = nothing, upper_bound = nothing,
                    scale = nothing, at = Boundary(), name = gensym(:con))
    if equals !== nothing
        (lower_bound === nothing && upper_bound === nothing) || throw(ArgumentError(
            "Constraint: give `equals`, or `lower_bound`/`upper_bound`, not both."))
        return _con(f, p, _bounds(equals), _bounds(equals), at, name, scale)
    end
    (lower_bound === nothing && upper_bound === nothing) && throw(ArgumentError(
        "Constraint needs `equals`, or a `lower_bound` or `upper_bound`."))
    n  = length(_bounds(something(lower_bound, upper_bound)))
    lb = lower_bound === nothing ? fill(-Inf, n) : _bounds(lower_bound)
    ub = upper_bound === nothing ? fill( Inf, n) : _bounds(upper_bound)
    return _con(f, p, lb, ub, at, name, scale)
end

# check_partials needs two things it cannot recover afterwards: the user's
# function (the phase keeps only the framework's wrapper) and the `at` that
# says which endpoint a boundary partial with respect to `state` refers to.
# `@partial(f, state)` at Initial and at Final are the same declaration; only
# the constraint distinguishes them. REVIEW 4.2 called this out as implicit.
const _phase_specs = Dict{UInt64, Vector{Any}}()

_record_spec!(p, f, at; kind::Symbol = :constraint) =
    push!(get!(() -> Any[], _phase_specs, objectid(p)), (f = f, at = at, kind = kind))

"A partial's owner reads better without the closure's `#` prefix."
_fname(f) = lstrip(_pname(f), '#')

phase_specs(p) = get(_phase_specs, objectid(p), ())

# Constraint scaling on a collocation phase. The solver sees each component divided by its scale,
# with its bounds divided the same way, which is the convention the shooting phases use for theirs
# (g_nlp = g / scale). The function and every declared partial are wrapped, so the analytic path
# and the differentiated one see the same rows. check_partials compares the function the user
# wrote, which is unscaled, so it is unaffected.
function _con_scale(scale, n)
    scale === nothing && return nothing
    s = scale isa Number ? fill(float(scale), n) : _vec(scale)
    length(s) == n || throw(ArgumentError(
        "Constraint: `scale` has $(length(s)) entries for a constraint with $n components; give " *
        "one per component, or a single number for all of them."))
    all(>(0), s) || throw(ArgumentError(
        "Constraint: every `scale` entry must be positive; got $s."))
    return s
end

_scaled(g, ::Nothing) = g
_scaled(g, s::AbstractVector) = (args...) -> g(args...) ./ s

function _con(f, p::CollocationPhase, lb, ub, at::Path, name, scale = nothing)
    _check_arity(f, at)
    _record_spec!(p, f, at)
    s = _con_scale(scale, length(lb))
    h = add_path_constraint!(_scaled(_ctx_path(f, p), s), p; name = Symbol(name),
                             lower_bounds = _scaled_bound(lb, s), upper_bounds = _scaled_bound(ub, s))
    _attach!(h, f, at, s)
    return h
end

function _con(f, p::CollocationPhase, lb, ub, at::Union{Initial,Final,Boundary}, name, scale = nothing)
    _record_spec!(p, f, at)
    s = _con_scale(scale, length(lb))
    h = add_boundary_constraint!(_scaled(_as_boundary(_ctx_bnd(f, p, at), at, p), s), p;
                                 name = Symbol(name), lower_bounds = _scaled_bound(lb, s),
                                 upper_bounds = _scaled_bound(ub, s))
    _attach!(h, f, at, s)
    return h
end

_scaled_bound(b, ::Nothing) = b
_scaled_bound(b, s) = b ./ s

#= gap: shooting phases have no add_path_constraint!. A per-node constraint has
   to be lowered to one boundary constraint that loops the segments itself —
   which is exactly the hand-written `unit_thrust_ctx` this was meant to
   remove. The user no longer writes the loop; the framework still needs one. =#

# --- path constraints on a Sims-Flanagan phase --------------------------------
#
# The physical limit on a throttle is |u| <= 1, a ball. `Vary` gives a box,
# which admits sqrt(3) at the corners — so before this, uc12 and uc13 were both
# handing the solver 73% more thrust than the vehicle has.
#
# A shooting phase has no `add_path_constraint!`, so the constraint is lowered
# to one boundary constraint that evaluates the user's function once per
# segment. The user does not write the loop; the framework does.
#
#     thrust_ball(c) = [dot(control(c), control(c))]
#     Constraint(thrust_ball, sf; upper_bound = 1.0, at = Path())
#
# The function takes a context, as everything else on the shooting side does.
# When SPEC section 3's second half settles on one calling convention, this
# comes along with the rest of that side.

"Block-diagonal Jacobian of a per-segment function: segment k touches only its own three columns."
function _sf_path_jac(f, u, half::Symbol, m::Int, nseg::Int, rows::Int, rowoff::Int)
    J = zeros(rows, 3 * nseg)
    for k in 1:nseg
        g = ForwardDiff.jacobian(uu -> f(SFPathContext(uu, k, half)), collect(u[:, k]))
        J[rowoff + (k - 1) * m .+ (1:m), (k - 1) * 3 .+ (1:3)] = g
    end
    return J
end

function _con(f, p::SimsFlanaganPhase, lb, ub, ::Path, name, scale = nothing)
    _record_spec!(p, f, Path())
    nf, nb = n_fwd(p), n_bwd(p)
    m      = length(_vec(lb))
    rows   = m * (nf + nb)

    (p.u_fwd_var === nothing || p.u_bwd_var === nothing) && throw(ArgumentError(
        "Constraint(..., at = Path()) on a Sims-Flanagan phase needs both control " *
        "blocks to be variables first — `Vary(forward_control, phase; ...)` and " *
        "`Vary(backward_control, phase; ...)`."))

    looped(ctx) = vcat(
        reduce(vcat, (_vec(f(SFPathContext(view(ctx.u_fwd, :, k), k, :forward)))
                      for k in 1:nf); init = Float64[]),
        reduce(vcat, (_vec(f(SFPathContext(view(ctx.u_bwd, :, k), k, :backward)))
                      for k in 1:nb); init = Float64[]))

    h = add_boundary_constraint!(looped, p; name = Symbol(name),
                                 lower_bounds = repeat(_vec(lb), nf + nb),
                                 upper_bounds = repeat(_vec(ub), nf + nb))

    # Read the throttles at evaluation time — the phase keeps them current.
    add_jacobian!(h, p.u_fwd_var,
                  _ -> _sf_path_jac(f, p._u_fwd, :forward,  m, nf, rows, 0))
    add_jacobian!(h, p.u_bwd_var,
                  _ -> _sf_path_jac(f, p._u_bwd, :backward, m, nb, rows, m * nf))
    return h
end

# An MGAnDSMs constraint is handed the phase's own boundary context, so the
# user's function reads it with the same quantity names Vary declares.
"""
Attach an MGAnDSMs partial, handing it the context its constraint reads.

The framework calls a Jacobian closure with the variable block alone, which is
enough only when the derivative depends on nothing else. It usually does — the
original script has partials reaching out to `arr2_var.value` and
`dep2_var.value`, module-level variable objects, because the closure was not
given anything better. That breaks with two missions in one session and
couples a derivative to a global name.

The context is reconstructible from the phase, which the framework keeps
current, so the user writes `do ctx` — the same argument the constraint takes
— and reads it with the same quantity names `Vary` declares.
"""
function _attach_mga!(h, f, p)
    for (q, jac) in partials_of(f)
        # A declared partial names the quantity it differentiates. If that
        # quantity was never varied the slot is `nothing`, and `variable_list`
        # builds itself with filter(!isnothing, ...) — so the derivative is
        # registered against nothing and the NLP is quietly one variable
        # smaller. Say which quantity, here, rather than letting add_jacobian!
        # report a MethodError on ::Nothing.
        slot_of(p, q) === nothing && throw(ArgumentError(
            "`$(_fname(f))` declares a partial with respect to " *
            "`$(nameof(q))`, but `$(nameof(q))` is not a variable on phase " *
            "$(p.name). Add `Vary($(nameof(q)), phase; ...)`, or drop the " *
            "partial — as written the derivative has nothing to act on."))
        add_jacobian!(h, slot_of(p, q), _ -> jac(_mga_boundary_context(p)))
    end
    return nothing
end

# Constraint scaling is not decoration either: a mass residual near 1000 and a
# C3 residual near 20 in the same NLP is the conditioning problem again.
function _con(f, p::MGAnDSMsPhase, lb, ub, at, name, scale = nothing)
    _record_spec!(p, f, at)
    h = scale === nothing ?
        add_boundary_constraint!(f, p; name = Symbol(name),
                                 lower_bounds = lb, upper_bounds = ub) :
        add_boundary_constraint!(f, p; name = Symbol(name), lower_bounds = lb,
                                 upper_bounds = ub, scale = _vec(scale))
    _attach_mga!(h, f, p)
    return h
end

# =============================================================================
# A ForceModel as a phase's right-hand side
#
# The propagator drives forces through `_build_odes!`, which allocates a
# 6-vector, calls `accel_eval!(force, t, posvel, acc, sc, p)` for each force,
# and writes six components. A transcription needs the same loop with three
# differences: the state may be longer than six (mass), the control has to
# reach the forces, and `t` is seconds from the phase start rather than an
# epoch.
#
# Nothing here reimplements a force. The state slice a quantity cares about is
# still 1:6, and anything past it is the user's.
# =============================================================================

"""
    forcemodel_rhs(forces, epoch0, nstates) -> f(dy, y, u, p, t, model)

Wrap a ForceModel so a phase can be flown by it.

`epoch0` is the phase's start epoch; `t` arrives in seconds from there.
The control is handed to the forces through the `params` slot the
`accel_eval!` interface already carries — it is the one thing a force needs
that lives nowhere else.
"""
function forcemodel_rhs(forces, epoch0, nstates::Int, sc)
    return function _fm_rhs!(dy, y, u, p, t, model)
        x = _flat(y, nstates)
        T = eltype(x)
        acc = zeros(T, nstates)
        epoch = epoch0 + t / 86400.0
        # Each force gets a fresh buffer and everything past the kinematic rows is summed. Built-in
        # forces overwrite what they write, so one shared buffer kept only the last force; a user
        # force that adds, like a thruster's mass flow, is summed correctly either way.
        @inbounds for i in 1:3
            dy[i] = x[i + 3]
        end
        @inbounds for i in 4:nstates
            dy[i] = zero(T)
        end
        for f in forces.forces
            fill!(acc, zero(T))
            accel_eval!(f, epoch, x, acc, sc, (control = u, params = p))
            @inbounds for i in 4:nstates
                dy[i] += acc[i]
            end
        end
        return nothing
    end
end

"""
A phase handed a ForceModel wraps it. The spacecraft supplies the start epoch
and is passed through to the forces; the state length comes from the state
type the constructor already derived.
"""
function set_dynamics!(phase::CollocationPhase, fm::ForceModel; model = nothing)
    model === nothing && throw(ArgumentError(
        "a ForceModel needs a spacecraft: give the phase `model = sat`, which " *
        "supplies the start epoch and is handed to each force."))
    f = forcemodel_rhs(fm, model.time, phase._n_states, model)
    set_dynamics!(phase, f; model = model)
end

"""Flatten a user state struct into the vector a force expects."""
_flat(y::AbstractVector, n) = y
_flat(y, n) = [getfield(y, i) for i in 1:n]

# =============================================================================
# Partials
#
# The framework has two calling conventions — a path function takes
# (y, u, p, t, model) and a boundary function takes (y0, yf, p, t0, tf, model).
# A user who writes one function and constrains it in both places should not
# have to write it twice, so the arity decides and the shorter one is presented
# at whichever endpoint it was asked for.
# =============================================================================

_arity(f) = minimum(m.nargs - 1 for m in methods(f))

"""
    subject_at(phase, y, t) -> Coordinate

Turn a phase endpoint into something an Epicycle quantity can read.

A quantity wants a subject with a state, a frame and an epoch. A phase
endpoint is a bare vector of numbers, but the phase's spacecraft already
carries the other two — so the frame is not a separate thing to declare, it is
the frame the spacecraft's state is in.

The first six components are the Cartesian state; anything past them is the
user's, and no quantity asks for it.
# Returns
A `Coordinate` holding the first six components as a Cartesian state, in the phase's coordinate
system, at the phase's epoch plus `t`.

# Example
<!-- doc-fragment -->
```julia
subject_at(phase, get_final_state(phase), get_final_time(phase))
```
"""
subject_at(phase, y, t) =
    Coordinate(CartesianState(collect(_first6(y))),
               phase.model.coord_sys,
               phase.model.time + t / 86400.0)

"""
    subject_at(phase::CollocationPhase, Final()) -> Coordinate
    subject_at(phase::CollocationPhase, Initial()) -> Coordinate

The spacecraft state at the end or the start of a solved phase, as a subject any Epicycle
quantity reads: `semi_major_axis(subject_at(phase, Final()))`.

# Notes
Applies to a phase flown by a force model, whose `model` is the spacecraft: the frame is the
spacecraft's coordinate system, and the epoch is the spacecraft's epoch plus the phase time,
which is in seconds. The first six state components are the Cartesian position and velocity,
in km and km/s.

# Returns
A `Coordinate` holding the Cartesian state, its frame and its epoch.
"""
subject_at(phase::CollocationPhase, ::Final)   = subject_at(phase, phase._Y[:, end], phase._tf)
subject_at(phase::CollocationPhase, ::Initial) = subject_at(phase, phase._Y[:, 1], phase._t0)

_first6(y::AbstractVector) = y[1:6]
_first6(y) = (getfield(y, i) for i in 1:6)

# `at` decides the calling convention; arity is checked against it rather than
# deciding it. That was the silent failure — a boundary function that happened
# to take five arguments was quietly evaluated at every node. Now it says so.
#
# `at` cannot separate the two remaining cases, because `at = Final()` fits both
# an Epicycle quantity and a user boundary function. That one stays on arity,
# but it is a narrow test with a loud error rather than a silent fallthrough.
# An Epicycle quantity and a context-shaped spec function both take one
# argument, so arity stopped telling them apart the moment collocation
# functions became `f(c)`. Method signatures cannot separate them either: a
# quantity like `eccentricity` accepts any subject, so it happily accepts a
# PhaseContext and fails deep inside `frame_of`.
#
# Ask the trait. A quantity declares a label — "Eccentricity", "Position
# magnitude" — where everything else gets the fallback. That is what the trait
# is for, and it is how the vocabulary already answers `setter` and
# `is_settable`.
_is_quantity(f) = EpicycleBase.label(f) != "quantity"

const _EXPECTED_ARITY = Dict(:path => 5, :boundary => 6)

# The check is one-directional, because the two mistakes are not symmetric.
# A path-shaped function CAN be evaluated at an endpoint — that is the wrap
# `_as_boundary` does, and it is what lets uc5 constrain one `speed` at a node
# and at the final time. A boundary-shaped function cannot be evaluated at a
# node: there is no y0 and yf there, only one state. So only that direction is
# an error.
function _check_arity(f, at)
    (_is_quantity(f) || !(at isa Path)) && return nothing
    got = _arity(f)
    got <= 5 && return nothing
    throw(ArgumentError(
        "$(nameof(typeof(f))) takes $got arguments, so it is a boundary " *
        "function — (y0, yf, p, t0, tf, model). `at = Path()` evaluates at a " *
        "single node, which has one state, and calls with " *
        "(y, u, p, t, model). Use at = Initial() or at = Final()."))
end

_as_boundary(f, at, phase) =
    (_check_arity(f, at); _as_boundary(f, at, phase, Val(_is_quantity(f))))

_as_boundary(f, at::Initial, phase, ::Val{true}) =
    (y0, yf, p, t0, tf, model) -> [f(subject_at(phase, y0, t0))]
_as_boundary(f, at::Union{Final,Boundary}, phase, ::Val{false}) =
    _arity(f) == 6 ? f : (y0, yf, p, t0, tf, model) -> f(yf, nothing, p, tf, model)
_as_boundary(f, at::Union{Final,Boundary}, phase, ::Val{true}) =
    (y0, yf, p, t0, tf, model) -> [f(subject_at(phase, yf, tf))]
_as_boundary(f, at::Initial, phase, ::Val{false}) =
    _arity(f) == 6 ? f : (y0, yf, p, t0, tf, model) -> f(y0, nothing, p, t0, model)

# =============================================================================
# One calling convention
#
# Collocation unpacked positionally and in two shapes — (y, u, p, t, model) at a
# path, (y0, yf, p, t0, tf, model) at an endpoint — while MGAnDSMs and
# Sims-Flanagan handed the function a context. So what a user wrote depended on
# the transcription they picked and on where the constraint sat, neither of
# which is about their problem.
#
# A collocation function now takes a context too:
#
#     speed(c)          = [state(c).v]
#     final_position(c) = [state(c).x, state(c).y]
#     @partial(speed, state) do c; [0.0 0.0 1.0]; end
#
# `state(c)` is the state where the constraint is evaluated — the node at a
# path, the relevant end at an endpoint. That is what lets one function serve
# both, which is uc5's `speed` and the case this turns on. Every endpoint
# function in the suite uses exactly the end its `at` names, so no reader for
# "the other end" is needed; add one when something needs it.
#
# Positional functions still work. They are wrapped, not rejected, so the
# migration can go file by file against recorded numbers.
# =============================================================================

"""
    PhaseContext

What a collocation constraint, objective or partial is handed. Read it with
quantity names: `state(c)`, `control(c)`, `initial_time(c)`, `final_time(c)`.
The model and any parameters are fields — `c.model`, `c.p` — because `model` is
a name users give their own variables.
"""
# The three times carry independent type parameters. Under ForwardDiff only the
# time being differentiated arrives as a Dual, so tying them to one parameter
# left no matching constructor exactly when the gradient was being taken.
struct PhaseContext{Y,U,T,T0,TF,P,M}
    y     ::Y     # state where this is evaluated
    u     ::U     # control there; nothing at an endpoint
    t     ::T     # time there
    t0    ::T0
    tf    ::TF
    p     ::P
    model ::M
end

state(c::PhaseContext)        = c.y
control(c::PhaseContext)      = c.u
parameter(c::PhaseContext)    = c.p
initial_time(c::PhaseContext) = c.t0
final_time(c::PhaseContext)   = c.tf

# An Epicycle quantity also takes one argument, and it is not a spec function —
# wrapping `eccentricity` in a context handed it a PhaseContext and it failed
# inside `frame_of`. One argument AND not a declared quantity.
_takes_context(f) = _arity(f) == 1 && !_is_quantity(f)

"Wrap a context-shaped function so the framework can call it at a path node."
_ctx_path(f, ph) = _takes_context(f) ?
    ((y, u, p, t, m) -> f(PhaseContext(y, u, t, ph._t0, ph._tf, p, m))) : f

"Wrap a context-shaped function so the framework can call it at an endpoint."
function _ctx_bnd(f, ph, at)
    _takes_context(f) || return f
    first_end = at isa Initial
    return (y0, yf, p, t0, tf, m) -> f(PhaseContext(first_end ? y0 : yf, nothing,
                                                    first_end ? t0 : tf, t0, tf, p, m))
end

# The quantity says what the derivative is with respect to; the constraint's
# `at` says where it is taken. Neither repeats the other, and the framework's
# own tags stay internal.
# `s` is the constraint's scale, or nothing; a partial is scaled by row as its function is.
_attach!(h, f, ::Path, s = nothing) =
    for (q, jac) in partials_of(f)
        tag = _path_tag(q)
        tag === nothing || add_jacobian!(h, tag, _scaled(_ctx_path(jac, h.phase), s))
    end

function _attach!(h, f, at::Union{Initial,Final}, s = nothing)
    for (q, jac) in partials_of(f)
        tag = _bnd_tag(q, at isa Initial)
        tag === nothing && continue
        j = _ctx_bnd(jac, h.phase, at)
        add_jacobian!(h, tag, _scaled(q === state && !_takes_context(jac) ?
                                      _as_boundary(j, at, h.phase) : j, s))
    end
end

_attach!(h, f, ::Boundary, s = nothing) = nothing

# =============================================================================
# Objective
# =============================================================================

# A shooting phase's Mayer term takes its own context, so nothing is wrapped.
function Objective(f::Function, p::MGAnDSMsPhase; sense = Min())
    _record_spec!(p, f, Final(); kind = :objective)
    h = set_objective!(p, sense)
    m = add_mayer!(f, h)
    _attach_mga!(m, f, p)
    return h
end

# The same for a Sims-Flanagan phase, whose Mayer term reads its own boundary context.
function Objective(f::Function, p::SimsFlanaganPhase; sense = Min())
    _record_spec!(p, f, Final(); kind = :objective)
    h = set_objective!(p, sense)
    m = add_mayer!(f, h)
    for (q, jac) in partials_of(f)
        slot_of(p, q) === nothing && throw(ArgumentError(
            "`$(_fname(f))` declares a partial with respect to `$(nameof(q))`, but " *
            "`$(nameof(q))` is not a variable on phase $(p.name). Add " *
            "`Vary($(nameof(q)), phase; ...)`, or drop the partial."))
        add_jacobian!(m, slot_of(p, q), _ -> jac(_sf_boundary_context(p)))
    end
    return h
end

# Quantities a Sims-Flanagan objective or constraint reads from the phase's boundary context.
state(c::SFBoundaryContext)            = c.x0
final_state(c::SFBoundaryContext)      = c.xf
departure_vinf(c::SFBoundaryContext)   = c.vinf_dep
arrival_vinf(c::SFBoundaryContext)     = c.vinf_arr
initial_time(c::SFBoundaryContext)     = c.t0
final_time(c::SFBoundaryContext)       = c.tf
initial_mass(c::SFBoundaryContext)     = c.m0
final_mass(c::SFBoundaryContext)       = c.mf
forward_control(c::SFBoundaryContext)  = c.u_fwd
backward_control(c::SFBoundaryContext) = c.u_bwd

# `maximize = true` said what `sense = Max()` already said. One name.
"""Wrap a context-shaped function so the framework can integrate it along the arc.

`add_lagrange!` calls its integrand as `fn(y, u, model, t)`, which is a
different shape from the path-constraint call, so this is not `_ctx_path`.
"""
_ctx_lag(f, ph) = _takes_context(f) ?
    ((y, u, m, t) -> f(PhaseContext(y, u, t, ph._t0, ph._tf, ph._params, m))) : f

"""
    Objective(f, phase; sense = Min(), at = Final())

The cost. `at = Final()` is a Mayer term, read at the end of the arc. `at =
Path()` is a Lagrange term, integrated along it.

A Bolza objective is both, written as two calls, the same way `uc5` constrains
`speed` twice with different `at`:

```julia
Objective(terminal_cost, phase; sense = Min())
Objective(running_cost,  phase; sense = Min(), at = Path())
```

Order does not matter. Whichever comes second keeps what the first declared.

# Returns
The objective handle it attached to the phase. A caller reads the achieved value off the solve
result rather than off this handle, so it is rarely kept.
"""
function Objective(f::Function, p; sense = Min(), at = Final())
    at isa Path && return _objective_path(f, p, sense)
    _is_quantity(f) || _record_spec!(p, f, Final(); kind = :objective)
    s = sense
    h = set_objective!(p, s)
    # An Epicycle quantity can be the objective too, read off the final state
    # the same way a terminal constraint reads it — except a Mayer term is a
    # scalar where a constraint is a vector.
    # A Mayer term is evaluated at the final point, so a context-shaped
    # objective is wrapped the same way a Final() constraint is.
    fo = _is_quantity(f) ?
         ((y0, yf, pp, t0, tf, model) -> f(subject_at(p, yf, tf))) :
         _ctx_bnd(f, p, Final())
    m = add_mayer!(fo, h)
    for (q, jac) in partials_of(f)
        tag = _bnd_tag(q, false)
        tag === nothing && continue
        if tag isa Union{InitialTime,FinalTime}
            # The framework calls a time gradient with the scalar alone; a
            # partial takes the same arguments as the function it comes from.
            jt = _ctx_bnd(jac, p, Final())
            add_jacobian!(m, tag, t -> jt(p._y0, p._yf, nothing,
                                          tag isa InitialTime ? t : p._t0,
                                          tag isa FinalTime   ? t : p._tf, p.model))
        else
            add_jacobian!(m, tag, _ctx_bnd(jac, p, Final()))
        end
    end
    return h
end

"An integral term. The integrand is evaluated at every node and weighted."
function _objective_path(f::Function, p, sense)
    _record_spec!(p, f, Path(); kind = :objective)
    h = set_objective!(p, sense)
    lag = add_lagrange!(_ctx_lag(f, p), h)
    for (q, jac) in partials_of(f)
        tag = _path_tag(q)
        tag === nothing || add_jacobian!(lag, tag, _ctx_lag(jac, p))
    end
    return h
end


# =============================================================================
# check_partials — a wrong partial is worse than no partial
# =============================================================================

"""
    check_partials(phase; step = 1e-6, verbose = true)

Compare every declared partial on a phase against central differences.

Covers the dynamics, and every constraint attached to the phase — path
constraints, endpoint constraints, and the derivatives with respect to phase
times. A wrong partial does not announce itself: the solver converges somewhere
else, or stops converging at all, which is what makes it worth an explicit
check rather than a debugging session.

Which convention a declared partial uses is read off its own arity — six
arguments is an endpoint partial, five a path partial. What a boundary partial
with respect to `state` differentiates is *not* recoverable that way, because
`@partial(f, state)` is the same declaration at either end; the constraint's
`at` settles it.
# Returns
`nothing`. The comparison is printed as a table, one row per function and variable block, giving the
largest absolute and relative disagreement and whether the block fell back to automatic
differentiation.

# Example
```julia
check_partials(phase)
```
"""
function check_partials(p; step = 1e-6, verbose = true)
    reg = get(_phase_registry, objectid(p), nothing)
    reg === nothing && return nothing
    rows = Tuple{String,String,Float64}[]
    bare = String[]

    _check_dynamics!(rows, p, reg, step)
    for spec in phase_specs(p)
        isempty(partials_of(spec.f)) && push!(bare, _fname(spec.f))
        _check_constraint!(rows, p, reg, spec, step)
    end

    # A declared partial the framework never picked up is invisible otherwise:
    # the answers stay right and the run is only slower. So say which
    # derivatives are being used, not only whether the declared ones are right.
    if verbose
        if isempty(rows) && isempty(bare)
            println("  check_partials: nothing declared on this phase.")
        else
            @printf("  %-28s %-16s %-7s %s
",
                    "function", "with respect to", "source", "max|analytic - fd|")
            println("  ", "-"^70)
            for (fname, qname, w) in rows
                @printf("  %-28s %-16s %-7s %.3e
", fname, qname, "user", w)
            end
            for f in unique(bare)
                @printf("  %-28s %-16s %-7s %s
",
                        f, "every block", "AD", "not compared")
            end
            println()
            println(isempty(bare) ? "  nothing fell back to AD" :
                    "  fell back to AD: " * join(unique(bare), ", "))
        end
    end
    worst = isempty(rows) ? 0.0 : maximum(r[3] for r in rows)
    worst > 1e-4 && @warn "check_partials: a declared partial disagrees with " *
                          "finite differences by $(round(worst, sigdigits = 3))."
    return worst
end

"Central difference of `ev` in the jth component of `base`."
function _fd_column(ev, base, j, step)
    hj = step * max(abs(base[j]), 1.0)
    up = copy(base); up[j] += hj
    dn = copy(base); dn[j] -= hj
    return (ev(up) .- ev(dn)) ./ (2hj)
end

function _worst_vs_fd(A, ev, base, step)
    worst = 0.0
    for j in eachindex(base)
        fd = _fd_column(ev, base, j, step)
        worst = max(worst, maximum(abs.(A[:, j] .- fd)))
    end
    return worst
end

# --- the dynamics, in place, as before ---------------------------------------

function _check_dynamics!(rows, p, reg, step)
    f = get(_raw_dynamics, objectid(p), nothing)
    f === nothing && return nothing

    ns, nc = reg.n_states, reg.n_controls
    y0 = p.state_var   === nothing ? zeros(ns) : collect(p.state_var.value)[1:ns]
    u0 = p.control_var === nothing ? Float64[] : collect(p.control_var.value)[1:nc]
    # A control-free phase has no control type to build a struct from, and a
    # phase with no parameters must still be handed a vector rather than
    # `nothing`, because a parameter partial indexes it.
    p0 = p.param_var   === nothing ? Float64[] : collect(p.param_var.value)
    _u(u) = reg.control_type === nothing ? u : reg.control_type(u...)

    ev(y, u, pv) = (dy = zeros(ns);
                    f(dy, reg.state_type(y...), _u(u), pv, p._t0, p.model); dy)

    for (q, jac) in partials_of(f)
        base = q === state     ? y0 :
               q === control   ? u0 :
               q === parameter ? p0 : continue
        isempty(base) && continue
        A = zeros(ns, length(base))
        jac(A, reg.state_type(y0...), _u(u0), p0, p._t0, p.model)
        vary = q === state     ? (y  -> ev(y,  u0, p0)) :
               q === control   ? (u  -> ev(y0, u,  p0)) :
                                 (pv -> ev(y0, u0, pv))
        w = _worst_vs_fd(A, vary, base, step)
        push!(rows, (_fname(f), string(nameof(q)), w))
    end
    return nothing
end

# --- constraints: path and endpoint ------------------------------------------
#
# The user's function and its partial take typed structs, so a perturbation
# happens on the underlying vector and the struct is rebuilt around it.

# The states a constraint is differentiated at. After `initialize!` they are the mesh's own first
# and last columns; before it, they are the guess the phase was given. Probing at zeros instead,
# which is what an uninitialised phase used to give, evaluates a user's function where the physics
# does not hold: an orbit-raising terminal condition reads `sqrt(mu / r)` at r = 0, and the finite
# difference then steps r negative and raises a DomainError inside `check_partials`.
function _probe_states(p, ns)
    size(p._Y, 2) > 0 && return collect(p._Y[1:ns, 1]), collect(p._Y[1:ns, end])
    p.state_var === nothing && return zeros(ns), zeros(ns)
    g = p.state_var.var.data
    return collect(g[1:ns, 1]), collect(g[1:ns, end])
end

function _check_constraint!(rows, p, reg, spec, step)
    f, at = spec.f, spec.at
    ns, nc = reg.n_states, reg.n_controls
    y0, yf = _probe_states(p, ns)
    u0 = p.control_var === nothing ? zeros(nc) : collect(p.control_var.value)[1:nc]
    t0, tf, m = p._t0, p._tf, p.model
    S, C = reg.state_type, reg.control_type

    for (q, jac) in partials_of(f)
        n_args = _arity(jac)

        if n_args == 1                       # context partial
            ctx(y, u) = at isa Path ?
                PhaseContext(S(y...), C(u...), float(t0), float(t0), float(tf), nothing, m) :
                PhaseContext(S(y...), nothing, float(at isa Initial ? t0 : tf),
                             float(t0), float(tf), nothing, m)
            here = at isa Final ? yf : y0
            base, ev = q === state   ? (here, y -> _vec(f(ctx(y, u0)))) :
                       q === control ? (u0,   u -> _vec(f(ctx(here, u)))) :
                       (nothing, nothing)
            if base === nothing && q in (initial_time, final_time)
                # A time partial moves the endpoint the quantity names.
                tsel = q === initial_time ? t0 : tf
                evt(v) = _vec(f(q === initial_time ?
                    PhaseContext(here, nothing, float(v[1]), float(v[1]), float(tf), nothing, m) :
                    PhaseContext(here, nothing, float(v[1]), float(t0), float(v[1]), nothing, m)))
                A = _mat(jac(ctx(here, u0)))
                push!(rows, (_fname(f) * " @" * _atname(at), string(nameof(q)),
                             _worst_vs_fd(A, evt, [tsel], step)))
                continue
            end
            base === nothing && continue
            A = _mat(jac(ctx(here, u0)))
            push!(rows, (_fname(f) * " @" * _atname(at), string(nameof(q)),
                         _worst_vs_fd(A, ev, base, step)))

        elseif n_args == 5                   # path partial: (y, u, p, t, model)
            base, ev = q === state  ? (y0, y -> _vec(f(S(y...), C(u0...), nothing, t0, m))) :
                       q === control ? (u0, u -> _vec(f(S(y0...), C(u...), nothing, t0, m))) :
                       (nothing, nothing)
            base === nothing && continue
            A = _mat(jac(S(y0...), C(u0...), nothing, t0, m))
            push!(rows, (_fname(f) * " @" * _atname(at), string(nameof(q)),
                     _worst_vs_fd(A, ev, base, step)))

        elseif n_args == 6                   # endpoint partial
            bf(a, b, s, e) = _vec(f(S(a...), S(b...), nothing, s, e, m))
            base, ev =
                q === state && at isa Initial ? (y0, y -> bf(y,  yf, t0, tf)) :
                q === state                   ? (yf, y -> bf(y0, y,  t0, tf)) :
                q === initial_time            ? ([t0], t -> bf(y0, yf, t[1], tf)) :
                q === final_time              ? ([tf], t -> bf(y0, yf, t0, t[1])) :
                (nothing, nothing)
            base === nothing && continue
            A = _mat(jac(S(y0...), S(yf...), nothing, t0, tf, m))
            push!(rows, (_fname(f) * " @" * _atname(at), string(nameof(q)),
                     _worst_vs_fd(A, ev, base, step)))
        end
    end
    return nothing
end

# --- partials written against a context --------------------------------------
#
# MGAnDSMs and Sims-Flanagan hand a constraint the phase's boundary context, so
# their partials take one argument. To difference one, rebuild the context with
# a single quantity moved. The map from quantity to field is the same one the
# readers above express in the other direction; it is short and explicit
# because guessing it from the reader would be worse.

const _CTX_FIELD = Dict{Any,Symbol}(
    state          => :x0,       final_state    => :xf,
    departure_vinf => :vinf_dep, arrival_vinf   => :vinf_arr,
    initial_time   => :t0,       final_time     => :tf,
    initial_mass   => :m0,       final_mass     => :mf,
    deep_space_dv  => :dv,       arc_fractions  => :alpha)

_ctx_of(p::MGAnDSMsPhase)     = _mga_boundary_context(p)
_ctx_of(p::SimsFlanaganPhase) = _sf_boundary_context(p)

"Rebuild an immutable context with one field replaced."
_ctx_with(c, field::Symbol, v) =
    (typeof(c).name.wrapper)((f === field ? v : getfield(c, f)
                              for f in fieldnames(typeof(c)))...)

"Put a perturbed flat vector back into the shape the context field had."
_reshape_like(orig::AbstractMatrix, v) = reshape(v, size(orig))
_reshape_like(orig::AbstractVector, v) = v
_reshape_like(orig::Number,         v) = v[1]

function _check_ctx_spec!(rows, p, spec, step)
    f, at = spec.f, spec.at
    for (q, jac) in partials_of(f)
        _arity(jac) == 1 || continue
        field = get(_CTX_FIELD, q, nothing)
        field === nothing && continue

        c0   = _ctx_of(p)
        orig = getfield(c0, field)
        base = _vec(orig)
        # A scalar-valued function (a Mayer objective) declares its derivative
        # as a gradient vector; a vector-valued one declares a Jacobian matrix.
        A    = f(c0) isa Number ? reshape(_vec(jac(c0)), 1, :) : _mat(jac(c0))
        ev(v) = _vec(f(_ctx_with(c0, field, _reshape_like(orig, v))))
        push!(rows, (_fname(f) * " @" * _atname(at), string(nameof(q)),
                     _worst_vs_fd(A, ev, base, step)))
    end
    return nothing
end

"""
    check_partials(phase::Union{MGAnDSMsPhase,SimsFlanaganPhase}; step = 1e-6)

The shooting-side counterpart: every partial declared against a context, on
every constraint and objective attached to the phase, against central
differences of the same function.
"""
function check_partials(p::Union{MGAnDSMsPhase,SimsFlanaganPhase};
                        step = 1e-6, verbose = true)
    rows = Tuple{String,String,Float64}[]
    bare = String[]
    for spec in phase_specs(p)
        isempty(partials_of(spec.f)) && push!(bare, _fname(spec.f))
        _check_ctx_spec!(rows, p, spec, step)
    end
    if verbose
        isempty(rows) && isempty(bare) &&
            println("  check_partials: nothing declared on this phase.")
        for (fname, qname, w) in rows
            @printf("  %-30s worst |analytic - fd| = %.3e
", fname * " wrt " * qname, w)
        end
    end
    # 078d0f3 and 1cfbdd4 replaced the zeros a shooting phase used to return for
    # an unregistered variable block. It differentiates now, so an undeclared
    # derivative costs time rather than being silently absent, and the message
    # says which functions are paying it.
    isempty(bare) || println("  fell back to AD: ", join(unique(bare), ", "))
    worst = isempty(rows) ? 0.0 : maximum(r[3] for r in rows)
    worst > 1e-4 && @warn "check_partials: a declared partial disagrees with " *
                          "finite differences by $(round(worst, sigdigits = 3))."
    return worst
end

_atname(::Path)    = "Path"
_atname(::Initial) = "Initial"
_atname(::Final)   = "Final"
_atname(::Any)     = "Boundary"

"A partial may return a matrix, a row, or a scalar; compare them all as matrices."
_mat(A::AbstractMatrix) = A
_mat(A::AbstractVector) = reshape(collect(float.(A)), :, 1)
_mat(A::Number)         = fill(float(A), 1, 1)


# =============================================================================
# MGAnDSMs: the same words, on quantities that are not a state vector
# =============================================================================

# A quantity read off the context a sequence constraint is handed. The same
# name that declares the variable also reads its value — `final_time(phase)`
# in a Vary, `final_time(ctx)` in a constraint.
# What a Link between two collocation phases is handed. Nothing read these
# before, because a Link always took the shooting route, so a condition written
# over two collocation phases had no vocabulary at all.
initial_state(c::BoundaryContext) = c.y0
final_state(c::BoundaryContext)   = c.yf
initial_time(c::BoundaryContext)  = c.t0
final_time(c::BoundaryContext)    = c.tf
# A link in a sequence of mixed phases is handed the unified context, which reports a shooting
# phase's mass as the last state component.
initial_state(c::UnifiedBoundaryCtx) = c.y0
final_state(c::UnifiedBoundaryCtx)   = c.yf
initial_time(c::UnifiedBoundaryCtx)  = c.t0
final_time(c::UnifiedBoundaryCtx)    = c.tf

initial_time(c::MGABoundaryContext)   = c.t0
final_time(c::MGABoundaryContext)     = c.tf
initial_mass(c::MGABoundaryContext)   = c.m0
final_mass(c::MGABoundaryContext)     = c.mf
departure_vinf(c::MGABoundaryContext) = c.vinf_dep
arrival_vinf(c::MGABoundaryContext)   = c.vinf_arr
deep_space_dv(c::MGABoundaryContext)  = c.dv
arc_fractions(c::MGABoundaryContext)  = c.alpha

"""
    slot_of(phase, quantity) -> the phase's variable for that quantity

Both phase families keep their decision variables in named slots, so a
quantity is a name for a slot and one lookup serves collocation and MGAnDSMs
alike. `initial_time` resolves to `t0_var` on either — the vocabulary was
common before anything was unified.
"""
slot_of(p, ::typeof(state))          = p.state_var
slot_of(p, ::typeof(control))        = p.control_var
slot_of(p, ::typeof(initial_time))   = p.t0_var
slot_of(p, ::typeof(final_time))     = p.tf_var
slot_of(p, ::typeof(departure_vinf)) = p.vinf_dep_var
slot_of(p, ::typeof(arrival_vinf))   = p.vinf_arr_var
slot_of(p, ::typeof(initial_mass))   = p.m0_var
slot_of(p, ::typeof(final_mass))     = p.mf_var
slot_of(p, ::typeof(deep_space_dv))  = p.dv_var
slot_of(p, ::typeof(arc_fractions))  = p.alpha_var
slot_of(p, ::typeof(forward_control))  = p.u_fwd_var
slot_of(p, ::typeof(backward_control)) = p.u_bwd_var

# Sims-Flanagan. Its variables are not a state and a control — forward and
# backward throttle blocks meeting at a match point — which is exactly why they
# get their own names rather than being forced into `state` and `control`. The
# verb does not change; only what it is applied to.
# A free endpoint state. With no ephemeris the transcription is just two half
# propagations meeting at a match point, so its ends can be decision variables
# like any other — which is what lets a Sims-Flanagan phase sit in a sequence
# beside a collocation phase instead of only between planets.
Vary(::typeof(state), p::SimsFlanaganPhase; guess, lower_bound, upper_bound,
     scale = nothing, shift = nothing, name = "x0") =
    (v = _dsv(SFStateBlock(), lower_bound, upper_bound, scale,
              isempty(name) ? "x0" : name; shift = shift, guess = guess);
     p._x0 .= _vec(guess); set_initial_state!(p, v); v)

Vary(::typeof(final_state), p::SimsFlanaganPhase; guess, lower_bound, upper_bound,
     scale = nothing, shift = nothing, name = "xf") =
    (v = _dsv(SFStateBlock(), lower_bound, upper_bound, scale,
              isempty(name) ? "xf" : name; shift = shift, guess = guess);
     p._xf .= _vec(guess); set_final_state!(p, v); v)

# The phase holds the values the solver starts from, so a guess is written there. A quantity
# pinned by equal bounds starts at that value; one with neither keeps what the phase has, which
# for the times is `tspan`.
for (q, tag, setter, field) in ((:forward_control,  :SFThrustBlock, :set_fwd_control!,  :_u_fwd),
                                (:backward_control, :SFThrustBlock, :set_bwd_control!,  :_u_bwd),
                                (:initial_mass,     :SFMassParam,   :set_initial_mass!, :_m0),
                                (:final_mass,       :SFMassParam,   :set_final_mass!,   :_mf),
                                (:initial_time,     :SFTime,        :set_initial_time!, :_t0),
                                (:final_time,       :SFTime,        :set_final_time!,   :_tf))
    @eval function Vary(::typeof($q), p::SimsFlanaganPhase;
                        lower_bound, upper_bound, scale = nothing,
                        shift = nothing, guess = nothing,
                        name = string($(QuoteNode(q))))
        v = _dsv($tag(), lower_bound, upper_bound, scale,
                 isempty(name) ? string($(QuoteNode(q))) : name;
                 shift = shift, guess = guess)
        start = guess !== nothing ? _vec(guess) :
                _vec(lower_bound) == _vec(upper_bound) ? _vec(lower_bound) : nothing
        start === nothing || _sf_start!(p, Val($(QuoteNode(field))), start)
        # The variable's value mirrors the phase from the start, as set_decision_vector! keeps
        # it afterwards, so a derivative taken before the first solver call starts from the
        # full throttle block rather than from one 3-vector.
        f = getfield(p, $(QuoteNode(field)))
        v.value = f isa AbstractArray ? vec(copy(f)) : [f]
        $setter(p, v)
        return v
    end
end

# A throttle guess is one direction for every segment, or all 3 × n segments of the half.
function _sf_start!(p::SimsFlanaganPhase, ::Val{F}, g) where {F}
    if F in (:_u_fwd, :_u_bwd)
        n = n_fwd(p)
        length(g) in (3, 3n) || throw(ArgumentError(
            "Vary: a throttle guess is one 3-vector for every segment or a 3 × $n matrix; " *
            "got $(length(g)) values"))
        getfield(p, F) .= reshape(length(g) == 3 ? repeat(g, n) : g, 3, n)
    else
        setfield!(p, F, float(first(g)))
    end
    return nothing
end

# The starting value defaults to the lower bound, which is what a shooting
# variable did before `guess` reached this side of the vocabulary.
_dsv(tag, lb, ub, scale, name; shift = nothing, guess = nothing) =
    (n = length(_vec(lb));
     DirectSolverVariable(tag, _vec(lb), _vec(ub),
                          _scale_vec(scale, n), _shift_vec(shift, n),
                          name, guess === nothing ? _vec(lb) : _vec(guess)))

# The guess rides on the Vary, next to the bounds it lives between, rather
# than in a separate block naming every quantity a second time.
for (q, tag, setter, field) in
        ((:departure_vinf, :MGAVInfinity3, :set_departure_vinf!, :vinf_dep),
         (:arrival_vinf,   :MGAVInfinity3, :set_arrival_vinf!,   :vinf_arr),
         (:initial_mass,   :MGAMassParam,  :set_initial_mass!,   :m0),
         (:final_mass,     :MGAMassParam,  :set_final_mass!,     :mf),
         (:initial_time,   :MGATime,       :set_initial_time!,   :t0),
         (:final_time,     :MGATime,       :set_final_time!,     :tf))
    @eval function Vary(::typeof($q), p::MGAnDSMsPhase;
                        lower_bound, upper_bound, scale = nothing,
                        shift = nothing, guess = nothing,
                        name = string($(QuoteNode(q))))
        v = _dsv($tag(), lower_bound, upper_bound, scale,
                 isempty(name) ? string($(QuoteNode(q))) : name; shift = shift)
        if guess !== nothing
            v.value .= _vec(guess)
            set_initial_guess!(p; $field = length(_vec(guess)) == 1 ?
                                            first(_vec(guess)) : _vec(guess))
        end
        $setter(p, v)
        return v
    end
end

# The burns and the arc fractions. Every piece below the vocabulary already handled them: the
# setters, the per-DSM bounds and scaling, the match-point Jacobian blocks and the alpha-sum row.
# What was missing was the verb, so the burns a transcription is named for could not be handed to
# a solver.
#
# A burn's bounds, scale and shift are one DSM's three components and apply to every DSM, as the
# NLP layout already tiles them. The guess may be one burn, repeated, or all 3 × n_dsm values.
function Vary(::typeof(deep_space_dv), p::MGAnDSMsPhase;
              lower_bound, upper_bound, scale = nothing, shift = nothing, guess = nothing,
              name = "deep_space_dv")
    n = n_dsm(p)
    (length(_vec(lower_bound)) == 3 && length(_vec(upper_bound)) == 3) || throw(ArgumentError(
        "Vary(deep_space_dv, phase): bounds are one burn's three components, applied to each of " *
        "the $n DSMs; got $(length(_vec(lower_bound))) and $(length(_vec(upper_bound)))"))
    g = guess === nothing ? vec(p._dv) : _vec(guess)
    length(g) == 3 && (g = repeat(g, n))
    length(g) == 3n || throw(ArgumentError(
        "Vary(deep_space_dv, phase): guess must be one burn (3 values) or every burn " *
        "(3 × $n = $(3n) values); got $(length(g))"))
    v = _dsv(MGADVBlock(), lower_bound, upper_bound, scale, isempty(name) ? "deep_space_dv" : name;
             shift = shift, guess = g)
    set_dsm_control!(p, v)
    set_initial_guess!(p; dv = reshape(g, 3, n))
    return v
end

# One fraction per arc. Varying them adds the equality that they sum to one.
function Vary(::typeof(arc_fractions), p::MGAnDSMsPhase;
              lower_bound, upper_bound, scale = nothing, shift = nothing, guess = nothing,
              name = "arc_fractions")
    na = n_alphas(p)
    (length(_vec(lower_bound)) == na && length(_vec(upper_bound)) == na) || throw(ArgumentError(
        "Vary(arc_fractions, phase): bounds need one value per arc, $na; got " *
        "$(length(_vec(lower_bound))) and $(length(_vec(upper_bound)))"))
    g = guess === nothing ? copy(p._alpha) : _vec(guess)
    length(g) == na || throw(ArgumentError(
        "Vary(arc_fractions, phase): guess needs one value per arc, $na; got $(length(g))"))
    v = _dsv(MGAAlphaBlock(), lower_bound, upper_bound, scale, isempty(name) ? "arc_fractions" : name;
             shift = shift, guess = g)
    set_alpha!(p, v)
    set_initial_guess!(p; alpha = g)
    return v
end

# =============================================================================
# Link — the thing a condition between phases is about
# =============================================================================

"""
    Link(phases...; body = nothing)

What ties one phase to the next: a flyby, a match, a handover.

A `Link` is a subject, so `Constraint` keeps taking exactly one and the four
conditions at a flyby are four constraints on the same thing rather than four
repetitions of a pair. It carries what the link is — `body` is what the
periapsis condition measures against.

Which side a quantity means is in its own name: `arrival_vinf`, `final_time`
and `final_mass` can only be the phase that is ending, `departure_vinf`,
`initial_time` and `initial_mass` only the one starting.

Two collocation phases that meet, state and time, are joined with
`Constraint(continuity, link)`. See [`continuity`](@ref).

# Fields
- `phases::Tuple`: the phases the link joins, in the order they are flown.
- `body`: the body a flyby is about, which is what the periapsis condition measures against, or
  `nothing` for a link that is not a flyby.
- `name::Symbol`: the link's name, which labels it in reports.

# Example
<!-- doc-fragment -->
```julia
link = Link(phase_in, phase_out; body = venus, name = :venus_flyby)

Constraint(continuity, link)
```
"""
struct Link{P<:Tuple}
    phases ::P
    body   ::Any
    name   ::Symbol
end

Link(phases...; body = nothing, name = :link) = Link(phases, body, name)

# Constraints are collected here and registered when the Sequence is
# initialised, because a Link is written before the Sequence exists.
const _PENDING = Vector{Any}()

function Constraint(f::Function, l::Link; equals = nothing, lower_bound = nothing,
                    upper_bound = nothing, scale = nothing,
                    name = gensym(:link), at = nothing)
    # A link joins two phases where they meet. There is no other point it could
    # be evaluated at, so `at` is accepted and refused rather than missing.
    at === nothing || throw(ArgumentError(
        "Constraint: `at` has no meaning on a Link — a link is evaluated where " *
        "the two phases meet, which is the only point it has."))
    if equals !== nothing
        lb = ub = _bounds(equals)
    else
        n  = length(_bounds(something(lower_bound, upper_bound)))
        lb = lower_bound === nothing ? fill(-Inf, n) : _bounds(lower_bound)
        ub = upper_bound === nothing ? fill( Inf, n) : _bounds(upper_bound)
    end
    push!(_PENDING, (l, f, lb, ub, String(name)))
    return nothing
end

"""
    continuity

The condition that two linked phases meet: the state where the first ends equals the state where
the second begins, and so does the time. The phases can be collocation, shooting, or one of each.

# Example
```julia
halfway = Link(first_half, second_half)
Constraint(continuity, halfway)
```

# Notes
Written as `Constraint(continuity, link)`, with no bounds, since the residual is zero by
definition. A shooting phase keeps mass beside its state, and continuity counts it as the last
state component, so a Sims-Flanagan phase with six states joins a collocation phase with seven,
the seventh being mass. Throws `ArgumentError` when the two phases carry different numbers of
components. A link that joins only some components, or joins them with an offset such as a mass
dropped at staging, is written as a `Constraint` with a function of its own, which reads the two
ends with `final_state`, `initial_state`, `final_time` and `initial_time`.
"""
function continuity(c1, c2)
    y1, y2 = final_state(c1), initial_state(c2)
    n = y1 isa AbstractVector ? length(y1) : fieldcount(typeof(y1))
    return vcat([_flat(y1, n)[i] - _flat(y2, n)[i] for i in 1:n],
                final_time(c1) - initial_time(c2))
end

function Constraint(::typeof(continuity), l::Link; name = :continuity)
    length(l.phases) == 2 || throw(ArgumentError(
        "Constraint(continuity, link): continuity joins two phases; the link has " *
        "$(length(l.phases))."))
    push!(_PENDING, (l, continuity, nothing, nothing, String(name)))
    return nothing
end

function _continuity_size(p1, p2)
    n1, n2 = p1._n_states, p2._n_states
    (n1 > 0 && n1 == n2) || throw(ArgumentError(
        "Constraint(continuity, link): phases $(p1.name) and $(p2.name) must have states of " *
        "the same length; got $n1 and $n2 components. Join phases with different states with " *
        "a function of your own."))
    return n1
end

"""Register every Link constraint written so far against `seq`.

Two phases meeting is one idea with two mechanisms underneath. Shooting phases
join through `add_sequence_constraint!` and the shooting manager; collocation
phases join through `add_linkage!` and the collocation manager. A `Link` writes
the same way over either, and this is where that is made true — until now it
always took the shooting route, so continuity between two collocation phases
could not be written in this vocabulary at all.
"""
function drain_links!(seq)
    # A sequence holding both collocation and shooting phases is assembled by the unified manager,
    # and every link in it goes there. Its boundary context reports each phase's state, with a
    # shooting phase's mass appended, as one vector, so continuity is one equality over it.
    mixed = !isempty(sequence_phases(seq)) && !isempty(sequence_sf_phases(seq))
    for (l, f, lb, ub, nm) in _PENDING
        p1, p2 = l.phases[1], l.phases[2]
        both_coll = p1 isa CollocationPhase && p2 isa CollocationPhase
        if mixed || (f === continuity && !both_coll)
            f === continuity ? add_continuity!(seq, p1, p2) :
                add_continuity!((c1, c2) -> f(c1, c2), seq, p1, p2;
                                lower_bounds = lb, upper_bounds = ub, name = nm)
            continue
        end
        if f === continuity
            # The residual length is only known once both phases have their state declared, which
            # a script may do after writing the link, so the bounds are sized here.
            n = _continuity_size(p1, p2)
            lb = ub = zeros(n + 1)
        end
        if p1 isa CollocationPhase && p2 isa CollocationPhase
            add_linkage!((c1, c2) -> f(c1, c2), seq, p1, p2;
                         lower_bounds = lb, upper_bounds = ub, name = nm)
        else
            add_sequence_constraint!((c1, c2) -> f(c1, c2), seq, p1, p2;
                                     lower_bounds = lb, upper_bounds = ub, name = nm)
        end
    end
    empty!(_PENDING)
    return nothing
end

# =============================================================================
# Estimation: the same words, on a spacecraft
#
# `Vary(state, sc)` is what `Vary(delta_v, toi)` is — the quantity names what
# is estimated, the subject names what owns it. SolveFor is what Vary means, so
# the role disappears, and CartesianStateVar wrapped in ModelVariable was
# saying "the state, of that spacecraft" twice.
#
# The estimation stack is reached lazily so this file does not depend on it.
# =============================================================================

# The estimators are AstroSolve's submodules now, not modules in Main.
_est(m) = getfield(AstroSolve, m)

# Vary(state, spacecraft) is AstroSolve's own now. It reaches the state through
# the `state` quantity like any other Vary, and it rejects a box bound given
# alongside a covariance, which is what this method used to exist for.

# =============================================================================
# One solve verb
#
# Seven verbs said the same thing and differed only in which subsystem you had
# reached for: solve!, solve_trajectory!, solve_trajectory_oc!,
# solve_batch_ls!, run_iterated_rts!, and init_ekf plus a loop. Which one you
# needed depended on the transcription — the exact thing the vocabulary exists
# to stop mattering.
#
#     solve!(seq;               method = Optimize(max_iter = 500))
#     solve!(problem, tracking; method = Batch(n_iters = 15))
#     solve!(problem, tracking; method = Sequential(smoother = RTS()))
#
# The method says how. The subject carries its own variables. The names avoid
# IPOPT, BatchLeastSquares and ExtendedKalmanFilter, which are already taken by
# SNOW and by the estimation modules — this codebase has been bitten twice by
# two subsystems claiming one name.
# =============================================================================

"""
    Optimize(; max_iter = 500, tol = 1e-6, print_level = 5, derivatives = :user,
               record_iterations = false, extra = Dict{String,Any}())

`Optimize` says how an NLP is solved, and is passed as `solve!(seq; method = Optimize(...))`.

# Fields
- `max_iter::Int`: the iteration limit. A solve that reaches it returns
  `Maximum_Iterations_Exceeded` rather than a converged answer, so the status is worth reading.
- `tol::Float64`: the tolerance the solver's own error measure must fall below before it reports
  success.
- `print_level::Int`: how much the solver prints, from 0 for silence.
- `derivatives::Symbol`: `:user` to use the partials declared with `@partial`, or `:fd` for finite
  differences throughout.
- `record_iterations::Bool`: keeps the iterate history where the solver supports it.
- `extra::Dict{String,Any}`: solver-specific options passed straight through, which is the escape
  hatch for what only one solver has rather than a keyword per knob.

# Example
```julia
Optimize(max_iter = 500, print_level = 0)
```
"""
struct Optimize
    max_iter           ::Int
    tol                ::Float64
    print_level        ::Int
    derivatives        ::Symbol
    record_iterations  ::Bool
    extra              ::Dict{String,Any}
end
Optimize(; max_iter = 500, tol = 1e-6, print_level = 5, derivatives = :user,
           record_iterations = false, extra = Dict{String,Any}()) =
    Optimize(max_iter, tol, print_level, derivatives, record_iterations,
             Dict{String,Any}(extra))

"""
    Batch(; n_iters = 15, tol = 1e-9, verbose = false)

Batch least squares: every observation is processed at once and the fit is iterated to convergence.
Passed as `solve!(problem, records; method = Batch(...))`.

# Fields
- `n_iters::Int`: the maximum number of least-squares iterations.
- `tol::Float64`: the convergence tolerance on the iteration.
- `verbose::Bool`: prints per-iteration progress.

# Example
```julia
Batch(n_iters = 10, tol = 1e-9)
```
"""
struct Batch
    n_iters ::Int
    tol     ::Float64
    verbose ::Bool
end
Batch(; n_iters = 15, tol = 1e-9, verbose = false) = Batch(n_iters, tol, verbose)

"""
    RTS()

The Rauch-Tung-Striebel smoother, a backward pass over an arc a filter has already run forward.
Passed to [`Sequential`](@ref) as `smoother = RTS()`, and it revisits each epoch with everything the
filter learned afterwards, so the smoothed estimate at an early epoch uses later data.

# Example
<!-- doc-fragment -->
```julia
Sequential(iterations = 1, smoother = RTS())
```
"""
struct RTS end

"""
    Sequential(; iterations = 5, tol = 1e-6, step_size = 1.0, smoother = nothing)

An extended Kalman filter: observations are processed one at a time, in order. Passed as
`solve!(problem, records; method = Sequential(...))`.

# Fields
- `iterations::Int`: how many times the filter is re-run over the arc, each pass starting from the
  last pass's estimate.
- `tol::Float64`: the convergence tolerance across those passes.
- `step_size::Float64`: scales the state update applied at each measurement.
- `smoother`: [`RTS`](@ref) to add a backward smoothing pass, or `nothing` for the filter alone.

# Example
```julia
Sequential(iterations = 1, smoother = RTS())
```
"""
struct Sequential
    iterations ::Int
    tol        ::Float64
    step_size  ::Float64
    smoother   ::Any
end
Sequential(; iterations = 5, tol = 1e-6, step_size = 1.0, smoother = nothing) =
    Sequential(iterations, tol, step_size, smoother)

"""A log path no other solve in this session will ask for.

IPOPT opens its output file and does not release it, so a second solve that
names the same path fails with "Couldn't open output file". Naming none is the
same thing by another route: IPOPT falls back to `ipopt.out` in the working
directory, which is why solving twice in one REPL used to die on the second and
why the use case folder filled with empty logs.

The timestamp is for a human reading a directory listing; `time_ns()` is what
makes it unique. It goes in the temp directory so a run leaves the repo alone,
and `extra` overrides it for anyone who wants the log kept somewhere.
"""
_ipopt_log() = joinpath(tempdir(),
    "ipopt_" * Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS") *
    "_" * string(time_ns()) * ".out")

function _snow(m::Optimize)
    opts = Dict{String,Any}("max_iter"    => m.max_iter,
                            "tol"         => m.tol,
                            "print_level" => m.print_level,
                            "output_file" => _ipopt_log())
    merge!(opts, m.extra)
    return SNOW.Options(
        derivatives = m.derivatives === :fd ? SNOW.ForwardFD() : SNOW.UserDeriv(),
        solver = SNOW.IPOPT(opts))
end

"""
    solve!(seq; method = Optimize())

Solve a sequence, whatever is in it. A sequence carrying continuity links goes
through the unified manager; one without goes through the event-graph path,
which itself routes to the shooting solver when it holds shooting phases. The
caller does not choose.

# Returns
The solver's result. `result.info` is the solver's exit status, and is worth reading rather than
assuming: a run that stops on the iteration limit reports `Maximum_Iterations_Exceeded` and still
returns the last iterate, which may be neither optimal nor feasible. `result.objective` is the
achieved objective value. The solved values themselves are read off the phase, with
[`get_final_state`](@ref) and the quantity accessors.

# Example
```julia
result = solve!(Sequence(phase); method = Optimize(max_iter = 500, print_level = 0))
```
"""
# =============================================================================
# check_configuration — a misconfigured element errors
#
# Neither of the two mistakes below announces itself. Both present as a hard
# problem, which is the worst possible disguise, and both have cost days.
#
#   * On a shooting phase an undeclared derivative is a zero, not a finite
#     difference: objective_gradient_chunk and jacobian_chunk return
#     zeros(nlp_length(p, var)) for a variable block with no registered
#     Jacobian. uc7 spent months maximising delivered mass against an
#     identically zero gradient and reporting Maximum_Iterations_Exceeded.
#
#   * A forgotten `Vary` silently shrinks the NLP. Dropping one from uc7 takes
#     it from 20 decision variables to 19, no error, and prints 3.115 kg where
#     the intact case gives 91.549.
#
# Something may legitimately want a zero derivative. Say so —
# `@partial(f, q) do c; zeros(...); end` — rather than leaving it to silence.
# =============================================================================

"""
    check_configuration(seq)

Reject a sequence that cannot mean what was written. Called by `solve!`.
"""
function check_configuration(seq::Sequence)
    problems = String[]
    for p in sequence_sf_phases(seq)
        _check_shooting_config!(problems, p)
    end
    isempty(problems) && return nothing
    throw(ArgumentError("this problem is not fully specified:
  " *
                        join(problems, "
  ") *
                        "

See SPEC section 8. A derivative that is " *
                        "genuinely zero should say so with `@partial`."))
end

"""
Every constraint and objective on a shooting phase needs its derivatives.

Ask the framework what it holds, not what the user typed. A path constraint
declares no `@partial` and is still fully differentiated — the framework builds
its block-diagonal Jacobian with ForwardDiff — so asking about `@partial` would
reject correct code. The question is whether a Jacobian is registered.
"""
# A phase that differentiates its constraints with ForwardDiff is exempt: an
# undeclared derivative there is computed, not zeroed. The trap belongs to the
# phases that return zeros(nlp_length(p, var)) instead.
_check_shooting_config!(problems, p) = nothing

"""
A declared partial names the quantity it differentiates. If that quantity is
not a variable on this phase, the derivative was registered against nothing —
`variable_list` builds itself with `filter(!isnothing, ...)`, so an empty slot
does not error, it disappears, and the NLP is quietly one variable smaller.
Dropping a single `Vary` from uc7 took it from 20 decision variables to 19 and
printed 3.115 kg where the intact case gives 91.549.
"""
function _check_dead_partials!(problems, p)
    for spec in phase_specs(p), (q, _) in partials_of(spec.f)
        slot_of(p, q) === nothing || continue
        push!(problems, "`$(_fname(spec.f))` on phase $(p.name) declares a " *
              "partial with respect to `$(nameof(q))`, which is not a variable " *
              "on that phase — add `Vary($(nameof(q)), phase; ...)`, or drop " *
              "the partial. As written the derivative is registered against " *
              "nothing.")
    end
    return nothing
end

"""
Check that a registered derivative is the shape the framework will use.

A Mayer objective's derivative is a gradient vector; a boundary constraint's is
a Jacobian matrix — for the identical expression. Handing the objective a 1xN
matrix broadcasts it against the N-element scale vector into an NxN and takes
the process down inside Ipopt with no stacktrace at all.
"""
function _check_jacobian_shapes!(problems, p)
    byid = Dict(objectid(v) => v for v in variable_list(p))
    say(name, got, want) = push!(problems,
        "`$name` on phase $(p.name) declares a derivative of size $got where " *
        "the framework will use $want.")

    for c in _constraints_of(p)
        c isa BoundaryConstraint || continue
        rows = length(c.lower_bounds)
        for (vid, fn) in registered_jacobians(c.calc)
            v = get(byid, vid, nothing); v === nothing && continue
            J = try fn() catch; continue end
            want = (rows, nlp_length(p, v))
            size(J) == want || say(c.calc.name, size(J), want)
        end
    end

    obj = p.objective
    obj === nothing && return nothing
    for m in (obj isa MayerObjective ? (obj,) : _mayers_of(obj))
        m isa MayerObjective || continue
        for (vid, fn) in registered_jacobians(m)
            v = get(byid, vid, nothing); v === nothing && continue
            g = try fn() catch; continue end
            n = nlp_length(p, v)
            g isa AbstractVector && length(g) == n && continue
            say("the objective", size(g),
                "a gradient vector of length $n (a Mayer term is a vector, " *
                "not a 1x$n matrix)")
        end
    end
    return nothing
end

function _check_shooting_config!(problems,
                                 p::Union{MGAnDSMsPhase,SimsFlanaganPhase})
    _check_dead_partials!(problems, p)
    _check_jacobian_shapes!(problems, p)
    _blank(name, what) = push!(problems,
        "`$name` on phase $(p.name) has no registered $what, so the framework " *
        "differentiates it. That is correct but slower, and it requires the " *
        "function to accept dual numbers: annotate its arguments `Real`, not " *
        "`Float64`.")

    # A constraint written without `name` carries a generated one, `##con#280`, which tells the
    # user nothing. Boundary constraints and the endpoint specs recorded as they were written come
    # in the same order, so the function the user wrote is found by position.
    written = [sp.f for sp in phase_specs(p) if sp.kind === :constraint && !(sp.at isa Path)]
    for (i, c) in enumerate(filter(c -> c isa BoundaryConstraint, collect(_constraints_of(p))))
        isempty(registered_jacobians(c.calc)) || continue
        cname = string(c.calc.name)
        label = startswith(cname, "##") && i <= length(written) ? _fname(written[i]) : cname
        _blank(label, "Jacobian")
    end
    obj = p.objective
    obj === nothing && return nothing
    # The framework holds its own wrapper closure, whose name is unreadable.
    # Report the function the user actually wrote.
    ospec = findfirst(sp -> sp.kind === :objective, collect(phase_specs(p)))
    oname = ospec === nothing ? "the objective" :
            _fname(collect(phase_specs(p))[ospec].f)
    for m in (obj isa MayerObjective ? (obj,) : _mayers_of(obj))
        m isa MayerObjective || continue
        isempty(registered_jacobians(m)) && _blank(oname, "gradient")
    end
    return nothing
end

"Phases spell the field differently, `constraints` or `_constraints`."
_constraints_of(p) = hasproperty(p, :constraints)  ? p.constraints  :
                     hasproperty(p, :_constraints) ? p._constraints : ()

"An objective handle may hold its Mayer terms rather than being one."
_mayers_of(obj) = hasproperty(obj, :mayers) ? obj.mayers :
                  hasproperty(obj, :mayer)  ? (obj.mayer,) : ()

"""
    initialize_sequence!(seq)

Prepare a sequence for solving, whichever kind it is, and reject a
configuration that cannot mean what was written.

There are three initialise verbs — `initialize!`, `initialize_oc!`,
`initialize_shooting!` — and which one applied depended on what the sequence
held, so seven scripts picked by hand. `solve!` calls this, so no ordinary path
picks anything. The underlying three stay reachable: uc6_jacobian_check builds
the NLP to check a Jacobian before solving, and uc12 asks the shooting manager
for its size. This is deliberately not named `initialize!` — that name is
already taken in Main by the collocation one, and a second definition with the
same positional signature would silently replace it, which is a mistake this
codebase has made twice.
"""
function initialize_sequence!(seq::Sequence)
    # uc1/uc2/uc3 build an event graph and hold no phases at all; there is
    # nothing to initialise and `initialize!` would demand a CollocationPhase.
    if !isempty(sequence_phases(seq)) || !isempty(sequence_sf_phases(seq))
        # Register the Link constraints first: a link between a collocation and a shooting phase
        # is what makes the sequence mixed, and that decides how it is initialised.
        drain_links!(seq)
        isempty(sequence_oc_links(seq)) ? initialize!(seq) : initialize_oc!(seq)
    end
    check_configuration(seq)
    return nothing
end

function solve!(seq::Sequence; method::Optimize = Optimize())
    initialize_sequence!(seq)
    if isempty(sequence_oc_links(seq))
        # Two functions are called solve_trajectory!. AstroSolve's drives an
        # event graph and takes record_iterations; OptControlStubs defines its
        # own at top level, which shadows it in Main, and drives phases. They
        # are for different things, so pick by what the sequence holds.
        #
        # Both phase registries have to be asked. Collocation phases and
        # shooting phases (MGAnDSMs, Sims-Flanagan) are kept in separate
        # side-channel dicts, and a sequence of shooting phases looked empty
        # when only the first was checked — it went down the event-graph path
        # and died in Ipopt on untyped bounds.
        if isempty(sequence_phases(seq)) && isempty(sequence_sf_phases(seq))
            return AstroSolve.solve_trajectory!(seq, _snow(method);
                       record_iterations = method.record_iterations)
        end
        return solve_trajectory!(seq, _snow(method))
    end
    return solve_trajectory_oc!(seq, _snow(method))
end

"""
    solve!(problem, data; method = Batch())

Estimate, from a tracking file or already-parsed records. The problem carries
what is being solved for, so the call says only how.

`data` is either a [`TrackingDataFile`](@ref), which is read here, or a vector of
[`ObservationRecord`](@ref) already in hand.

# Returns
The fit. `fit.X_hat` is the estimated state and `fit.sigma` its formal standard deviations, and a
`Sequential` fit also carries `fit.ekf.records`, one per observation, holding the prefit and postfit
residuals in the order the filter processed them.

Throws an `ArgumentError` when the problem's `solve_for` list is empty, since there is then nothing
to estimate.

# Example
```julia
fit = solve!(problem, records; method = Batch(n_iters = 10, tol = 1e-9))
```
"""
function solve!(problem, data; method = Batch())
    BLS = _est(:BatchLeastSquares)
    EKF = _est(:ExtendedKalmanFilter)
    # solve_for is Vector{Any} so a problem can hold variables of any
    # kind; the estimators dispatch on element type, so narrow it here.
    vars = identity.(problem.solve_for)
    isempty(vars) && throw(ArgumentError(
        "solve!: nothing to estimate — give the problem `solve_for = [Vary(...)]`."))
    # A tracking data file is read here, so a script can pass the file it was given.
    TDIO = _est(:TrackingDataIO)
    data isa TDIO.TrackingDataFile && (data = first(TDIO.read_records(data)))

    if method isa Batch
        return BLS.solve_batch_ls!(vars, data; model = problem,
                                   n_iters = method.n_iters, tol = method.tol,
                                   verbose = method.verbose)
    elseif method isa Sequential
        method.smoother === nothing &&
            return EKF.run_ekf!(vars, data; model = problem)
        return EKF.run_iterated_rts!(vars, data; model = problem,
                                     max_iters = method.iterations,
                                     tol = method.tol,
                                     step_size = method.step_size,
                                     verbose = false)
    end
    throw(ArgumentError("solve!: unknown method $(typeof(method))"))
end

push!(_INIT_HOOKS, drain_links!)
