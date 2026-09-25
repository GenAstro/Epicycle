# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A targeting problem written as a script, in the order the mission flies.
#
# GMAT users write a Target block: propagate, vary a burn, apply it, propagate, achieve a goal.
# `target!` reads the same way. The block runs once with a recorder installed, and in that mode the
# verbs do not act: `propagate!` and `maneuver!` record a deferred action, `Vary` attaches to the
# action that follows it, and `Constraint` attaches to the action before it. The recording becomes
# an event sequence, chained in the order written, and is solved the way an assembled `Sequence`
# is. Outside a block the same verbs act immediately.
#
# The hooks for `propagate!` and `maneuver!` live in AstroProp and AstroManeuvers, since those are
# the functions a script calls. `Vary` and `Constraint` are AstroSolve's own and register directly.

import AstroProp, AstroManeuvers

mutable struct _TargetBlock
    events::Vector{Event}
    pending_vars::Vector{Any}
end

const _BLOCK = Ref{Union{Nothing, _TargetBlock}}(nothing)

function _record_action(name::AbstractString, action::Function)
    b = _BLOCK[]
    push!(b.events, Event(name = String(name), event = action,
                          vars = copy(b.pending_vars), funcs = Any[]))
    empty!(b.pending_vars)
    return nothing
end

"""Attach a variable to the next recorded action, when a block is recording."""
function _register_vary!(sv)
    b = _BLOCK[]
    b === nothing || push!(b.pending_vars, sv)
    return sv
end

"""Attach a constraint to the last recorded action, when a block is recording."""
function _register_constraint!(con)
    b = _BLOCK[]
    b === nothing && return con
    isempty(b.events) && throw(ArgumentError(
        "target!: a Constraint must follow a propagate! or maneuver! in the block, since it " *
        "checks what that step produced; got a Constraint before any step"))
    push!(b.events[end].funcs, con)
    return con
end

"""
    target!(block; method = Optimize(derivatives = :fd))

Solve a targeting problem written as a script, in the order the mission flies.

Inside `block`, `propagate!` and `maneuver!` are recorded rather than run. A `Vary` applies to the
step written after it, and a `Constraint` checks the step written before it. The recorded steps
are chained in the order written and solved as one sequence, so a block reads like a GMAT Target
sequence.

# Arguments
- `block`: A function of no arguments, usually written as a `do` block.
- `method`: How the sequence is solved. See [`Optimize`](@ref). Finite differences are the
  default, since the steps in a block carry no declared partials.

# Notes
Throws `ArgumentError` if blocks are nested, if the block records no step, if a `Constraint` comes
before any step, or if a `Vary` is not followed by a step for it to apply to.

# Returns
The solver result, as `solve!` returns it. The varied quantities hold their solved values, and
the spacecraft is left where the last solver pass put it.

# Example
<!-- doc-fragment -->
```julia
result = target!() do
    Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [2.5, 0.0, 0.0])
    maneuver!(sat, toi)

    propagate!(prop, sat, StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1))

    Vary(delta_v, moi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [3.0, 0.0, 0.0])
    maneuver!(sat, moi)
    Constraint(position_magnitude, sat; equals = 45000.0)
    Constraint(eccentricity, sat; equals = 0.0)
end
```
"""
function target!(block::Function; method::Optimize = Optimize(derivatives = :fd))
    _BLOCK[] === nothing || throw(ArgumentError(
        "target!: blocks cannot be nested; got a target! inside another target!"))

    b = _TargetBlock(Event[], Any[])
    _BLOCK[] = b
    AstroProp._RECORDER[]      = _record_action
    AstroManeuvers._RECORDER[] = _record_action
    try
        block()
    finally
        _BLOCK[] = nothing
        AstroProp._RECORDER[]      = nothing
        AstroManeuvers._RECORDER[] = nothing
    end

    isempty(b.events) && throw(ArgumentError(
        "target!: the block must contain at least one propagate! or maneuver!; it recorded none"))
    isempty(b.pending_vars) || throw(ArgumentError(
        "target!: every Vary must be followed by a propagate! or maneuver! that it applies to; " *
        "got $(length(b.pending_vars)) Vary after the last step"))

    seq = Sequence()
    for (i, ev) in enumerate(b.events)
        add_events!(seq, ev, i == 1 ? Event[] : [b.events[i-1]])
    end
    return solve!(seq; method = method)
end
