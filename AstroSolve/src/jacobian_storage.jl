# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# One storage layer for every registered Jacobian, keyed by target and tag.
#

# ─────────────────────────────────────────────────────────────────────────────
# Registered Jacobians — one storage layer
#
# A Jacobian is registered against a (problem function, solver variable) pair.
# That pairing is the mechanism the framework is committed to, and it is what the
# sparsity pattern is assembled from.
#
# Six kinds of target carry such a table, and each stored, keyed and looked up
# identically through its own accessor triple. This is that, once. Adding a
# seventh means one `_jac_table` method, not three more accessors.
#
# This layer always returns what was stored. What a caller does with it differs,
# and most of the difference is real. A boundary or objective Jacobian is a
# zero-arg closure evaluated at the current point, so its getter calls it. A path
# or dynamics Jacobian is a mutating `fn!(dJ, y, ctx, t)` the caller invokes once
# per node, so its getter returns it. Those are two shapes of derivative, not two
# styles.
#
# 🔴 One difference is not real. `add_jacobian!(sc::SequenceConstraint, var, fn)`
# takes `fn(x)` and closes over `var.value` for the user; the other three take a
# closure the user wrote themselves. That is a user-facing signature difference,
# so unifying it is an interface decision rather than a cleanup.
# ─────────────────────────────────────────────────────────────────────────────

_jac_table(bf::BoundaryFunction)    = bf.jac_fns
_jac_table(pc::PathConstraint)      = pc.jac_fns
_jac_table(o::MayerObjective)       = o.jac_fns
_jac_table(o::BolzaObjective)       = o.jac_fns
_jac_table(p::CollocationPhase)     = p.dynamics_jac
_jac_table(sc::SequenceConstraint)  = sc.jac_fns

# What a Jacobian may be registered against: a solver variable, or the time tag,
# which has no variable of its own.
_jac_key(var::DirectSolverVariable) = objectid(var)
_jac_key(::TimeTag)                 = _TIME_JAC_KEY

register_jacobian!(target, key, fn::Function) =
    (_jac_table(target)[_jac_key(key)] = fn; nothing)

has_registered_jacobian(target, key) = haskey(_jac_table(target), _jac_key(key))

registered_jacobian(target, key) = _jac_table(target)[_jac_key(key)]

# Every Jacobian a target has, for a caller that reports on a whole table rather
# than looks one up. Iterating yields the key as stored, so a caller that walks
# this sees whatever `_jac_key` produced.
registered_jacobians(target) = _jac_table(target)

const _sequence_mga_constraints = Dict{UInt64, Vector{SequenceConstraint}}()

"""
    add_sequence_constraint!(fn, seq, phase1, phase2; lower_bounds, upper_bounds, name="")

Register a cross-phase constraint between two shooting phases.  Do-block form:

    sc = add_sequence_constraint!(seq, phase1, phase2;
             lower_bounds=[0.,0.,0.,r_min], upper_bounds=[0.,0.,0.,Inf],
             name="earth_flyby") do ctx1, ctx2
        [ctx1.tf - ctx2.t0,
         ctx1.mf - ctx2.m0,
         dot(ctx1.vinf_arr, ctx1.vinf_arr) - dot(ctx2.vinf_dep, ctx2.vinf_dep),
         flyby_rp(ctx1.vinf_arr, ctx2.vinf_dep, EARTH)]
    end

    add_jacobian!(sc, arr_var_1) do val
        ...          # val = arr_var_1.value at eval time; return n_con × 3 matrix
    end

Equality constraints use identical lower and upper bounds (typically zeros).
Unregistered Jacobians fall back to ForwardDiff AD.
"""
function add_sequence_constraint!(fn::Function, seq::Sequence,
                                   phase1, phase2;
                                   lower_bounds::AbstractVector{<:Real},
                                   upper_bounds::AbstractVector{<:Real},
                                   name::String = "")
    length(lower_bounds) == length(upper_bounds) ||
        throw(ArgumentError(
            "add_sequence_constraint!: lower_bounds and upper_bounds must have the same " *
            "length; got $(length(lower_bounds)) and $(length(upper_bounds))"))
    sc = SequenceConstraint(phase1, phase2, fn,
                             Float64.(lower_bounds), Float64.(upper_bounds),
                             name, Dict{UInt64,Function}())
    key = objectid(seq)
    if !haskey(_sequence_mga_constraints, key)
        _sequence_mga_constraints[key] = SequenceConstraint[]
    end
    push!(_sequence_mga_constraints[key], sc)
    _invalidate_sm(seq)   # force ShootingManager rebuild with new constraint count
    return sc
end

"""
    add_jacobian!(sc::SequenceConstraint, var, fn)
    add_jacobian!(fn, sc, var)    # do-block sugar

Register an analytic Jacobian for the (SequenceConstraint, variable) pair.
`fn(val)` receives a copy of `var.value` and returns a matrix of shape
(n_con × nlp_length(var)).
"""
# Wraps at registration rather than at lookup. See the 🔴 on the storage layer.
add_jacobian!(sc::SequenceConstraint, var::DirectSolverVariable, fn::Function) =
    register_jacobian!(sc, var, () -> fn(copy(var.value)))

add_jacobian!(fn::Function, sc::SequenceConstraint, var::DirectSolverVariable) =
    add_jacobian!(sc, var, fn)

sequence_mga_constraints(seq::Sequence) =
    get(_sequence_mga_constraints, objectid(seq), SequenceConstraint[])



# add_continuity! lives in oc_manager.jl. The definition here was the same
# body, and the one there carries it plus the mixed-sequence case, which its
# own comment noted as a replication of this one.
