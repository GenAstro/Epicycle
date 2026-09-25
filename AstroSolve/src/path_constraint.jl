# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# PathConstraint, BoundaryConstraint, and the Constraint form taking a BoundaryFunction.
#

# ─────────────────────────────────────────────────────────────────────────────
# PathConstraint
#
# A constraint evaluated at every mesh node: g(y_k, u_k, t_k) ∈ [lb, ub]
# fn!: (g::Vector, y::Vector, p, t::Float64) → nothing    (in-place, g pre-zeroed)
# jac_fns: Dict{UInt64, Function} keyed by objectid(var)
#   each fn!: (dg::Matrix, y::Vector, p, t::Float64) → nothing  (in-place, dg pre-zeroed)
#   size of dg: (n_pc, n_states) for state_var, (n_pc, n_controls) for control_var
# ─────────────────────────────────────────────────────────────────────────────

struct PathConstraint
    fn          ::Function                # fn!(g, y, p, t)
    n_pc        ::Int                     # number of constraint outputs
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    jac_fns     ::Dict{UInt64, Function}  # keyed by objectid(var)
    name        ::String
end

# do-block sugar: PathConstraint(n, lb, ub; name="") do g, y, p, t ... end
PathConstraint(fn::Function, n_pc::Int;
               lower_bounds::Vector{Float64},
               upper_bounds::Vector{Float64},
               name::String = "") =
    PathConstraint(fn, n_pc, lower_bounds, upper_bounds, Dict{UInt64,Function}(), name)

struct PathConstraintBlock
    pc::PathConstraint
end

set_path_jacobian!(pc::PathConstraint, var::DirectSolverVariable, fn::Function) =
    register_jacobian!(pc, var, fn)
set_path_jacobian!(fn::Function, pc::PathConstraint, var::DirectSolverVariable) =
    set_path_jacobian!(pc, var, fn)
has_path_jacobian(pc::PathConstraint, var::DirectSolverVariable) =
    has_registered_jacobian(pc, var)
get_path_jacobian(pc::PathConstraint, var::DirectSolverVariable) =
    registered_jacobian(pc, var)

"""
    set_path_jacobian!(fn, pc, TimeTag())
    set_path_jacobian!(pc, TimeTag(), fn)

Register an analytic ∂g/∂t for the path constraint.  Physics layer only.

Signature:  fn(dgt::Vector, y, ctx, t) → nothing  (fills n_pc-Vector in-place)

The transcription assembles ∂g_k/∂t0 and ∂g_k/∂tf from this via the mesh
chain rule — identical contract to add_jacobian!(fn, phase, TimeTag()).
"""
function set_path_jacobian!(pc::PathConstraint, ::TimeTag, fn::Function)
    register_jacobian!(pc, TimeTag(), fn)
    return nothing
end
set_path_jacobian!(fn::Function, pc::PathConstraint, ::TimeTag) =
    set_path_jacobian!(pc, TimeTag(), fn)
has_path_jac_t(pc::PathConstraint) = has_registered_jacobian(pc, TimeTag())
get_path_jac_t(pc::PathConstraint) = registered_jacobian(pc, TimeTag())

# ─────────────────────────────────────────────────────────────────────────────
# BoundaryConstraint: a constraint on a phase's endpoints, from a BoundaryFunction
# ─────────────────────────────────────────────────────────────────────────────

struct BoundaryConstraint
    calc         ::BoundaryFunction
    lower_bounds ::Vector{Float64}
    upper_bounds ::Vector{Float64}
    scale        ::Vector{Float64}   # g_nlp = (g_phys - shift) / scale  (default: 1)
    shift        ::Vector{Float64}   # g_nlp = (g_phys - shift) / scale  (default: 0)
end

# Backwards-compatible 3-arg constructor (no scaling)
BoundaryConstraint(calc, lb, ub) =
    BoundaryConstraint(calc, lb, ub, ones(length(lb)), zeros(length(lb)))

# The keyword form builds a BoundaryConstraint from a BoundaryFunction. It is
# distinct from `Constraint(quantity, subject, deps...)`, which constrains a
# quantity on a subject and takes positional arguments.
#
# It used to fall through to a `Constraint(; calc = <a Calc>)` method for
# anything that was not a BoundaryFunction. That method was retired, so the
# branch could only raise; it says so now instead of raising somewhere else.
function Constraint(; calc, lower_bounds=nothing, upper_bounds=nothing,
                      scale=nothing, shift=nothing)
    calc isa BoundaryFunction || throw(ArgumentError(
        "Constraint(; calc = ...) takes a BoundaryFunction. To constrain a " *
        "quantity on a subject, write Constraint(quantity, subject; equals = ...)."))
    n  = length(something(lower_bounds, upper_bounds))
    lb = lower_bounds === nothing ? fill(-Inf, n) : Float64.(lower_bounds)
    ub = upper_bounds === nothing ? fill( Inf, n) : Float64.(upper_bounds)
    sc = scale === nothing ? ones(n)  : Float64.(scale)
    sh = shift === nothing ? zeros(n) : Float64.(shift)
    return BoundaryConstraint(calc, lb, ub, sc, sh)
end
