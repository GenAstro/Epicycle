# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Automatic differentiation of the functions a user writes, used wherever no
# analytic partial was declared. Two seeding modes:
#
#   Argument-seeded — differentiate with respect to a direct argument (y, u, p).
#     The argument's struct must be parameterized on its element type, so a Dual
#     can flow through it. Generators: _make_dynamics_jac_y_ad and its siblings.
#
#   Model-field-seeded — differentiate with respect to one field nested inside a
#     model, identified by a lens recorded when the partial was declared.
#     Primitive: _ad_field_deriv.
#

using Accessors

# ─────────────────────────────────────────────────────────────────────────────
# EvalContext accessors
#
# The type is defined in transcription.jl, because a transcription in any package
# has to build one. What is here is how the framework reads it: `p.model` and
# `get_control(p, t)`.
#
# U and P are parameterized (not fixed to Float64) so ForwardDiff Dual numbers
# can flow through u and params without type-barrier copying. The user never
# constructs EvalContext directly — transcription code builds them internally.
# ─────────────────────────────────────────────────────────────────────────────

get_control(p::EvalContext, t) = p.u
get_param(p::EvalContext)      = p.params

# BoundaryContext and _state_named are in element_interface.jl, beside the
# boundary function machinery that calls them.

# ─────────────────────────────────────────────────────────────────────────────
# AD Jacobian generators
#
# Each generator accepts:
#   dynamics!  — user ODE function f!(dy, y, p::EvalContext, t)
#   ns         — number of states
#
# Returns a jac!(dF, y, ctx, t) closure ready to pass to any transcription's
# defect Jacobian assembly.
#
# Layering contract
# ─────────────────
# Physics layer (these generators):
#   ∂f/∂y  →  ns×ns Matrix        _make_dynamics_jac_y_ad
#   ∂f/∂u  →  ns×nu Matrix        _make_dynamics_jac_u_ad
#   ∂f/∂p  →  ns×np Matrix        _make_dynamics_jac_p_ad
#   ∂f/∂t  →  ns-Vector           _make_dynamics_jac_t_ad  ← per-node time partial
#
# Shooting arc (propagator-level):
#   ∂y_f/∂t_f →  ns-Vector        _make_shooting_jac_tf_ad
#   ∂y_f/∂t_0 →  ns-Vector        _make_shooting_jac_t0_ad  ← non-autonomous correct
#
# Transcription layer (NOT here):
#   ∂defect/∂t0, ∂defect/∂tf assembled from ∂f/∂t via mesh chain rule.
#   Each transcription (LGL, H-S, ZOH) owns this chain rule exactly once.
#   Shooting transcription uses the propagator-level generators directly.
#   Users never write it.
#
# Same pattern for path constraints:
#   ∂g/∂y  →  n_pc×ns Matrix      _make_path_jac_y_ad
#   ∂g/∂u  →  n_pc×nu Matrix      _make_path_jac_u_ad
#   ∂g/∂p  →  n_pc×np Matrix      _make_path_jac_p_ad
#   ∂g/∂t  →  n_pc-Vector         _make_path_jac_t_ad
# ─────────────────────────────────────────────────────────────────────────────

"""
    _make_dynamics_jac_y_ad(dynamics!, ns) → jac!(dF, y, ctx, t)

Returns a Jacobian function ∂f/∂y computed via ForwardDiff.
`dF` is overwritten in-place (ns × ns).
"""
function _make_dynamics_jac_y_ad(dynamics!, ns::Int)
    function jac!(dF, y, ctx, t)
        f = y_ad -> begin
            dy_ad = Vector{eltype(y_ad)}(undef, ns)
            dynamics!(dy_ad, y_ad, ctx, t)
            dy_ad
        end
        dF .= ForwardDiff.jacobian(f, y)
    end
    return jac!
end

"""
    _make_dynamics_jac_u_ad(dynamics!, ns) → jac!(dF, y, ctx, t)

Returns a Jacobian function ∂f/∂u computed via ForwardDiff.
`dF` is overwritten in-place (ns × nu).
If the phase has no controls the function is a no-op (dF stays zero).
"""
function _make_dynamics_jac_u_ad(dynamics!, ns::Int)
    function jac!(dF, y, ctx, t)
        u0 = get_control(ctx, t)
        isempty(u0) && return
        f = u_ad -> begin
            dy_ad = Vector{eltype(u_ad)}(undef, ns)
            ctx_ad = EvalContext(ctx.model, u_ad, ctx.params)
            dynamics!(dy_ad, y, ctx_ad, t)
            dy_ad
        end
        dF .= ForwardDiff.jacobian(f, Vector{Float64}(u0))
    end
    return jac!
end

"""
    _make_dynamics_jac_p_ad(dynamics!, ns) → jac!(dF, y, ctx, t)

Returns a Jacobian function ∂f/∂p (free parameters) computed via ForwardDiff.
`dF` is overwritten in-place (ns × np).
If the phase has no free parameters the function is a no-op (dF stays zero).
"""
function _make_dynamics_jac_p_ad(dynamics!, ns::Int)
    function jac!(dF, y, ctx, t)
        p0 = get_param(ctx)
        isempty(p0) && return
        f = p_ad -> begin
            dy_ad = Vector{eltype(p_ad)}(undef, ns)
            ctx_ad = EvalContext(ctx.model, ctx.u, p_ad)
            dynamics!(dy_ad, y, ctx_ad, t)
            dy_ad
        end
        dF .= ForwardDiff.jacobian(f, Vector{Float64}(p0))
    end
    return jac!
end

"""
    _make_dynamics_jac_t_ad(dynamics!, ns) → jac!(dft, y, ctx, t)

Returns a function computing ∂f/∂t (the per-node time partial) via ForwardDiff.
`dft` is overwritten in-place (ns-Vector, not a Matrix).

Physics layer only.  The transcription assembles ∂defect/∂t0 and ∂defect/∂tf
from this partial via the mesh chain rule  ∂t_k/∂t0 = (1-τ_k),  ∂t_k/∂tf = τ_k.
Users register this function; transcriptions call it — the chain rule is never
the user's responsibility.
"""
function _make_dynamics_jac_t_ad(dynamics!, ns::Int)
    function jac!(dft, y, ctx, t)
        # Seed t as a 1-element vector so jacobian handles the scalar→ns-vector
        # case unambiguously (derivative() expects scalar output).
        J = ForwardDiff.jacobian([t]) do tv
            dy = Vector{eltype(tv)}(undef, ns)
            dynamics!(dy, y, ctx, tv[1])
            dy
        end
        dft .= J[:, 1]
    end
    return jac!
end

# ─────────────────────────────────────────────────────────────────────────────
# Path constraint AD Jacobian generators
#
# Parallel structure to the dynamics generators above.
# g!: (dg::Vector, y::Vector, ctx::EvalContext, t) → nothing
# ─────────────────────────────────────────────────────────────────────────────

"""
    _make_path_jac_y_ad(g!, n_pc) → jac!(dg, y, ctx, t)

AD fallback for ∂g/∂y.  `dg` overwritten in-place (n_pc × ns).
"""
function _make_path_jac_y_ad(g!, n_pc::Int)
    function jac!(dg, y, ctx, t)
        dg .= ForwardDiff.jacobian(y) do y_ad
            g_ad = Vector{eltype(y_ad)}(undef, n_pc)
            g!(g_ad, y_ad, ctx, t)
            g_ad
        end
    end
    return jac!
end

"""
    _make_path_jac_u_ad(g!, n_pc) → jac!(dg, y, ctx, t)

AD fallback for ∂g/∂u.  `dg` overwritten in-place (n_pc × nu).
No-op when nu = 0.
"""
function _make_path_jac_u_ad(g!, n_pc::Int)
    function jac!(dg, y, ctx, t)
        u0 = get_control(ctx, t)
        isempty(u0) && return
        dg .= ForwardDiff.jacobian(Vector{Float64}(u0)) do u_ad
            ctx_ad = EvalContext(ctx.model, u_ad, ctx.params)
            g_ad = Vector{eltype(u_ad)}(undef, n_pc)
            g!(g_ad, y, ctx_ad, t)
            g_ad
        end
    end
    return jac!
end

"""
    _make_path_jac_p_ad(g!, n_pc) → jac!(dg, y, ctx, t)

AD fallback for ∂g/∂p.  `dg` overwritten in-place (n_pc × np).
No-op when np = 0.
"""
function _make_path_jac_p_ad(g!, n_pc::Int)
    function jac!(dg, y, ctx, t)
        p0 = get_param(ctx)
        isempty(p0) && return
        dg .= ForwardDiff.jacobian(Vector{Float64}(p0)) do p_ad
            ctx_ad = EvalContext(ctx.model, ctx.u, p_ad)
            g_ad = Vector{eltype(p_ad)}(undef, n_pc)
            g!(g_ad, y, ctx_ad, t)
            g_ad
        end
    end
    return jac!
end

"""
    _make_path_jac_t_ad(g!, n_pc) → jac!(dgt, y, ctx, t)

AD fallback for ∂g/∂t.  `dgt` overwritten in-place (n_pc-Vector).
Same layering contract as _make_dynamics_jac_t_ad: physics layer only.
"""
function _make_path_jac_t_ad(g!, n_pc::Int)
    function jac!(dgt, y, ctx, t)
        J = ForwardDiff.jacobian([t]) do tv
            g = Vector{eltype(tv)}(undef, n_pc)
            g!(g, y, ctx, tv[1])
            g
        end
        dgt .= J[:, 1]
    end
    return jac!
end

# ─────────────────────────────────────────────────────────────────────────────
# Mode 2 — Model-Field-Seeded AD  (Accessors.jl)
#
# Differentiates f w.r.t. one scalar field inside a carrier struct, identified
# by an Accessors.jl lens.  The lens is stored once at registration time; this
# function is called at every Jacobian evaluation.
#
# Only the targeted field becomes Dual{…}.  All sibling fields in the carrier
# stay Float64 — "surgical" substitution via Accessors.set.
#
# Requirement: every struct on the path from carrier root to the target field
# must be fully parameterized {T<:Real} with a positional inner constructor,
# so ConstructionBase can rebuild it with a mixed-type field tuple.
#
# Usage:
#   lens = @optic _.gravity.mu          # stored once at registration
#   B_mu = _ad_field_deriv(model, lens) do model_d
#       f!(dy, y, u, p, t, model_d)
#       dy
#   end
# ─────────────────────────────────────────────────────────────────────────────

"""
    _ad_field_deriv(f, carrier, lens) → value matching return type of f

Mode 2 AD primitive.  Computes ∂f/∂(carrier.field) via ForwardDiff.derivative
with a surgical Dual substitution via Accessors.set.

  carrier — root struct containing the target field (e.g. model, sc)
  lens    — Accessors.jl lens identifying the field path (@optic _.sub.field)
  f       — closure receiving carrier_dual, returning a scalar or Vector
"""
function _ad_field_deriv(f, carrier, lens)
    v0 = lens(carrier)                          # lenses are callable: lens(obj) → value
    ForwardDiff.derivative(v0) do v_dual
        carrier_dual = Accessors.set(carrier, lens, v_dual)   # one field → Dual; rest Float64
        f(carrier_dual)
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Mode 2 — framework-level generators
#
# These wrap _ad_field_deriv into the same jac!(dF, y, ctx, t) interface used
# by all Mode 1 generators, so transcriptions call them identically.
#
# The lens is captured at generator-call time (registration).  At Jacobian
# evaluation time the current ctx.model is used as the carrier.
#
# Each generator returns a column for ONE model field (scalar lens).
# For a length-nf vector field use a separate lens per component or extend
# to _make_dynamics_jac_modelfield_vec_ad (future).
#
# Result shapes:
#   dynamics: ns-Vector   (∂f/∂model_field, one column)
#   path:     n_pc-Vector (∂g/∂model_field, one column)
#
# These are "column generators" — one per solve-for model field.  The
# transcription assembles them into a matrix column during defect Jacobian
# assembly (exactly as it does for scalar parameter columns).
# ─────────────────────────────────────────────────────────────────────────────

"""
    _make_dynamics_jac_modelfield_ad(dynamics!, ns, lens) → jac!(dF_col, y, ctx, t)

Mode 2 generator for ∂f/∂(model.field) — one scalar model field identified
by `lens` (an Accessors.jl optic, e.g. `@optic _.sub.field`).

Returns a closure `jac!(dF_col, y, ctx, t)` where `dF_col` is an ns-Vector
(one column of the full ∂f/∂model Jacobian).  The carrier is `ctx.model`.

Requirement: every struct on the lens path must be fully `{T<:Real}`-parameterized.
"""
function _make_dynamics_jac_modelfield_ad(dynamics!, ns::Int, lens)
    function jac!(dF_col, y, ctx, t)
        col = _ad_field_deriv(ctx.model, lens) do model_d
            dy = Vector{eltype(typeof(model_d))}(undef, ns)
            ctx_d = EvalContext(model_d, ctx.u, ctx.params)
            dynamics!(dy, y, ctx_d, t)
            dy
        end
        dF_col .= col
    end
    return jac!
end

"""
    _make_path_jac_modelfield_ad(g!, n_pc, lens) → jac!(dg_col, y, ctx, t)

Mode 2 generator for ∂g/∂(model.field) — one scalar model field.

Returns `jac!(dg_col, y, ctx, t)` where `dg_col` is an n_pc-Vector.
"""
function _make_path_jac_modelfield_ad(g!, n_pc::Int, lens)
    function jac!(dg_col, y, ctx, t)
        col = _ad_field_deriv(ctx.model, lens) do model_d
            g = Vector{eltype(typeof(model_d))}(undef, n_pc)
            ctx_d = EvalContext(model_d, ctx.u, ctx.params)
            g!(g, y, ctx_d, t)
            g
        end
        dg_col .= col
    end
    return jac!
end
