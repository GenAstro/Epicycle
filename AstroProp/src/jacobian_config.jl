# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    JacobianConfig(; partial_y = false, partial_p = ModelVariable[])

Declares which Jacobian blocks to compute.  Does not allocate any result storage.

# Fields
- `partial_y::Bool`                     — request df/dy (n×n state Jacobian)
- `partial_p::Vector{ModelVariable}`    — list of parameter variables to differentiate
"""
struct JacobianConfig
    partial_y::Bool
    partial_p::Vector{ModelVariable}

    function JacobianConfig(; partial_y::Bool = false, partial_p = ModelVariable[])
        return new(partial_y, collect(ModelVariable, partial_p))
    end
end

# ---------------------------------------------------------------------------
# DispatchPlan — built once at JacobianResult construction via hasmethod.
# Encodes which force types have registered analytic state_jac! / param_jac!
# methods.  Checked at setup time so the inner ODE loop is branch-free.
# ---------------------------------------------------------------------------

struct DispatchPlan
    analytic_state::Set{DataType}                 # force types with state_jac! defined
    analytic_param::Set{Tuple{DataType,DataType}} # (force_type, tag_type) pairs
end

function DispatchPlan(forces::ForceModel, cfg::JacobianConfig,
                      sc::Spacecraft, t_example)
    analytic_state = Set{DataType}()
    analytic_param = Set{Tuple{DataType,DataType}}()

    TT = typeof(t_example)
    TS = typeof(sc)

    for force in forces.forces
        FT = typeof(force)

        if cfg.partial_y
            if hasmethod(state_jac!,
                    Tuple{Matrix{Float64}, FT, TT, AbstractVector, TS})
                push!(analytic_state, FT)
            end
        end

        for mv in cfg.partial_p
            GT = typeof(mv.tag)
            if hasmethod(param_jac!,
                    Tuple{Vector{Float64}, FT, GT, TT, AbstractVector, TS})
                push!(analytic_param, (FT, GT))
            end
        end
    end

    return DispatchPlan(analytic_state, analytic_param)
end

"""
    JacobianResult

Pre-allocated container for one Jacobian evaluation.  Allocate once via
`JacobianResult(cfg; n_state, forces, sc, t_example)` at problem setup, then
reuse across all ODE steps via `eval_jacobian!`.

The `DispatchPlan` is built at construction time via `hasmethod` — zero cost
per ODE step.

# Fields
- `config::JacobianConfig`
- `partial_y::Matrix{Float64}`                       — df/dy (n×n)
- `partial_p::Dict{ModelVariable, Vector{Float64}}`  — df/dp_i (n×1) per variable
- `plan::DispatchPlan`                               — analytic vs. AD routing
"""
struct JacobianResult
    config::JacobianConfig
    partial_y::Matrix{Float64}
    partial_p::Dict{ModelVariable, Vector{Float64}}
    plan::DispatchPlan
end

"""
    JacobianResult(cfg; n_state, forces, sc, t_example) -> JacobianResult

Allocate a `JacobianResult` and work out which partial each force contributes. Call it once at
problem setup and reuse the result by passing it to `eval_jacobian!`, which then evaluates into
the buffers rather than allocating new ones.
"""
function JacobianResult(cfg::JacobianConfig;
                         n_state::Int,
                         forces::ForceModel,
                         sc::Spacecraft,
                         t_example)
    partial_y = cfg.partial_y ? zeros(n_state, n_state) :
                                Matrix{Float64}(undef, 0, 0)
    partial_p = Dict{ModelVariable, Vector{Float64}}(
        mv => zeros(n_state) for mv in cfg.partial_p)
    plan = DispatchPlan(forces, cfg, sc, t_example)
    return JacobianResult(cfg, partial_y, partial_p, plan)
end

# ---------------------------------------------------------------------------
# AD fallbacks — called when no analytic method is registered
# ---------------------------------------------------------------------------

# State fallback: ForwardDiff.jacobian on accel_eval! for this force only.
function _state_jac_ad!(out, force, t, y::AbstractVector, sc::Spacecraft)
    f_ = y_ -> begin
        acc = zeros(eltype(y_), length(y_))
        accel_eval!(force, t, y_, acc, sc, [])
        acc
    end
    out .+= ForwardDiff.jacobian(f_, y)
    return nothing
end

# Param fallback: central FD on accel_eval! for this force/tag pair only.
# Used when no param_jac! method is registered and no kernel is available to
# inject a Dual.  Correct but approximate — register param_jac! for production.
function _param_jac_fd!(out, force, mv::ModelVariable, t, y::AbstractVector,
                         sc::Spacecraft)
    out .+= _central_difference(mv.model, mv.tag) do
        acc = zeros(length(y))
        accel_eval!(force, t, y, acc, sc, [])
        acc
    end
    return nothing
end

# Central difference of `f()` with respect to the scalar field `tag` on `model`. The one
# implementation behind the parameter-Jacobian fallback and `fd_differentiate_wrt`.
#
# h = cbrt(eps) · max(1, |p₀|) balances truncation against round-off for a second-order central
# difference on a well-scaled parameter. The field is restored in `finally`, so a model that
# throws mid-evaluation does not leave the caller's model perturbed by h.
function _central_difference(f, model, tag::AbstractVarTag)
    p0 = get_field(model, tag)
    h  = cbrt(eps(Float64)) * max(1.0, abs(p0))
    try
        set_field!(model, tag, p0 + h);  f_plus  = f()
        set_field!(model, tag, p0 - h);  f_minus = f()
        return (f_plus .- f_minus) ./ (2h)
    finally
        set_field!(model, tag, p0)
    end
end

# ---------------------------------------------------------------------------
# eval_jacobian! — primary in-place entry point
# ---------------------------------------------------------------------------

"""
    eval_jacobian!(result, forces, y, sc, t) -> result

In-place Jacobian evaluation.  For each force:
  • Calls `state_jac!`  if registered (analytic), else `ForwardDiff.jacobian`.
  • Calls `param_jac!`  if registered (analytic), else central-FD fallback.

Registration is detected once at `JacobianResult` construction via `hasmethod`, so no method
lookup happens per ODE step. The result arrays are reused; the automatic-differentiation and
finite-difference fallbacks still allocate on each call.
"""
function eval_jacobian!(result::JacobianResult, fm::ForceModel,
                         y::AbstractVector, sc::Spacecraft, t)
    cfg  = result.config
    plan = result.plan

    if cfg.partial_y
        fill!(result.partial_y, 0.0)
    end
    for mv in cfg.partial_p
        fill!(result.partial_p[mv], 0.0)
    end

    for force in fm.forces
        FT = typeof(force)

        # A = ∂f/∂y
        if cfg.partial_y
            if FT in plan.analytic_state
                state_jac!(result.partial_y, force, t, y, sc)
            else
                _state_jac_ad!(result.partial_y, force, t, y, sc)
            end
        end

        # B_i = ∂f/∂p_i
        for mv in cfg.partial_p
            GT  = typeof(mv.tag)
            out = result.partial_p[mv]
            if (FT, GT) in plan.analytic_param
                param_jac!(out, force, mv.tag, t, y, sc)
            else
                _param_jac_fd!(out, force, mv, t, y, sc)
            end
        end
    end

    # Each force's contribution, analytic or AD, carries the kinematic rows ṙ = v as well as its
    # acceleration, so summing N forces would give ∂ṙ/∂v = N·I and ∂ṙ/∂p = N·0. The kinematics
    # appear once in the dynamics, so they are written once here.
    if cfg.partial_y
        result.partial_y[1:3, :] .= 0.0
        for i in 1:3
            result.partial_y[i, i + 3] = 1.0
        end
    end
    for mv in cfg.partial_p
        result.partial_p[mv][1:3] .= 0.0
    end

    return result
end

"""
    eval_jacobian(cfg, forces, y, sc, t; n_state) -> JacobianResult

Allocating convenience wrapper.  For REPL use and tests.
"""
function eval_jacobian(cfg::JacobianConfig, fm::ForceModel,
                        y::AbstractVector, sc::Spacecraft, t;
                        n_state::Int = length(y))
    result = JacobianResult(cfg; n_state = n_state, forces = fm,
                             sc = sc, t_example = t)
    eval_jacobian!(result, fm, y, sc, t)
    return result
end

# ---------------------------------------------------------------------------
# fd_differentiate_wrt — AstroProp-level FD oracle for testing
# ---------------------------------------------------------------------------

function _eval_all!(fm::OrbitODE, t, y, acc, sc)
    fill!(acc, 0.0)
    accel_eval!(fm, t, y, acc, sc, [])
end

# Accelerations summed and kinematics written once, as `propagate!` assembles them. Each force's
# `accel_eval!` overwrites the rows it writes, so calling them into one buffer keeps only the last.
function _eval_all!(fm::ForceModel, t, y, acc, sc)
    fill!(acc, 0.0)
    acc[1:3] .= y[4:6]
    one = zeros(eltype(acc), length(acc))
    for force in fm.forces
        fill!(one, 0.0)
        accel_eval!(force, t, y, one, sc, [])
        acc[4:6] .+= one[4:6]
    end
    return acc
end

"""
    fd_differentiate_wrt(rhs_fn, model, tag::AbstractVarTag) -> Vector{Float64}
    fd_differentiate_wrt(fm, model, tag::AbstractVarTag, y, sc::Spacecraft, t) -> Vector{Float64}

Central finite-difference derivative with respect to the scalar field on `model` identified by
`tag`, used as an oracle for checking an analytic `param_jac!` method.

The first form differentiates a zero-argument closure `rhs_fn`. The second differentiates the
total acceleration of force model `fm` (a `ForceModel` or `OrbitODE`) at state `y` [km, km/s] and
epoch `t` for spacecraft `sc`.

# Arguments
- `rhs_fn`: Zero-argument closure returning the quantity to differentiate, as a vector.
- `fm`: Force model whose acceleration is differentiated.
- `model`: The object owning the field, such as a `CelestialBody` for `Mu()`.
- `tag`: Tag identifying the field.
- `y`: Cartesian state [km, km/s].
- `sc`: Spacecraft the forces act on.
- `t`: Epoch of the evaluation.

# Notes
The field is perturbed by h = cbrt(eps) · max(1, |p₀|) either side and restored afterwards, also
when the evaluation throws. `model` must have `get_field` and `set_field!` methods for `tag`.

# Returns
The derivative as a vector, the length of `rhs_fn()` or of `y`.

# Examples
```julia
using AstroProp, AstroUniverse
body = CelestialBody("Earth", 398600.4415, 6378.137, 1/298.257223563, 399)
fd_differentiate_wrt(() -> [2.0 * body.mu], body, Mu())    # ≈ [2.0]
```
"""
fd_differentiate_wrt(rhs_fn, model, tag::AbstractVarTag) = _central_difference(rhs_fn, model, tag)

function fd_differentiate_wrt(fm, model, tag::AbstractVarTag,
                               y::AbstractVector, sc::Spacecraft, t)
    return _central_difference(model, tag) do
        acc = zeros(length(y))
        _eval_all!(fm, t, y, acc, sc)
        acc
    end
end
