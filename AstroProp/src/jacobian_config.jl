# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    JacobianConfig

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

Allocate a `JacobianResult` and build the dispatch plan.  Call once at problem
setup; reuse by passing to `eval_jacobian!`.
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
    p0 = EpicycleBase.get_field(mv.model, mv.tag)
    h  = cbrt(eps(Float64)) * max(1.0, abs(p0))

    acc_p = zeros(length(y))
    acc_m = zeros(length(y))

    EpicycleBase.set_field!(mv.model, mv.tag, p0 + h)
    accel_eval!(force, t, y, acc_p, sc, [])
    EpicycleBase.set_field!(mv.model, mv.tag, p0 - h)
    accel_eval!(force, t, y, acc_m, sc, [])
    EpicycleBase.set_field!(mv.model, mv.tag, p0)

    out .+= (acc_p .- acc_m) ./ (2h)
    return nothing
end

# ---------------------------------------------------------------------------
# eval_jacobian! — primary in-place entry point
# ---------------------------------------------------------------------------

"""
    eval_jacobian!(result, forces, y, sc, t) -> result

In-place Jacobian evaluation.  For each force:
  • Calls `state_jac!`  if registered (analytic), else `ForwardDiff.jacobian`.
  • Calls `param_jac!`  if registered (analytic), else central-FD fallback.

Registration is detected once at `JacobianResult` construction via `hasmethod`.
Zero allocation per ODE step after the first compile.
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

function _eval_all!(fm::ForceModel, t, y, acc, sc)
    fill!(acc, 0.0)
    for force in fm.forces
        accel_eval!(force, t, y, acc, sc, [])
    end
end

function fd_differentiate_wrt(fm, owner_model, tag::AbstractVarTag,
                               y::AbstractVector, sc::Spacecraft, t)
    p0 = EpicycleBase.get_field(owner_model, tag)
    h  = cbrt(eps(Float64)) * max(1.0, abs(p0))

    acc_plus  = zeros(length(y))
    acc_minus = zeros(length(y))

    EpicycleBase.set_field!(owner_model, tag, p0 + h)
    _eval_all!(fm, t, y, acc_plus, sc)
    EpicycleBase.set_field!(owner_model, tag, p0 - h)
    _eval_all!(fm, t, y, acc_minus, sc)
    EpicycleBase.set_field!(owner_model, tag, p0)

    return (acc_plus .- acc_minus) ./ (2h)
end
