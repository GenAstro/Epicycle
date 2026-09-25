# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# Batch least-squares estimator specialized to the spring-mass problem of Tapley,
# Schutz & Born, *Statistical Orbit Determination*, §4.8.2.
#
# Algorithm: Tapley/Schutz/Born Eq. 4.6.4, batch with a priori, identity weighting (R = I).
#
# It takes `SolverVariable`s through the prototype contract and never inspects what they wrap,
# which is what lets the direct-style and model-style use cases reuse this core verbatim.
#
#     Λ_k  = Σ Hᵢᵀ Hᵢ + P̄₀⁻¹
#     N_k  = Σ Hᵢᵀ yᵢ + P̄₀⁻¹ δX̄₀
#     δX̂₀  = Λ_k⁻¹ N_k
#     X*₀ ← X*₀ + δX̂₀         (shift reference)
#     δX̄₀ ← δX̄₀  − δX̂₀         (shift a priori deviation)
#
# Closed-form linear-system propagation; analytic Φ(t,0) and H̃ from the book.
#
# Inputs assume a 2-element solve-for list ordered [x0, v0], each with
# `covariance::Real` interpreted as a *variance* (so P̄₀ = diag(cov_x0,
# cov_v0)).

module SpringMassEstimator

using LinearAlgebra
using ..AstroSolve: SolverVariable
using ..AstroSolve: MeasurementFunction, DynamicsFunction, current_value,
                  assign!, length_of, SolveFor, State, get_jacobian

export solve_spring_mass_batch!

"""
    SpringMassResult

The estimated initial state, covariance, residuals, and iteration count from the spring-mass fit.

# Fields
- `X_hat::Vector{Float64}`: Converged `[x̂₀, v̂₀]`.
- `P_hat::Matrix{Float64}`: Formal covariance, `inv(Λ_final)`.
- `sigma::Vector{Float64}`: Formal standard deviations from `sqrt.(diag(P_hat))`.
- `corr::Float64`: Correlation `P_hat[1,2] / (sigma[1] * sigma[2])`.
- `residuals::Vector{Vector{Float64}}`: Length-two residual at each observation time, evaluated
                 at the converged reference.
- `iters::Int`: Number of iterations performed.
"""
struct SpringMassResult
    X_hat     :: Vector{Float64}
    P_hat     :: Matrix{Float64}
    sigma     :: Vector{Float64}
    corr      :: Float64
    residuals :: Vector{Vector{Float64}}
    iters     :: Int
end

"""
    solve_spring_mass_batch!(svs, obs_times, obs_data, dyn_fn, meas_fn;
                             params = (;), n_iters = 4)
        -> SpringMassResult

Fit the initial state of a linear spring-mass system to range and range-rate observations.

# Arguments
- `svs::Vector{<:SolverVariable}` — exactly `[sv_x0, sv_v0]`.  Each must be
  `SolveFor` and have a scalar `covariance` (interpreted as variance).
- `obs_times::Vector{Float64}`    — observation epochs (s).
- `obs_data::Vector{Vector{Float64}}` — each entry is `[ρ, ρ̇]` at that epoch.
- `dyn_fn::DynamicsFunction`      — dynamics bundle.
  Closure: `(y, u, p, t, model) -> ẋ::AbstractVector` (universal signature).
  Must have an analytic `State()` Jacobian: `(y, u, p, t, model) -> A::Matrix`.
  Spring-mass is LTI so `A` is constant; the estimator evaluates it once
  and uses `Φ = exp(t·A)` for both reference propagation and STM.
- `meas_fn::MeasurementFunction`  — measurement bundle.
  Closure: `(y, u, p, t, model) -> z::AbstractVector`.
  Must have analytic `State()` Jacobian: `(y, u, p, t, model) -> H̃::Matrix`.
  This problem has no controls or solve-for parameters, so the estimator
  passes `u = (;)` and `p = (;)`.
- `model`                         — the live model object (or any struct /
                                    NamedTuple holding the constants the
                                    closures need).  Passed to every
                                    closure as the second positional arg.
- `n_iters::Int = 4`              — iteration count. Fixed rather than convergence-tested,
                                    which matches the worked example in the reference.

# Returns
A `SpringMassResult` carrying the converged state, its formal covariance, the standard
deviations and correlation read off that covariance, the postfit residuals at every
observation epoch, and the iteration count.

# Notes
The system is linear and time-invariant, so the dynamics Jacobian `A` is evaluated once and
`Φ(t) = exp(tA)` serves as both the reference trajectory propagator and the state transition
matrix. That is exact here and is why this harness has analytic truth: noise-free observations
generated from a known initial state are recovered to numerical precision rather than
approximately.

Weighting is identity, so `R = I` and the observations are treated as equally precise. A
weighted fit is the general `solve_batch_ls!`, not this one.

Throws `ArgumentError` unless `svs` is exactly two scalar `SolveFor` variables.

# Side effects
After convergence, the estimator calls `assign!(sv, x̂₀_block)` on each `SolverVariable`, so the
variable structs — and any model objects bound to them through `ModelVariable` — end up holding
the converged estimate.

# Example
<!-- doc-fragment -->
```julia
# svs are [x0, v0] with a priori variances; dyn_fn and meas_fn carry analytic Jacobians.
result = solve_spring_mass_batch!(svs, obs_times, obs_data, dyn_fn, meas_fn; model = m)
result.X_hat        # converged [x̂₀, v̂₀]
result.sigma        # formal standard deviations
```
"""
function solve_spring_mass_batch!(svs::AbstractVector{<:SolverVariable},
                                   obs_times::AbstractVector{<:Real},
                                   obs_data::AbstractVector{<:AbstractVector{<:Real}},
                                   dyn_fn::DynamicsFunction,
                                   meas_fn::MeasurementFunction;
                                   model, n_iters::Int = 4)

    length(svs) == 2 || throw(ArgumentError(
        "solve_spring_mass_batch!: svs must be exactly [x0, v0], so length 2; " *
        "got $(length(svs))"))
    all(s -> s.role isa SolveFor, svs) || throw(ArgumentError(
        "solve_spring_mass_batch!: every entry of svs must have role SolveFor; got " *
        "$(join([string(typeof(s.role)) for s in svs], ", "))"))
    all(s -> length_of(s) == 1, svs) || throw(ArgumentError(
        "solve_spring_mass_batch!: each solve-for must be scalar, so length 1; got " *
        "$(join([string(length_of(s)) for s in svs], ", "))"))

    A_of = get_jacobian(dyn_fn,  State())
    H̃_of = get_jacobian(meas_fn, State())

    # Universal signature: every closure is called as f(y, u, p, t, model).
    # This problem has no controls and no parameter solve-fors, so the
    # framework passes empty NamedTuples for u and p.
    u_empty = (;)
    p_empty = (;)

    # ── A priori (variances → P̄₀⁻¹) ────────────────────────────────────
    P_bar_inv = zeros(2, 2)
    P_bar_inv[1, 1] = 1.0 / Float64(svs[1].covariance)
    P_bar_inv[2, 2] = 1.0 / Float64(svs[2].covariance)

    # ── Initial reference state X*₀ from the variables ──────────────────
    X_star_0      = Float64[ current_value(svs[1]), current_value(svs[2]) ]
    delta_X_bar_0 = zeros(2)

    # Spring-mass is LTI: A is independent of state and time.
    # Evaluate once at the initial reference, t = 0.
    A = A_of((x = X_star_0[1], v = X_star_0[2]), u_empty, p_empty, 0.0, model)

    Λ = zeros(2, 2)
    N = zeros(2)

    for _ in 1:n_iters
        Λ .= P_bar_inv
        N .= P_bar_inv * delta_X_bar_0

        for i in eachindex(obs_times)
            t = Float64(obs_times[i])

            # State transition matrix and reference trajectory at t,
            # both from the user-supplied dynamics Jacobian.
            Φ      = exp(t * A)
            x_vec  = Φ * X_star_0
            state  = (x = x_vec[1], v = x_vec[2])

            # Predicted measurement and its state Jacobian (user-supplied).
            G  = meas_fn(state, u_empty, p_empty, t, model)
            H̃  = H̃_of(state, u_empty, p_empty, t, model)

            # H = H̃ · Φ — maps initial-state perturbations to obs perturbations.
            H = H̃ * Φ

            yobs = obs_data[i]
            r    = Float64.(yobs) .- G

            # Λ += Hᵀ H, N += Hᵀ r
            Λ .+= H' * H
            N .+= H' * r
        end

        δX̂_0 = Λ \ N

        X_star_0      .+= δX̂_0
        delta_X_bar_0 .-= δX̂_0
    end

    # ── Push converged estimate back through the variable contract ─────
    assign!(svs[1], X_star_0[1])
    assign!(svs[2], X_star_0[2])

    # ── Postfit residuals at the converged reference ───────────────────
    residuals = Vector{Vector{Float64}}(undef, length(obs_times))
    for i in eachindex(obs_times)
        t     = Float64(obs_times[i])
        Φ     = exp(t * A)
        x_vec = Φ * X_star_0
        state = (x = x_vec[1], v = x_vec[2])
        G     = meas_fn(state, u_empty, p_empty, t, model)
        residuals[i] = Float64.(obs_data[i]) .- G
    end

    P_hat = inv(Λ)
    σ     = sqrt.(diag(P_hat))
    corr  = P_hat[1, 2] / (σ[1] * σ[2])

    return SpringMassResult(copy(X_star_0), P_hat, σ, corr, residuals, n_iters)
end

end # module SpringMassEstimator
