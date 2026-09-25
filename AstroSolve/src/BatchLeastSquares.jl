# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# BatchLeastSquares.jl
#
# Generalized batch least-squares estimator implementing the
# Tapley/Schutz/Born §4.6 algorithm.  Generalization of
# SpringMassEstimator: arbitrary state dimension, real ODE-based
# variational propagation (no LTI shortcut), Cholesky linear solve.
#
# Step-1 scope:
#   · State solve-fors only — but each SV may carry a scalar OR a vector
#     of state components.  The total scalar count across all state-tagged
#     SVs is the state dimension; parameter-tagged SVs are appended via
#     the `parameter_names` keyword (a Tuple of Symbols giving the
#     NamedTuple keys for the closure's `p` argument).
#   · No control solve-fors (u = (;)).
#   · Identity measurement weight by default; user may pass R.
#   · ForwardDiff fallback for state Jacobian if not registered.
#
# Algorithm (Tapley/Schutz/Born Eq. 4.6.4, batch with a priori, weighted):
#
#     Λ_k  = P̄₀⁻¹ + Σ Hᵢᵀ Rᵢ⁻¹ Hᵢ
#     N_k  = P̄₀⁻¹ δX̄₀ + Σ Hᵢᵀ Rᵢ⁻¹ yᵢ
#     δX̂₀  = Λ_k⁻¹ N_k                    (Cholesky solve)
#     X*₀ ← X*₀ + δX̂₀
#     δX̄₀ ← δX̄₀ − δX̂₀
#
# Hᵢ = H̃ᵢ · Φ(tᵢ, t₀); Φ propagated with the dynamics state-Jacobian.

module BatchLeastSquares

using LinearAlgebra
# The variable is AstroSolve's. What an estimator needs of it — how many
# scalars, what it reads now, how to write an iterate back, and its role — is
# the same interface an optimizer needs, so there is one type for both.
using ..AstroSolve: SolverVariable, current_value, assign!, length_of,
                  SolveFor, Fixed
using ..AstroSolve: MeasurementFunction, DynamicsFunction, State, Parameter, has_jacobian, get_jacobian, add_jacobian!
import AstroSolve: length_of
using AstroProp: propagate_with_sensitivities_to_times
using ..Measurements:           SignalPath, TwoWayRange, TwoWayDoppler,
                                AbstractMeasurement,
                                AbstractMeasurementNoise, MeasurementNoise,
                                variance
using ..TrackingDataIO:         ObservationRecord
using AstroModels: GroundStation
using ForwardDiff

import EpicycleBase
using EpicycleBase: AbstractState
using AstroModels:  Spacecraft, get_state
import AstroProp
using AstroProp:    OrbitPropagator, ForceModel, accel_eval!, state_jac!,
                    to_posvel, set_posvel!, JacobianConfig, JacobianResult, eval_jacobian!
using AstroEpochs:  Time

export solve_batch_ls!, BatchLSResult, ODProblem, CartesianStateVar,
       build_od_closures

# ── Resolving a Jacobian ─────────────────────────────────────────────────────
#
# The propagator takes a Jacobian rather than looking one up, so resolving it is
# the caller's job. An analytic partial if one was declared, automatic
# differentiation if not. These are the only two places in estimation that ask
# a dynamics bundle what it knows.

"Analytic state Jacobian if the bundle carries one, otherwise an AD closure."
function _state_jacobian_callable(dyn_fn::DynamicsFunction)
    has_jacobian(dyn_fn, State()) && return get_jacobian(dyn_fn, State())
    return function (y, u, p, t, model)
        ForwardDiff.jacobian(yv -> collect(dyn_fn.f(yv, u, p, t, model)), y)
    end
end

"""
Analytic parameter Jacobian if the bundle carries one, otherwise an AD closure
that perturbs the values of `p_nt` while keeping its field names.
"""
function _param_jacobian_callable(dyn_fn::DynamicsFunction,
                                  p_names::Tuple{Vararg{Symbol}})
    if has_jacobian(dyn_fn, Parameter())
        jac = get_jacobian(dyn_fn, Parameter())
        return function (y, u, p_nt, t, model)
            length(p_names) == 0 && return zeros(length(y), 0)
            return jac(y, u, p_nt, t, model)
        end
    end
    return function (y, u, p_nt, t, model)
        length(p_names) == 0 && return zeros(length(y), 0)
        p0 = collect(values(p_nt))
        return ForwardDiff.jacobian(p0) do pv
            collect(dyn_fn.f(y, u, NamedTuple{p_names}(Tuple(pv)), t, model))
        end
    end
end


"""
    BatchLSResult

The converged state, covariance, residuals, and status from a batch least squares solve.

# Fields
- `X_hat::Vector{Float64}`: Converged state and parameter vector.
- `P_hat::Matrix{Float64}`: Formal covariance of the solve-for components, `inv(Λ_final)`.
- `sigma::Vector{Float64}`: Formal standard deviations, `sqrt.(diag(P_hat))`.
- `corr::Matrix{Float64}`: Correlation matrix derived from `P_hat`.
- `residuals::Vector{Vector{Float64}}`: Postfit residual for each observation, at the converged
  reference trajectory.
- `iters::Int`: Number of iterations executed.
- `converged::Bool`: Whether the relative correction fell below `tol`.
- `solve_for_idx::Vector{Int}`: Indices of estimated components in `X_hat`.
"""
struct BatchLSResult
    X_hat         :: Vector{Float64}        # full state+param vector, length NS+NP
    P_hat         :: Matrix{Float64}        # n_solve × n_solve
    sigma         :: Vector{Float64}        # length n_solve
    corr          :: Matrix{Float64}        # n_solve × n_solve
    residuals     :: Vector{Vector{Float64}}
    iters         :: Int
    converged     :: Bool
    solve_for_idx :: Vector{Int}            # indices into X_hat that were estimated
end

"""
    solve_batch_ls!(svs, obs_times, obs_data, dyn_fn, meas_fn;
                    model, parameter_names = (),
                    R = nothing,
                    n_iters = 10, tol = 1e-10, t0 = 0.0,
                    verbose = false)
        -> BatchLSResult

Estimate state and optional model parameters with batch least squares.

# Arguments
- `svs::AbstractVector{<:SolverVariable}`: Variables in declaration
  order.  Each must be `SolveFor` or `Fixed`.  The flat global vector
  `X` is the concatenation of `current_value(sv)` (length
  `length_of(sv)` each).  The first `NS = sum(length_of) - NP` slots
  are state; the trailing `NP = length(parameter_names)` slots are
  parameters.
- `obs_times`, `obs_data`: Sorted observation times (`>= t0`) and per-epoch
  measurement vectors.
- `dyn_fn::DynamicsFunction`: Dynamics function with signature
  `f(y::AbstractVector, u, p, t, model)`.
- `meas_fn::MeasurementFunction`: Measurement function with the same signature.

# Keyword arguments
- `model`: Model passed to the dynamics and measurement functions.
- `parameter_names::Tuple{Vararg{Symbol}}`: NamedTuple keys for the
  `p` argument; their values come from the parameter slots of `X`.
- `R`: Shared measurement covariance. Defaults to the identity matrix.
- `R_per_obs`: Optional covariance matrix for each observation.
- `n_iters::Int`: Maximum iteration count.
- `tol::Real`: Relative correction tolerance.
- `t0::Real`: Initial integration time.
- `verbose::Bool`: Whether to log iteration progress.

# Notes
Every variable must have role `SolveFor` or `Fixed`. The observation times and
data must have equal lengths. `R` and `R_per_obs` are alternatives and cannot
be supplied together.

# Returns
A `BatchLSResult` containing the estimate, covariance, residuals, and convergence status.

# Example
<!-- doc-fragment -->
```julia
result = solve_batch_ls!(variables, obs_times, obs_data, dynamics, measurement;
                         model = model, R = covariance)
```
"""
function solve_batch_ls!(svs::AbstractVector{<:SolverVariable},
                          obs_times::AbstractVector{<:Real},
                          obs_data::AbstractVector{<:AbstractVector{<:Real}},
                          dyn_fn::DynamicsFunction,
                          meas_fn::MeasurementFunction;
                          model,
                          parameter_names::Tuple{Vararg{Symbol}} = (),
                          R = nothing,
                          R_per_obs::Union{Nothing, AbstractVector} = nothing,
                          n_iters::Int = 10,
                          tol::Real = 1e-10,
                          t0::Real = 0.0,
                          verbose::Bool = false)

    # ── Validation ──────────────────────────────────────────────────────────
    NP = length(parameter_names)

    for s in svs
        (s.role isa SolveFor || s.role isa Fixed) ||
            throw(ArgumentError(
                "solve_batch_ls!: every variable must have role SolveFor or Fixed; " *
                "$(s.name) has $(typeof(s.role))"))
    end
    sv_lens = Int[length_of(s) for s in svs]
    NT      = sum(sv_lens)
    NS      = NT - NP
    NS >= 0 || throw(ArgumentError(
        "solve_batch_ls!: parameter_names names $(NP) parameters, which exceeds the " *
        "$(NT) scalars the variables carry in total"))
    length(obs_times) == length(obs_data) || throw(ArgumentError(
        "solve_batch_ls!: obs_times and obs_data must be the same length; got " *
        "$(length(obs_times)) and $(length(obs_data))"))

    # SV slot ranges in the full (state+param) global vector, length NT.
    sv_ranges = Vector{UnitRange{Int}}(undef, length(svs))
    let off = 0
        for i in eachindex(svs)
            sv_ranges[i] = (off + 1):(off + sv_lens[i])
            off += sv_lens[i]
        end
    end

    # Partition into solve-for vs fixed.
    is_solve = Bool[s.role isa SolveFor for s in svs]
    # Global (NT-length) indices that participate in the estimator.
    solve_for_idx = Int[]
    for i in eachindex(svs)
        if is_solve[i]
            append!(solve_for_idx, collect(sv_ranges[i]))
        end
    end
    n_solve = length(solve_for_idx)
    n_solve > 0 || throw(ArgumentError(
        "solve_batch_ls!: at least one variable must have role SolveFor; all " *
        "$(length(svs)) supplied are Fixed, so there is nothing to estimate"))

    # For each SolveFor SV, its range inside the reduced solve-for vector
    # (length n_solve), used for assembling P_bar_inv and assign!.
    sv_solve_ranges = Vector{UnitRange{Int}}(undef, length(svs))
    let off = 0
        for i in eachindex(svs)
            if is_solve[i]
                sv_solve_ranges[i] = (off + 1):(off + sv_lens[i])
                off += sv_lens[i]
            else
                sv_solve_ranges[i] = 1:0   # empty
            end
        end
    end

    # The arc is a coast from t0 through the last observation: the dynamics take no input.
    u_coast = (;)

    m = length(obs_data[1])
    R_mat = R === nothing ? Matrix{Float64}(I, m, m) : Matrix{Float64}(R)
    R_inv = inv(R_mat)

    # Per-observation R (overrides `R` when given).  Length must match
    # `obs_data`; each entry is m × m and gets inverted once up-front.
    R_inv_per_obs = if R_per_obs === nothing
        nothing
    else
        length(R_per_obs) == length(obs_data) ||
            throw(ArgumentError(
                "solve_batch_ls!: R_per_obs must have one entry per observation, so " *
                "length $(length(obs_data)); got $(length(R_per_obs))"))
        [inv(Matrix{Float64}(Ri)) for Ri in R_per_obs]
    end

    # State Jacobian for the measurement (∂G/∂y) — registered or AD fallback.
    H̃y_of = if has_jacobian(meas_fn, State())
        get_jacobian(meas_fn, State())
    else
        function (y, u, p, t, model)
            ForwardDiff.jacobian(yv -> collect(meas_fn.f(yv, u, p, t, model)), y)
        end
    end

    # Parameter Jacobian for the measurement (∂G/∂p) — registered or AD fallback.
    H̃p_of = if NP > 0 && has_jacobian(meas_fn, Parameter())
        let jac = get_jacobian(meas_fn, Parameter())
            (y, u, p_nt, t, model) -> jac(y, u, p_nt, t, model)
        end
    else
        function (y, u, p_nt, t, model)
            if NP == 0
                return zeros(m, 0)
            end
            p0 = collect(values(p_nt))
            return ForwardDiff.jacobian(p0) do pv
                pn = NamedTuple{parameter_names}(Tuple(pv))
                collect(meas_fn.f(y, u, pn, t, model))
            end
        end
    end

    # ── A priori (variances → P̄₀⁻¹) — only over solve-for slots ─────────────
    P_bar_inv = zeros(n_solve, n_solve)
    for i in eachindex(svs)
        is_solve[i] || continue
        a   = svs[i].covariance
        rng = sv_solve_ranges[i]
        if a === nothing
            continue
        elseif a isa Real
            for k in rng
                P_bar_inv[k, k] = 1.0 / Float64(a)
            end
        elseif a isa AbstractVector
            length(a) == length(rng) || throw(ArgumentError(
                "solve_batch_ls!: covariance for $(svs[i].name) must have one entry per " *
                "component, so length $(length(rng)); got $(length(a))"))
            for (kk, k) in enumerate(rng)
                P_bar_inv[k, k] = 1.0 / Float64(a[kk])
            end
        elseif a isa AbstractMatrix
            size(a) == (length(rng), length(rng)) || throw(ArgumentError(
                "solve_batch_ls!: covariance for $(svs[i].name) must be " *
                "($(length(rng)), $(length(rng))) to match its components; got $(size(a))"))
            P_bar_inv[rng, rng] = inv(Matrix{Float64}(a))
        else
            throw(ArgumentError(
                "solve_batch_ls!: covariance for $(svs[i].name) must be a Real variance, a " *
                "vector of variances, or a covariance matrix; got $(typeof(a))"))
        end
    end

    # ── Reference state from ALL variables (Fixed + SolveFor) ──────────────
    # X_star_0 has length NT and feeds the propagator's full y0/p initial
    # values.  Only the entries indexed by `solve_for_idx` are updated by the
    # iteration; Fixed slots are held constant.
    X_star_0 = zeros(Float64, NT)
    for i in eachindex(svs)
        v = current_value(svs[i])
        if v isa Real
            X_star_0[sv_ranges[i][1]] = Float64(v)
        else
            X_star_0[sv_ranges[i]] .= Float64.(collect(v))
        end
    end
    delta_X_bar_0 = zeros(n_solve)

    Λ     = zeros(n_solve, n_solve)
    N_vec = zeros(n_solve)

    state_slots = 1:NS
    param_slots = (NS + 1):NT

    # Map global slot indices into the reduced solve-for vector for column
    # selection from full Φ_y / Φ_p.
    solve_for_state_cols = [k for k in solve_for_idx if k ≤ NS]
    solve_for_param_cols = [k - NS for k in solve_for_idx if k > NS]

    converged = false
    iters_run = 0

    for iter in 1:n_iters
        iters_run = iter
        Λ      .= P_bar_inv
        N_vec  .= P_bar_inv * delta_X_bar_0

        y0    = collect(@view X_star_0[state_slots])
        p_nt  = NamedTuple{parameter_names}(Tuple(@view X_star_0[param_slots]))

        ys, Φys, Φps = propagate_with_sensitivities_to_times(
            dyn_fn.f, _state_jacobian_callable(dyn_fn),
            _param_jacobian_callable(dyn_fn, parameter_names),
            y0, t0, obs_times, u_coast, p_nt, model)

        for i in eachindex(obs_times)
            t       = Float64(obs_times[i])
            y_vec   = ys[i]
            Φy      = Φys[i]
            Φp      = Φps[i]
            u_active = u_coast

            G   = meas_fn(y_vec, u_active, p_nt, t, model)
            H̃y = H̃y_of(y_vec, u_active, p_nt, t, model)
            H̃p = H̃p_of(y_vec, u_active, p_nt, t, model)

            # Full per-iter Jacobian wrt full (state, param) vector, then
            # slice columns that correspond to solve-for slots only.
            H_full = zeros(m, NT)
            H_full[:, state_slots] .= H̃y * Φy
            if NP > 0
                H_full[:, param_slots] .= H̃y * Φp .+ H̃p
            end
            H = H_full[:, solve_for_idx]

            r   = Float64.(obs_data[i]) .- collect(G)
            R_inv_i = R_inv_per_obs === nothing ? R_inv : R_inv_per_obs[i]
            HtR = H' * R_inv_i
            Λ      .+= HtR * H
            N_vec  .+= HtR * r
        end

        # Cholesky solve.  Falls back to `\` if Λ is not PD numerically.
        local δX̂_0
        try
            C = cholesky(Symmetric(Λ))
            δX̂_0 = C \ N_vec
        catch
            δX̂_0 = Λ \ N_vec
        end

        rel = norm(δX̂_0) / max(norm(@view X_star_0[solve_for_idx]), 1.0)
        if verbose
            @info "iter $(iter): ‖x̂‖=$(norm(δX̂_0))  rel=$(rel)"
        end

        # Update only solve-for slots; Fixed slots stay at their initial value.
        @views X_star_0[solve_for_idx] .+= δX̂_0
        delta_X_bar_0 .-= δX̂_0

        if rel < tol
            converged = true
            break
        end
    end

    # ── Push converged estimate back through the variable contract ─────────
    # Only SolveFor SVs receive `assign!`; Fixed SVs are left untouched.
    for i in eachindex(svs)
        is_solve[i] || continue
        rng = sv_ranges[i]
        if length(rng) == 1
            assign!(svs[i], X_star_0[rng[1]])
        else
            assign!(svs[i], X_star_0[rng])
        end
    end

    # ── Postfit residuals at the converged reference ───────────────────────
    y0    = collect(@view X_star_0[state_slots])
    p_nt  = NamedTuple{parameter_names}(Tuple(@view X_star_0[param_slots]))
    ys, _, _ = propagate_with_sensitivities_to_times(
        dyn_fn.f, _state_jacobian_callable(dyn_fn),
        _param_jacobian_callable(dyn_fn, parameter_names),
        y0, t0, obs_times, u_coast, p_nt, model)
    residuals = Vector{Vector{Float64}}(undef, length(obs_times))
    for i in eachindex(obs_times)
        t = Float64(obs_times[i])
        u_active = u_coast
        G        = meas_fn(ys[i], u_active, p_nt, t, model)
        residuals[i] = Float64.(obs_data[i]) .- collect(G)
    end

    P_hat = Matrix(inv(Symmetric(Λ)))
    σ     = sqrt.(diag(P_hat))
    corr  = zeros(n_solve, n_solve)
    for j in 1:n_solve, i in 1:n_solve
        corr[i, j] = P_hat[i, j] / (σ[i] * σ[j])
    end

    return BatchLSResult(copy(X_star_0), P_hat, σ, corr,
                         residuals, iters_run, converged, copy(solve_for_idx))
end

# ─────────────────────────────────────────────────────────────────────────────
# Model-based interface — Spacecraft / OrbitPropagator / Measurements
# ─────────────────────────────────────────────────────────────────────────────
#
# A `CartesianStateVar` is a value-carrying variable struct whose identity
# means "the Cartesian (pos, vel) state of a Spacecraft".  Wrap it in a
# `ModelVariable(var, sc)` and the solver's iterate flows back into the
# spacecraft via `set_field!` → `set_posvel!`.

"""
    CartesianStateVar{T}(value::Vector{T})

Six-element Cartesian state variable bound to a `Spacecraft` via
`ModelVariable`.  Reads/writes round-trip through `to_posvel` /
`set_posvel!`.
"""
mutable struct CartesianStateVar{T<:Real} <: AbstractState
    value::Vector{T}
end

length_of(::CartesianStateVar) = 6

# Bind to a Spacecraft via the EpicycleBase get/set protocol.
EpicycleBase.get_field(sc::Spacecraft, ::CartesianStateVar)         = to_posvel(sc)
EpicycleBase.set_field!(sc::Spacecraft, ::CartesianStateVar, x::AbstractVector) =
    (set_posvel!(sc, x); nothing)

"""
    ODProblem(; spacecraft, propagator, measurements)

The spacecraft, propagation model, measurements, and variables used for orbit determination.

# Fields
- `spacecraft::Spacecraft`: the state lives here, and the solver updates it.
- `propagator::OrbitPropagator`: the dynamics live here, under any force model.
- `measurements::Vector{AbstractMeasurement}`: the measurement specs the estimator dispatches over
  to build predictors. Each spec carries its own `noise::MeasurementNoise`, and the estimator reads
  the measurement covariance as the square of `meas.noise.sigma`.
- `solve_for::Vector{Any}`: what is estimated, as the `Vary` declarations. An `ODProblem` holds its
  own variables, the way a `Sequence` holds them in its phases and events.

# Example
<!-- doc-fragment -->
```julia
range_meas   = TwoWayRange(SignalPath(gs, sc, gs);   noise = MeasurementNoise(1.0e-3))
doppler_meas = TwoWayDoppler(SignalPath(gs, sc, gs); noise = MeasurementNoise(1.0e-5))

problem = ODProblem(
    spacecraft   = sc,
    propagator   = prop,
    measurements = [range_meas, doppler_meas],
)
```
"""
struct ODProblem
    spacecraft   :: Spacecraft
    propagator   :: OrbitPropagator
    measurements :: Vector{AbstractMeasurement}
    solve_for    :: Vector{Any}      # what is estimated; a Sequence holds its
                                     # own variables in phases and events, and
                                     # an ODProblem should hold its own too
end

ODProblem(; spacecraft::Spacecraft,
            propagator::OrbitPropagator,
            measurements::AbstractVector,
            solve_for::AbstractVector = Any[]) =
    ODProblem(spacecraft, propagator,
              collect(AbstractMeasurement, measurements),
              collect(Any, solve_for))

# ── Measurement predictors — model-based, dispatch on spec type ───────────

"""
    _predict(spec, y, t::Time, sc) -> Real

Predicted observable for measurement `spec` given the current
spacecraft state `y` (length-6 Cartesian, ICRF) at epoch `t`.
The `sc` template is used only for type-based identification of the
satellite participant inside a `SignalPath`.
"""
function _predict(spec::TwoWayRange, y::AbstractVector, t::Time, sc::Spacecraft)
    path = spec.path
    n    = length(path)
    n >= 2 || throw(ArgumentError(
        "range prediction needs a path with at least two participants; got $(n)"))
    r_prev = _participant_pos(path.participants[1], y, t)
    total  = zero(eltype(r_prev))
    @inbounds for i in 2:n
        r_next = _participant_pos(path.participants[i], y, t)
        total += norm(r_next - r_prev)
        r_prev = r_next
    end
    return total
end

"""
    _predict(::TwoWayDoppler, y, t::Time, sc) -> Real

Time derivative of the `TwoWayRange` predictor along the same path
— two-way Doppler in km/s.  Each leg `r_a -> r_b` contributes
`(r_b - r_a)·(v_b - v_a) / |r_b - r_a|`; legs are summed.
"""
function _predict(spec::TwoWayDoppler, y::AbstractVector, t::Time, sc::Spacecraft)
    path = spec.path
    n    = length(path)
    n >= 2 || throw(ArgumentError(
        "Doppler prediction needs a path with at least two participants; got $(n)"))
    r_prev = _participant_pos(path.participants[1], y, t)
    v_prev = _participant_vel(path.participants[1], y, t)
    total  = zero(eltype(r_prev))
    @inbounds for i in 2:n
        r_next = _participant_pos(path.participants[i], y, t)
        v_next = _participant_vel(path.participants[i], y, t)
        Δr = r_next - r_prev
        Δv = v_next - v_prev
        total += dot(Δr, Δv) / norm(Δr)
        r_prev = r_next
        v_prev = v_next
    end
    return total
end

# Per-participant position at time `t` — ground stations come from
# their geodetic frame; the spacecraft's position comes from the
# integrator's current `y` (so AD perturbations through the iterate
# reach the predictor correctly).
_participant_pos(p::GroundStation, y::AbstractVector, t::Time) = get_state(p, t)[1]
_participant_pos(::Spacecraft,     y::AbstractVector, ::Time) = @inbounds (@view y[1:3])

_participant_vel(p::GroundStation, y::AbstractVector, t::Time) = get_state(p, t)[2]
_participant_vel(::Spacecraft,     y::AbstractVector, ::Time) = @inbounds (@view y[4:6])

# ── High-level solve_batch_ls! ────────────────────────────────────────────

"""
    build_od_closures(records, model::ODProblem)
        -> (dyn_fn, meas_fn, obs_times, obs_data, R_per_obs)

Build the universal-signature dynamics and measurement closures, the
observation schedule (seconds since `model.spacecraft.time`), the
per-record scalar observations, and a per-observation measurement
covariance vector `R_per_obs[i] = [σ_i^2;;]` (1 × 1).  σ is read off
the matching spec's `meas.noise.sigma` (per F4 — noise lives on the
function).

Multiple measurement kinds are allowed in `records`; for each kind
appearing in the record stream there must be exactly one matching
spec in `model.measurements`.

Used by both `solve_batch_ls!` and `ExtendedKalmanFilter.run_ekf!`.

# Returns
The five pieces an estimator steps over, as a tuple `(dyn, meas, obs_times, obs_data, R_per_obs)`:
the dynamics closure, the measurement closure, the observation times as seconds from the
spacecraft's epoch, the observed values, and the measurement covariance for each observation.

Throws an `ArgumentError` when `records` is empty. The records are sorted in place by receive time,
since TDM segments arrive grouped by station and a filter needs them interleaved.

# Example
<!-- doc-fragment -->
```julia
dyn, meas, obs_times, obs_data, R_per_obs = build_od_closures(records, problem)
```
"""
function build_od_closures(records::AbstractVector{<:ObservationRecord},
                            model::ODProblem)

    isempty(records) && throw(ArgumentError(
        "build_od_closures: at least one observation record is required; got none"))

    sc     = model.spacecraft
    epoch0 = sc.time

    # TDM segments arrive grouped by participant_1 (e.g. all Goldstone, then
    # all Madrid, then all Canberra), so `t_receive` is non-monotone across
    # segment boundaries.  The EKF requires strictly increasing propagator
    # times; sort the records in place by `t_receive` so per-station streams
    # interleave correctly and the downstream ε-disambiguation only collapses
    # genuine ties (e.g. range + Doppler at the same epoch).
    if records isa AbstractVector
        sort!(records; by = r -> Float64(r.t_receive - epoch0))
    end

    # Resolve, per (kind, participant_1) appearing in the stream, the
    # spec instance and σ² — both pulled off the spec itself (F4: noise
    # lives on the function).  `variance(spec.noise)` lets the abstract
    # noise hierarchy supply R without each estimator knowing the
    # concrete subtype.  Keying on `participant_1` is what makes
    # multi-station OD work: a TDM with Goldstone- and Madrid-range
    # records picks the matching `TwoWayRange` spec by station name.
    keys_seen = unique((r.measurement_type, r.participant_1) for r in records)
    spec_for     = Dict{Tuple{Symbol,String}, AbstractMeasurement}()
    variance_for = Dict{Tuple{Symbol,String}, Float64}()
    for key in keys_seen
        kind, part1 = key
        spec_type = _spec_type_for_tdm(kind)
        idx = findfirst(model.measurements) do m
            m isa spec_type || return false
            # Empty participant_1 on a record (legacy / untagged) matches
            # any spec of the right type — preserves single-station
            # back-compat for callers that haven't yet started tagging.
            isempty(part1) && return true
            length(m.path.names) >= 1 && String(m.path.names[1]) == part1
        end
        idx === nothing && throw(ArgumentError(
            "build_od_closures: the data has TDM kind $(repr(kind)), which needs a " *
            "$(spec_type) whose first participant is $(repr(part1)); model.measurements " *
            "has none. Add that measurement to the model, or drop the records."))
        spec = model.measurements[idx]
        spec_for[key]     = spec
        variance_for[key] = variance(spec.noise)
    end

    obs_times = Float64[(r.t_receive - epoch0) * 86400.0 for r in records]
    obs_data  = [Float64[r.observed]                     for r in records]

    # Per-observation 1×1 covariance.
    R_per_obs = [reshape([variance_for[(r.measurement_type, r.participant_1)]],
                         1, 1) for r in records]

    # Co-located observations of different kinds (e.g. range + Doppler at the
    # same epoch) need distinct keys in the per-time spec lookup below, but we
    # still want the propagator to see strictly increasing times.  Disambiguate
    # by adding a sub-microsecond ε per duplicate; ε ≪ any meaningful sample
    # rate, so the dynamics see effectively the same instant.
    let last_t = -Inf, ε = 1.0e-9
        for i in eachindex(obs_times)
            if obs_times[i] <= last_t
                obs_times[i] = last_t + ε
            end
            last_t = obs_times[i]
        end
    end

    # Per-time spec lookup so the universal-signature meas closure can
    # dispatch the right predictor without an obs-index argument.  Times
    # are unique by construction after the disambiguation above.
    spec_by_time = Dict{Float64, AbstractMeasurement}()
    for (i, r) in pairs(records)
        spec_by_time[obs_times[i]] = spec_for[(r.measurement_type, r.participant_1)]
    end

    # Both go through AstroProp, which sums the forces' accelerations and writes the kinematic rows
    # once. Calling each force into one buffer kept only the last force's acceleration, and summing
    # each force's state_jac! counted ∂ṙ/∂v once per force; a force with no analytic Jacobian
    # takes the automatic-differentiation fallback instead of throwing.
    dyn = DynamicsFunction(name = :epicycle_forces) do y, _u, _p, t_rel, m
        t  = epoch0 + t_rel / 86400.0
        return AstroProp._eval_all!(m.propagator.forces, t, y, zeros(eltype(y), 6), m.spacecraft)
    end
    jac_result = JacobianResult(JacobianConfig(partial_y = true); n_state = 6,
                                forces = model.propagator.forces, sc = model.spacecraft,
                                t_example = epoch0)
    add_jacobian!(dyn, State()) do y, _u, _p, t_rel, m
        t = epoch0 + t_rel / 86400.0
        eval_jacobian!(jac_result, m.propagator.forces, y, m.spacecraft, t)
        return copy(jac_result.partial_y)
    end

    meas = MeasurementFunction(name = :epicycle_measurement) do y, _u, _p, t_rel, m
        spec = spec_by_time[Float64(t_rel)]
        t    = epoch0 + t_rel / 86400.0
        return [_predict(spec, y, t, m.spacecraft)]
    end

    return dyn, meas, obs_times, obs_data, R_per_obs
end

"""
    solve_batch_ls!(svs, records; model::ODProblem,
                    n_iters = 10, tol = 1e-10, verbose = false)
        -> BatchLSResult

Model-based entry point.  Pulls dynamics off `model.propagator.forces`
and predictors off `model.measurements` — the use case writes no
physics or geometry.

`records` is a vector of `ObservationRecord` (e.g. straight from
`read_records`).  Time-since-epoch is computed from
`r.t_receive - model.spacecraft.time`.  Each record contributes one
scalar observation; `R = σ²` is taken from the matching spec's
`meas.noise` field.
"""
function solve_batch_ls!(svs::AbstractVector{<:SolverVariable},
                          records::AbstractVector{<:ObservationRecord};
                          model::ODProblem,
                          n_iters::Int = 10,
                          tol::Real    = 1e-10,
                          verbose::Bool = false)

    dyn, meas, obs_times, obs_data, R_per_obs = build_od_closures(records, model)

    return solve_batch_ls!(svs, obs_times, obs_data, dyn, meas;
                           model     = model,
                           R_per_obs = R_per_obs,
                           n_iters   = n_iters,
                           tol       = tol,
                           t0        = 0.0,
                           verbose   = verbose)
end

# TDM measurement-type Symbol → measurement spec type.  Extend as new
# measurement kinds come online.
const _TDM_KIND_TO_SPEC = Dict{Symbol, Type}(
    :RANGE   => TwoWayRange,
    :DOPPLER => TwoWayDoppler,
)

# §9.4: a small closed set enumerates its options in the message rather than sending the
# engineer back to the docstring.
_spec_type_for_tdm(s::Symbol) = get(_TDM_KIND_TO_SPEC, s) do
    throw(ArgumentError(
        "TDM measurement type must be one of " *
        "$(join(sort(string.(collect(keys(_TDM_KIND_TO_SPEC)))), ", ")); got $(repr(s))"))
end

end # module BatchLeastSquares
