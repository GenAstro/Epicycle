# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# Continuous-discrete Extended Kalman Filter using UDUᵀ covariance factorization.
# Tapley, Schutz & Born, *Statistical Orbit Determination*, §4.7 and §5.7.
# State variables may be `SolveFor` or `Fixed`; parameter solve-fors are not supported.
# Process noise is attached to each solve-for variable, and measurement covariance may be dense.

module ExtendedKalmanFilter

using LinearAlgebra
# The variable is AstroSolve's. What an estimator needs of it — how many
# scalars, what it reads now, how to write an iterate back, and its role — is
# the same interface an optimizer needs, so there is one type for both.
using ..AstroSolve: SolverVariable, current_value, assign!, length_of,
                  SolveFor, Fixed
using ..AstroSolve: MeasurementFunction, DynamicsFunction, State, has_jacobian, get_jacobian
using AstroProp: propagate_with_stm
using ..BatchLeastSquares: _state_jacobian_callable
using ..UDU: udu_from_P, udu_to_P,
             thornton_time_update!, bierman_measurement_update!
using ..ProcessNoiseModels: ProcessNoiseModel, NoNoise, discretize
using ..BatchLeastSquares:    ODProblem, build_od_closures
using ..TrackingDataIO:       ObservationRecord
using ForwardDiff

export run_ekf!, EKFRecord, EKFResult,
       EKFState, init_ekf, time_update!, measurement_update!,
       current_state, current_covariance, current_time, commit_to_svs!,
       RTSResult, run_rts, run_iterated_rts!, IteratedRTSResult

"""
    EKFRecord

Per-observation record from `run_ekf!`.

# Fields
- `t::Float64`: Filter time.
- `y_pre::Vector{Float64}`: Propagated full state before the update.
- `y_post::Vector{Float64}`: Full state after the update.
- `P_pre::Matrix{Float64}`: Solve-for covariance before the update.
- `P_post::Matrix{Float64}`: Solve-for covariance after the update.
- `Phi::Matrix{Float64}`: State transition matrix from the
                previous record's epoch to `t` (identity for the
                first record and for any zero-duration step).
- `prefit::Vector{Float64}`: Prefit residual `z - z_hat_pre`.
- `postfit::Vector{Float64}`: Postfit residual `z - z_hat_post`.
- `K::Matrix{Float64}`: Kalman gain.
"""
struct EKFRecord
    t       :: Float64
    y_pre   :: Vector{Float64}
    y_post  :: Vector{Float64}
    P_pre   :: Matrix{Float64}
    P_post  :: Matrix{Float64}
    Phi     :: Matrix{Float64}
    prefit  :: Vector{Float64}
    postfit :: Vector{Float64}
    K       :: Matrix{Float64}
end

"""
    EKFResult

The final estimate, covariance, and update history from `run_ekf!`.

# Fields
- `X_hat::Vector{Float64}`: Final full state vector.
- `P_hat::Matrix{Float64}`: Final solve-for covariance.
- `sigma::Vector{Float64}`: Formal standard deviations, `sqrt.(diag(P_hat))`.
- `records::Vector{EKFRecord}`: Update history, one record per observation.
- `solve_for_idx::Vector{Int}`: Indices of estimated components in `X_hat`.
"""
struct EKFResult
    X_hat         :: Vector{Float64}
    P_hat         :: Matrix{Float64}
    sigma         :: Vector{Float64}
    records       :: Vector{EKFRecord}
    solve_for_idx :: Vector{Int}
end

# ─────────────────────────────────────────────────────────────────────────────
# Stateful EKF API
# ─────────────────────────────────────────────────────────────────────────────

"""
Mutable filter state used by the stateful API (`init_ekf`,
`time_update!`, `measurement_update!`, …).  Carries the SV layout, the
UDU factors, the running full state vector, and the current epoch.
Constructed by `init_ekf`; the user normally only sees it as an opaque
handle.
"""
mutable struct EKFState
    dyn_fn          :: DynamicsFunction
    meas_fn         :: MeasurementFunction
    model           :: Any
    u_const         :: Any
    R               :: Matrix{Float64}
    L_R             :: Matrix{Float64}
    svs             :: Vector{SolverVariable}
    sv_ranges       :: Vector{UnitRange{Int}}
    sv_solve_ranges :: Vector{UnitRange{Int}}
    is_solve        :: Vector{Bool}
    solve_for_idx   :: Vector{Int}
    n_solve         :: Int
    NS              :: Int
    pn_models       :: Vector{ProcessNoiseModel}
    pn_ranges       :: Vector{UnitRange{Int}}
    H̃y_of           :: Any
    y               :: Vector{Float64}
    U               :: Matrix{Float64}
    D               :: Vector{Float64}
    t               :: Float64
    Phi_last        :: Matrix{Float64}
end

"""
    init_ekf(svs, dyn_fn, meas_fn; model, R, u = (;), t0 = 0.0) -> EKFState

Build an `EKFState` from the SV declaration, the dynamics/measurement
function bundles, the model, and the (default) measurement covariance
`R`.  Validates the SV roles and covariance, builds the initial state
vector and UDU factorisation, and caches the Cholesky factor of `R`.

After this call the user drives the filter explicitly:

    time_update!(ekf, t)
    info = measurement_update!(ekf, z)

# Arguments
- `svs::AbstractVector{<:SolverVariable}`: State variables in declaration order.
- `dyn_fn::DynamicsFunction`: Dynamics function used for propagation.
- `meas_fn::MeasurementFunction`: Measurement prediction function.
- `model`: Model passed to the dynamics and measurement functions.
- `R`: Default measurement covariance matrix.
- `u`: Constant control passed to the dynamics. Defaults to an empty named tuple.
- `t0::Real`: Initial filter time.

# Notes
Every variable must have role `SolveFor` or `Fixed` and carry a valid covariance.
The measurement covariance must be positive definite.

# Returns
The filter state, which [`time_update!`](@ref) and [`measurement_update!`](@ref) then advance. It
carries the estimate, its covariance in factored form, and the current epoch.

# Example
<!-- doc-fragment -->
```julia
dyn, meas, obs_times, obs_data, R_per_obs = build_od_closures(records, problem)

ekf = init_ekf([y0], dyn, meas; model = problem, R = R_per_obs[1], t0 = 0.0)
```
"""
function init_ekf(svs::AbstractVector{<:SolverVariable},
                  dyn_fn::DynamicsFunction,
                  meas_fn::MeasurementFunction;
                  model,
                  R,
                  u = (;),
                  t0::Real = 0.0)

    for s in svs
        (s.role isa SolveFor || s.role isa Fixed) ||
            throw(ArgumentError(
                "init_ekf: every variable must have role SolveFor or Fixed; " *
                "$(s.name) has $(typeof(s.role))"))
    end
    sv_lens = Int[length_of(s) for s in svs]
    NS      = sum(sv_lens)

    sv_ranges = Vector{UnitRange{Int}}(undef, length(svs))
    let off = 0
        for i in eachindex(svs)
            sv_ranges[i] = (off + 1):(off + sv_lens[i])
            off += sv_lens[i]
        end
    end

    is_solve = Bool[s.role isa SolveFor for s in svs]
    solve_for_idx = Int[]
    for i in eachindex(svs)
        if is_solve[i]
            append!(solve_for_idx, collect(sv_ranges[i]))
        end
    end
    n_solve = length(solve_for_idx)
    n_solve > 0 || throw(ArgumentError(
        "init_ekf: at least one variable must have role SolveFor; all " *
        "$(length(svs)) supplied are Fixed, so the filter has nothing to estimate"))

    sv_solve_ranges = Vector{UnitRange{Int}}(undef, length(svs))
    let off = 0
        for i in eachindex(svs)
            if is_solve[i]
                sv_solve_ranges[i] = (off + 1):(off + sv_lens[i])
                off += sv_lens[i]
            else
                sv_solve_ranges[i] = 1:0
            end
        end
    end

    # Initial covariance from a-priori.
    P0 = zeros(n_solve, n_solve)
    for i in eachindex(svs)
        is_solve[i] || continue
        a   = svs[i].covariance
        rng = sv_solve_ranges[i]
        if a === nothing
            throw(ArgumentError(
                "init_ekf: $(svs[i].name) has no covariance, and a filter needs one to " *
                "start from. Give it as `covariance` on the Vary that declares it"))
        elseif a isa Real
            for k in rng
                P0[k, k] = Float64(a)
            end
        elseif a isa AbstractVector
            length(a) == length(rng) || throw(ArgumentError(
                "init_ekf: covariance for $(svs[i].name) must have one entry per " *
                "component, so length $(length(rng)); got $(length(a))"))
            for (kk, k) in enumerate(rng)
                P0[k, k] = Float64(a[kk])
            end
        elseif a isa AbstractMatrix
            size(a) == (length(rng), length(rng)) || throw(ArgumentError(
                "init_ekf: covariance for $(svs[i].name) must be " *
                "($(length(rng)), $(length(rng))) to match its components; got $(size(a))"))
            P0[rng, rng] = Matrix{Float64}(a)
        else
            throw(ArgumentError(
                "init_ekf: covariance for $(svs[i].name) must be a Real variance, a vector " *
                "of variances, or a covariance matrix; got $(typeof(a))"))
        end
    end

    # Process-noise models on SolveFor SVs.
    pn_models = ProcessNoiseModel[]
    pn_ranges = UnitRange{Int}[]
    for i in eachindex(svs)
        is_solve[i] || continue
        pn = svs[i].process_noise
        if pn === nothing
            push!(pn_models, NoNoise())
        elseif pn isa ProcessNoiseModel
            push!(pn_models, pn)
        else
            throw(ArgumentError(
                "init_ekf: process_noise on $(svs[i].name) must be a ProcessNoiseModel, " *
                "such as NoNoise() or DiagonalSNC(psd); got $(typeof(pn))"))
        end
        push!(pn_ranges, sv_solve_ranges[i])
    end

    # Initial state vector from ALL SVs.
    y = zeros(Float64, NS)
    for i in eachindex(svs)
        v = current_value(svs[i])
        if v isa Real
            y[sv_ranges[i][1]] = Float64(v)
        else
            y[sv_ranges[i]] .= Float64.(collect(v))
        end
    end

    # Measurement Jacobian ∂G/∂y.
    H̃y_of = if has_jacobian(meas_fn, State())
        get_jacobian(meas_fn, State())
    else
        function (yv, uv, pv, tv, mv)
            ForwardDiff.jacobian(z -> collect(meas_fn.f(z, uv, pv, tv, mv)), yv)
        end
    end

    R_mat  = Matrix{Float64}(R)
    L_R    = Matrix(cholesky(Symmetric(R_mat)).L)
    U, D   = udu_from_P(P0)

    return EKFState(dyn_fn, meas_fn, model, u, R_mat, L_R,
                    Vector{SolverVariable}(svs), sv_ranges, sv_solve_ranges,
                    is_solve, solve_for_idx, n_solve, NS,
                    pn_models, pn_ranges, H̃y_of,
                    y, U, D, Float64(t0),
                    Matrix{Float64}(LinearAlgebra.I, n_solve, n_solve))
end

"""
    time_update!(ekf, t_new) -> ekf

Propagate the full state from `ekf.t` to `t_new` via
`propagate_with_stm`, advance the UDU factors with a Thornton update
(picking up process noise from each SolveFor SV's `ProcessNoiseModel`),
and stamp `ekf.t = t_new`.  No-op if `t_new ≤ ekf.t`.
# Returns
The filter state, propagated to `t_new` and mutated in place, so the return value is the same object
that was passed in.

# Example
<!-- doc-fragment -->
```julia
time_update!(ekf, obs_times[k])
```
"""
function time_update!(ekf::EKFState, t_new::Real)
    t_new = Float64(t_new)
    if t_new > ekf.t
        Δt = t_new - ekf.t
        y_new, Φ_full = propagate_with_stm(ekf.dyn_fn.f,
                                     _state_jacobian_callable(ekf.dyn_fn),
                                     ekf.y, ekf.t, t_new,
                                            ekf.u_const, (;), ekf.model)
        ekf.y .= y_new
        Φ_sub = Φ_full[ekf.solve_for_idx, ekf.solve_for_idx]
        ekf.Phi_last = Matrix{Float64}(Φ_sub)

        G_blocks = Matrix{Float64}[]
        q_blocks = Vector{Float64}[]
        for (mdl, rng) in zip(ekf.pn_models, ekf.pn_ranges)
            Gl, ql = discretize(mdl, Δt, length(rng))
            if !isempty(ql)
                Gfull = zeros(ekf.n_solve, length(ql))
                Gfull[rng, :] = Gl
                push!(G_blocks, Gfull)
                push!(q_blocks, ql)
            end
        end
        G_step = isempty(G_blocks) ? zeros(ekf.n_solve, 0) : hcat(G_blocks...)
        Q_step = isempty(q_blocks) ? Float64[]             : vcat(q_blocks...)
        thornton_time_update!(ekf.U, ekf.D, Φ_sub, G_step, Q_step)
        ekf.t = t_new
    else
        ekf.Phi_last = Matrix{Float64}(LinearAlgebra.I, ekf.n_solve, ekf.n_solve)
    end
    return ekf
end

"""
    measurement_update!(ekf, z; R = nothing)
        -> (; prefit, postfit, K, y_pre, y_post, P_pre, P_post)

Apply a measurement update at the current epoch using observation `z`.
If `R` is `nothing` (default) the cached `R` from `init_ekf` is used.
Returns a NamedTuple with prefit/postfit residuals, the equivalent
Kalman gain in the original (non-decorrelated) measurement space, and
pre/post snapshots of the full state and reduced covariance.

# Example
<!-- doc-fragment -->
```julia
time_update!(ekf, obs_times[k])
info = measurement_update!(ekf, obs_data[k]; R = R_per_obs[k])

info.postfit[1]
```
"""
function measurement_update!(ekf::EKFState,
                              z::AbstractVector;
                              R = nothing)
    L_R = if R === nothing
        ekf.L_R
    else
        Matrix(cholesky(Symmetric(Matrix{Float64}(R))).L)
    end

    P_pre = udu_to_P(ekf.U, ekf.D)
    y_pre = copy(ekf.y)

    ẑ        = collect(ekf.meas_fn(ekf.y, ekf.u_const, (;), ekf.t, ekf.model))
    H̃y_full  = ekf.H̃y_of(ekf.y, ekf.u_const, (;), ekf.t, ekf.model)
    H        = H̃y_full[:, ekf.solve_for_idx]

    z_f    = Float64.(z)
    prefit = z_f .- ẑ

    z̃ = L_R \ z_f
    ẑ̃ = L_R \ ẑ
    H̃ = L_R \ H
    ν̃ = z̃ .- ẑ̃

    x_solve     = view(ekf.y, ekf.solve_for_idx)
    x_pre_solve = copy(x_solve)
    m           = length(z_f)
    gain_cols   = Matrix{Float64}(undef, ekf.n_solve, m)
    for j in 1:m
        h_row  = vec(H̃[j, :])
        δx     = x_solve .- x_pre_solve
        ν_seq  = ν̃[j] - dot(h_row, δx)
        gain_j = bierman_measurement_update!(ekf.U, ekf.D, x_solve,
                                              h_row, 1.0, ν_seq)
        gain_cols[:, j] = gain_j
    end
    K_eff = gain_cols / L_R

    y_post  = copy(ekf.y)
    ẑ_post  = collect(ekf.meas_fn(ekf.y, ekf.u_const, (;), ekf.t, ekf.model))
    postfit = z_f .- ẑ_post
    P_post  = udu_to_P(ekf.U, ekf.D)

    return (; prefit, postfit, K = K_eff, y_pre, y_post, P_pre, P_post)
end

"""
    current_state(ekf) -> Vector{Float64}

The filter's current estimate, as of the last update applied to it.

# Returns
A copy of the estimated state, in the units of the dynamics the filter was built with. It is a copy,
so writing to it does not reach the filter.

# Example
<!-- doc-fragment -->
```julia
y = current_state(ekf)
```
"""
current_state(ekf::EKFState) = copy(ekf.y)

"""
    current_covariance(ekf) -> Matrix{Float64}

The covariance of the filter's current estimate, formed from the factors it carries internally.

# Returns
The covariance matrix over the estimated components, in the squared units of the corresponding
states, so the square root of its first diagonal entry is a one-sigma position uncertainty.

# Example
<!-- doc-fragment -->
```julia
sigma_position = sqrt(diag(current_covariance(ekf))[1])
```
"""
current_covariance(ekf::EKFState) = udu_to_P(ekf.U, ekf.D)

"Return the current filter epoch."
current_time(ekf::EKFState) = ekf.t

"""
    commit_to_svs!(ekf) -> ekf

Push the current filter estimate back into the SolveFor SVs via
`assign!`.  Fixed SVs are not touched.
"""
function commit_to_svs!(ekf::EKFState)
    for i in eachindex(ekf.svs)
        ekf.is_solve[i] || continue
        rng = ekf.sv_ranges[i]
        if length(rng) == 1
            assign!(ekf.svs[i], ekf.y[rng[1]])
        else
            assign!(ekf.svs[i], ekf.y[rng])
        end
    end
    return ekf
end

# ─────────────────────────────────────────────────────────────────────────────
# Driver
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_ekf!(svs, obs_times, obs_data, dyn_fn, meas_fn;
             model, R, u = (;), t0 = 0.0, verbose = false)
        -> EKFResult

Continuous-discrete EKF with UDUᵀ covariance.  Same SV interface as
`solve_batch_ls!`:

  · `svs`         — SolverVariables; `SolveFor` are estimated, `Fixed`
                    are propagated but not updated.  The flat state
                    vector `y` is the concatenation of `current_value`
                    over `svs`, in declaration order.
  · `obs_times`, `obs_data` — measurement schedule (sorted, ≥ t0).
  · `dyn_fn`, `meas_fn`     — universal-signature `(y, u, p, t, model)`,
                              with `y::AbstractVector`.
  · `R`           — measurement covariance (m × m, SPD).
  · `u`           — control held constant over the run (any type).
  · `t0`          — initial epoch.
"""
function run_ekf!(svs::AbstractVector{<:SolverVariable},
                   obs_times::AbstractVector{<:Real},
                   obs_data::AbstractVector{<:AbstractVector{<:Real}},
                   dyn_fn::DynamicsFunction,
                   meas_fn::MeasurementFunction;
                   model,
                   R,
                   R_per_obs::Union{Nothing, AbstractVector} = nothing,
                   u = (;),
                   t0::Real = 0.0,
                   verbose::Bool = false)

    length(obs_times) == length(obs_data) || throw(ArgumentError(
        "run_ekf!: obs_times and obs_data must be the same length; got " *
        "$(length(obs_times)) and $(length(obs_data))"))
    if R_per_obs !== nothing
        length(R_per_obs) == length(obs_data) ||
            throw(ArgumentError(
                "run_ekf!: R_per_obs must have one entry per observation, so length " *
                "$(length(obs_data)); got $(length(R_per_obs))"))
    end

    ekf = init_ekf(svs, dyn_fn, meas_fn; model = model, R = R, u = u, t0 = t0)

    records = Vector{EKFRecord}(undef, length(obs_times))
    for k in eachindex(obs_times)
        time_update!(ekf, obs_times[k])
        Φ_k = copy(ekf.Phi_last)
        info = if R_per_obs === nothing
            measurement_update!(ekf, obs_data[k])
        else
            measurement_update!(ekf, obs_data[k]; R = R_per_obs[k])
        end
        records[k] = EKFRecord(ekf.t, info.y_pre, info.y_post,
                               info.P_pre, info.P_post, Φ_k,
                               info.prefit, info.postfit, info.K)
        if verbose
            @info "k=$(k) t=$(ekf.t)  prefit_rms=$(sqrt(mean_sq(info.prefit)))  postfit_rms=$(sqrt(mean_sq(info.postfit)))"
        end
    end

    commit_to_svs!(ekf)
    P_hat = current_covariance(ekf)
    σ     = sqrt.(diag(P_hat))
    return EKFResult(current_state(ekf), P_hat, σ, records,
                     copy(ekf.solve_for_idx))
end

mean_sq(v) = isempty(v) ? 0.0 : sum(abs2, v) / length(v)

# ─────────────────────────────────────────────────────────────────────────────
# Model-based entry point — Spacecraft / OrbitPropagator / Measurements
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_ekf!(svs, records; model::ODProblem, verbose = false) -> EKFResult

Model-based EKF entry point, mirroring the high-level
`solve_batch_ls!(svs, records; model)`.  Pulls dynamics off
`model.propagator` and predictor + σ off each spec's `meas.noise`
field in `model.measurements` (F4) — the use case writes no physics
or geometry.

Drives the standard EKF loop over `records`: for each
`ObservationRecord`, time-update to `r.t_receive`, then measurement-
update with the scalar `r.observed`.  After the last record the
estimate is committed back to the SolveFor SVs via `assign!`.
"""
function run_ekf!(svs::AbstractVector{<:SolverVariable},
                   records::AbstractVector{<:ObservationRecord};
                   model::ODProblem,
                   verbose::Bool = false)

    dyn, meas, obs_times, obs_data, R_per_obs = build_od_closures(records, model)

    return run_ekf!(svs, obs_times, obs_data, dyn, meas;
                    model     = model,
                    R         = R_per_obs[1],   # init_ekf needs a placeholder; per-obs values override
                    R_per_obs = R_per_obs,
                    t0        = 0.0,
                    verbose   = verbose)
end

# ─────────────────────────────────────────────────────────────────────────────
# Rauch-Tung-Striebel smoother
# ─────────────────────────────────────────────────────────────────────────────

"""
Returned by `run_rts`.

Smoothed estimates evaluated at every observation epoch stored in the
forward `EKFResult`.

- `t`             — observation epochs (length N).
- `y_smooth`      — smoothed full state vector at each epoch (length-NS each).
- `P_smooth`      — smoothed n_solve × n_solve covariance at each epoch.
- `X_hat`         — smoothed full state at the first observation epoch.
- `P_hat`         — smoothed covariance at the first observation epoch.
- `sigma`         — `sqrt.(diag(P_hat))`.
- `solve_for_idx` — indices into each `y_smooth[k]` covered by `P_smooth[k]`.
"""
struct RTSResult
    t             :: Vector{Float64}
    y_smooth      :: Vector{Vector{Float64}}
    P_smooth      :: Vector{Matrix{Float64}}
    X_hat         :: Vector{Float64}
    P_hat         :: Matrix{Float64}
    sigma         :: Vector{Float64}
    solve_for_idx :: Vector{Int}
end

"""
    run_rts(result::EKFResult) -> RTSResult

Standard Rauch-Tung-Striebel backward smoother driven entirely off
the per-step records stored in `result`.  For k = N−1, …, 1:

    C_k        = P_post_k · Φ_{k+1}ᵀ · (P_pre_{k+1})⁻¹
    x_smooth_k = x_filt_k + C_k (x_smooth_{k+1} − x_pred_{k+1})
    P_smooth_k = P_post_k + C_k (P_smooth_{k+1} − P_pre_{k+1}) C_kᵀ

Implementation notes:

  * Operates on the *solve-for sub-state* only.  Non-solve-for slots
    of `y_smooth[k]` are copied verbatim from `records[k].y_post`.
  * `Φ_{k+1}` is read from `records[k+1].Phi`, populated by the
    forward EKF (`time_update!`).
  * The solve is `P_pre_{k+1} \\ M` rather than an explicit inverse,
    matching the Joseph-form style of the forward filter.
"""
function run_rts(result::EKFResult)
    records = result.records
    N       = length(records)
    N >= 1 || throw(ArgumentError(
        "run_rts: the filter result has no records, so there is nothing to smooth. " *
        "Run run_ekf! over at least one observation first"))
    sf      = result.solve_for_idx

    y_smooth = Vector{Vector{Float64}}(undef, N)
    P_smooth = Vector{Matrix{Float64}}(undef, N)

    y_smooth[end] = copy(records[end].y_post)
    P_smooth[end] = copy(records[end].P_post)

    for k in (N - 1):-1:1
        rec_k  = records[k]
        rec_kp = records[k + 1]
        Φ      = rec_kp.Phi
        Cᵀ     = rec_kp.P_pre \ (Φ * rec_k.P_post)
        C      = Cᵀ'

        x_filt_k       = rec_k.y_post[sf]
        x_pred_kp      = rec_kp.y_pre[sf]
        x_smooth_kp_sf = y_smooth[k + 1][sf]
        dx             = x_smooth_kp_sf .- x_pred_kp
        dP             = P_smooth[k + 1] .- rec_kp.P_pre

        x_smooth_k_sf = x_filt_k .+ C * dx
        Pk            = rec_k.P_post .+ C * dP * C'
        P_smooth[k]   = 0.5 .* (Pk .+ Pk')

        y_full = copy(rec_k.y_post)
        y_full[sf] .= x_smooth_k_sf
        y_smooth[k] = y_full
    end

    P_hat = P_smooth[1]
    return RTSResult([records[k].t for k in 1:N],
                     y_smooth, P_smooth,
                     copy(y_smooth[1]), P_hat, sqrt.(diag(P_hat)),
                     copy(sf))
end

# ─────────────────────────────────────────────────────────────────────────────
# Iterated EKF + RTS  (Gauss-Newton-on-the-MAP outer loop)
# ─────────────────────────────────────────────────────────────────────────────

"""
Result of `run_iterated_rts!`.

- `ekf`         — `EKFResult` from the final forward sweep.
- `rts`         — `RTSResult` from the final backward smoothing.
- `iters_run`   — number of outer sweeps actually executed.
- `converged`   — whether the relative-Δ test passed before `max_iters`.
- `delta_history` — per-iteration relative change in the smoothed initial
                    state (length `iters_run`; first entry is `NaN`).
"""
struct IteratedRTSResult
    ekf            :: EKFResult
    rts            :: RTSResult
    iters_run      :: Int
    converged      :: Bool
    delta_history  :: Vector{Float64}
end

# Walk svs in declaration order; copy SolveFor slots out of a full-state
# vector (concatenation in the same declaration order used by `init_ekf`)
# and `assign!` them back onto each SV.  Fixed slots are left alone.
function _assign_solve_for_from_flat!(svs::AbstractVector{<:SolverVariable},
                                       y_flat::AbstractVector{<:Real})
    off = 0
    for s in svs
        n = length_of(s)
        if s.role isa SolveFor
            if n == 1
                assign!(s, y_flat[off + 1])
            else
                assign!(s, collect(y_flat[(off + 1):(off + n)]))
            end
        end
        off += n
    end
    return svs
end

# Build the full initial state vector (all SVs, declaration order) from
# the current values on `svs`.
function _flatten_svs(svs::AbstractVector{<:SolverVariable})
    NS = sum(length_of(s) for s in svs)
    y  = zeros(Float64, NS)
    off = 0
    for s in svs
        n = length_of(s)
        v = current_value(s)
        if n == 1
            y[off + 1] = Float64(v)
        else
            y[(off + 1):(off + n)] .= Float64.(v)
        end
        off += n
    end
    return y
end

# Build the initial-prior covariance over the SolveFor sub-state, in the
# same row/column order `init_ekf` uses.
function _build_initial_P(svs::AbstractVector{<:SolverVariable})
    n_solve = sum(length_of(s) for s in svs if s.role isa SolveFor)
    P0 = zeros(Float64, n_solve, n_solve)
    off = 0
    for s in svs
        s.role isa SolveFor || continue
        n   = length_of(s)
        rng = (off + 1):(off + n)
        a   = s.covariance
        if a isa Real
            for k in rng
                P0[k, k] = Float64(a)
            end
        elseif a isa AbstractVector
            for (kk, k) in enumerate(rng)
                P0[k, k] = Float64(a[kk])
            end
        elseif a isa AbstractMatrix
            P0[rng, rng] = Matrix{Float64}(a)
        else
            throw(ArgumentError(
                "run_iterated_rts!: covariance for $(s.name) must be a Real variance, a " *
                "vector of variances, or a covariance matrix; got $(typeof(a))"))
        end
        off += n
    end
    return P0
end

"""
    run_iterated_rts!(svs, records; model::ODProblem,
                       max_iters::Int = 5,
                       tol::Real = 1e-6,
                       step_size::Real = 1.0,
                       verbose::Bool = false) -> IteratedRTSResult

Iterated forward-EKF + RTS smoother (Gauss-Newton on the MAP problem,
a.k.a. the Bell-Cattivelli iterated Kalman smoother).  At each outer
sweep:

  1. `run_ekf!(svs, records; model)` — forward filter from the current
     SolveFor initial condition (read off `svs` at the EKF reference
     epoch t = 0).
  2. `run_rts(...)` — backward RTS smoother.
  3. The smoothed estimate at the first observation epoch
     (`rts.y_smooth[1]` at `records[1].t > 0`) is back-propagated to
     t = 0 using the STM `records[1].Phi` and the linear RTS gain
     `C0 = P0 Φᵀ P_pre⁻¹`, then written onto the SolveFor SVs as the
     new t = 0 linearization point for the next sweep.

The prior covariance carried on the SVs is **not modified** —
iterating Gauss-Newton with a fixed prior is what gives MAP.

Keyword arguments:

- `max_iters` — maximum outer sweeps (>= 1).
- `tol`       — relative ‖Δx₀‖ convergence threshold on the t = 0
                solve-for slot.
- `step_size` — damping factor in (0, 1] applied to the back-propagated
                correction Δx₀ each iteration.  Use `1.0` for pure
                Gauss-Newton; choose `< 1.0` when the first sweep has
                a large startup transient that would otherwise drive
                the next sweep further from truth (Levenberg-Marquardt
                style).
- `verbose`   — log per-iteration Δ and step norms.

Notes:

- A single iteration (`max_iters=1`) gives the same forward-EKF + RTS
  result as calling `run_ekf!` + `run_rts` directly, plus the back-
  propagated Δx₀ written onto the SVs as a side effect.
- The algorithm helps when the EKF startup transient is dominated by
  a poor t = 0 linearization (initial guess far from truth).  It does
  not by itself address measurement-update linearization error at the
  first few observations; for that, use an iterated *measurement*
  update (IEKF) inside the forward filter.
"""
function run_iterated_rts!(svs::AbstractVector{<:SolverVariable},
                            records::AbstractVector{<:ObservationRecord};
                            model::ODProblem,
                            max_iters::Int = 5,
                            tol::Real = 1e-6,
                            step_size::Real = 1.0,
                            verbose::Bool = false)

    max_iters >= 1 || throw(DomainError(max_iters,
        "run_iterated_rts!: max_iters must be at least 1; got $(max_iters)"))
    (0.0 < step_size <= 1.0) || throw(DomainError(step_size,
        "run_iterated_rts!: step_size must be in (0, 1]; got $(step_size)"))

    P0 = _build_initial_P(svs)

    ekf_result    = nothing
    rts_result    = nothing
    x0_prev       = nothing
    delta_history = Float64[]
    converged     = false
    iter_done     = 0

    for iter in 1:max_iters
        y_init_full = _flatten_svs(svs)

        ekf_result = run_ekf!(svs, records; model = model, verbose = false)
        rts_result = run_rts(ekf_result)
        iter_done  = iter

        sf      = rts_result.solve_for_idx
        rec1    = ekf_result.records[1]
        Φ       = rec1.Phi
        P_pre1  = rec1.P_pre
        y_pre1  = rec1.y_pre
        y_sm1   = rts_result.y_smooth[1]

        dy1_sf  = y_sm1[sf] .- y_pre1[sf]
        C0      = (P0 * Φ') / P_pre1
        dy0_sf  = step_size .* (C0 * dy1_sf)

        y_at_t0       = copy(y_init_full)
        y_at_t0[sf]  .= y_init_full[sf] .+ dy0_sf
        _assign_solve_for_from_flat!(svs, y_at_t0)

        x0_new = y_at_t0[sf]
        if x0_prev === nothing
            push!(delta_history, NaN)
            verbose && @info "run_iterated_rts! iter=$iter (seed)" Δx0_step = norm(dy0_sf)
        else
            Δ   = norm(x0_new .- x0_prev)
            nx  = max(norm(x0_new), 1.0)
            rel = Δ / nx
            push!(delta_history, rel)
            verbose && @info "run_iterated_rts! iter=$iter" Δ rel Δx0_step = norm(dy0_sf)
            if rel < tol
                converged = true
                break
            end
        end
        x0_prev = x0_new
    end

    return IteratedRTSResult(ekf_result, rts_result,
                              iter_done, converged, delta_history)
end

end # module ExtendedKalmanFilter
