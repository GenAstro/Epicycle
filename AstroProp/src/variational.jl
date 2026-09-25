# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Propagating a state alongside its state transition matrix.
#
# The variational equation is Phi-dot = A(t) Phi, so an arc that carries Phi
# needs the Jacobian of its own dynamics. That Jacobian is an argument here
# rather than something looked up: a caller that declared an analytic partial
# passes it, a caller that did not passes an automatic-differentiation closure,
# and this file needs to know about neither. Which is what lets it live in the
# propagation package instead of the solver.


# Single-phase augmented RHS: Y = [y; vec(Φ_y)]
function _augmented_rhs!(dY, Y, ctx, t)
    n   = ctx.n
    y_v = @view Y[1:n]
    Φ   = reshape(@view(Y[n+1:n+n*n]), n, n)

    f_val = ctx.f(y_v, ctx.u, ctx.p, t, ctx.model)
    @inbounds for i in 1:n
        dY[i] = f_val[i]
    end

    A = ctx.A_of(y_v, ctx.u, ctx.p, t, ctx.model)
    dΦ = A * Φ
    @inbounds for j in 1:n, i in 1:n
        dY[n + i + (j - 1) * n] = dΦ[i, j]
    end
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Single-phase API
# ─────────────────────────────────────────────────────────────────────────────

"""
    propagate_with_stm(f, dfdy, y0, t0, t1, u, p, model;
                       solver = Tsit5(), reltol = 1e-12, abstol = 1e-12)
        -> (y_final::Vector, Φ::Matrix)

Integrate a state and its state transition matrix Φ(t, t0) from `t0` to `t1`, solving
ẏ = f and Φ̇ = A(t) Φ with Φ(t0) = I.

# Arguments
- `f`: Dynamics, called as `f(y, u, p, t, model)` and returning ẏ.
- `dfdy`: The Jacobian A = ∂f/∂y, called with the same arguments. Pass an analytic partial or an
  automatic-differentiation closure.
- `y0`: Initial state at `t0`.
- `t0`, `t1`: Start and end times, in the time units `f` uses. `t1 < t0` integrates backward.
- `u`, `p`, `model`: Passed through to `f` and `dfdy` unchanged.

# Returns
The state at `t1` and the n×n state transition matrix Φ(t1, t0).
"""
function propagate_with_stm(f, dfdy,
                            y0::AbstractVector,
                            t0::Real, t1::Real,
                            u, p, model;
                            solver = Tsit5(),
                            reltol = 1e-12, abstol = 1e-12)
    n  = length(y0)
    Y0 = zeros(Float64, n + n * n)
    @inbounds for i in 1:n
        Y0[i] = Float64(y0[i])
    end
    @inbounds for i in 1:n
        Y0[n + i + (i - 1) * n] = 1.0   # vec(I)
    end

    A_of = dfdy
    ctx  = (f = f, A_of = A_of, n = n, u = u, p = p, model = model)

    if t1 == t0
        return Y0[1:n], Matrix{Float64}(I, n, n)
    end

    prob = ODEProblem(_augmented_rhs!, Y0, (Float64(t0), Float64(t1)), ctx)
    sol  = solve(prob, solver; reltol = reltol, abstol = abstol,
                 save_everystep = false, save_start = false)

    Yf = sol.u[end]
    return Yf[1:n], reshape(Yf[n+1:n+n*n], n, n)
end

"""
    propagate_with_stm_to_times(f, dfdy, y0, t0, ts, u, p, model;
                                solver = Tsit5(),
                                reltol = 1e-12, abstol = 1e-12)
        -> (ys::Vector{Vector{Float64}}, Φs::Vector{Matrix{Float64}})

Integrate a state and its state transition matrix once from `t0` and sample both at each time in
`ts`. Arguments are those of `propagate_with_stm`.

# Returns
One state and one Φ(tₖ, t0) per entry of `ts`, in the order given. Throws `ArgumentError` when a
time in `ts` precedes `t0`.
"""
function propagate_with_stm_to_times(f, dfdy,
                                      y0::AbstractVector,
                                      t0::Real,
                                      ts::AbstractVector{<:Real},
                                      u, p, model;
                                      solver = Tsit5(),
                                      reltol = 1e-12, abstol = 1e-12)
    n  = length(y0)
    Y0 = zeros(Float64, n + n * n)
    @inbounds for i in 1:n
        Y0[i] = Float64(y0[i])
    end
    @inbounds for i in 1:n
        Y0[n + i + (i - 1) * n] = 1.0
    end

    A_of = dfdy
    ctx  = (f = f, A_of = A_of, n = n, u = u, p = p, model = model)

    t0f  = Float64(t0)
    tsf  = collect(Float64.(ts))
    # The solve runs forward from t0, so an earlier sample would be read off the
    # interpolant outside its span, which extrapolates rather than failing.
    all(>=(t0f), tsf) || throw(ArgumentError(
        "propagate_with_stm_to_times: every time in ts must be at or after t0 = $(t0f); " *
        "got $(minimum(tsf))"))
    tend = max(t0f, maximum(tsf))

    prob = ODEProblem(_augmented_rhs!, Y0, (t0f, tend), ctx)
    sol  = solve(prob, solver; reltol = reltol, abstol = abstol)

    ys = Vector{Vector{Float64}}(undef, length(tsf))
    Φs = Vector{Matrix{Float64}}(undef, length(tsf))
    @inbounds for (k, t) in enumerate(tsf)
        if t == t0f
            ys[k] = Y0[1:n]
            Φs[k] = Matrix{Float64}(I, n, n)
        else
            Y = sol(t)
            ys[k] = Y[1:n]
            Φs[k] = reshape(Y[n+1:n+n*n], n, n)
        end
    end
    return ys, Φs
end

# ─────────────────────────────────────────────────────────────────────────────
# Propagation with parameter sensitivity Φ_p = ∂y/∂p
# ─────────────────────────────────────────────────────────────────────────────

# Augmented RHS: Y = [y; vec(Φ_y); vec(Φ_p)]
function _augmented_rhs_sens!(dY, Y, ctx, t)
    n   = ctx.n
    n_p = ctx.n_p

    y_v = @view Y[1:n]
    Φy  = reshape(@view(Y[n+1:n+n*n]), n, n)
    Φp  = reshape(@view(Y[n+n*n+1:n+n*n+n*n_p]), n, n_p)

    f_val = ctx.f(y_v, ctx.u, ctx.p, t, ctx.model)
    @inbounds for i in 1:n
        dY[i] = f_val[i]
    end

    A = ctx.A_of(y_v, ctx.u, ctx.p, t, ctx.model)
    dΦy = A * Φy
    @inbounds for j in 1:n, i in 1:n
        dY[n + i + (j - 1) * n] = dΦy[i, j]
    end

    if n_p > 0
        B = ctx.B_of(y_v, ctx.u, ctx.p, t, ctx.model)
        dΦp = A * Φp .+ B
        base = n + n * n
        @inbounds for j in 1:n_p, i in 1:n
            dY[base + i + (j - 1) * n] = dΦp[i, j]
        end
    end
    return nothing
end

"""
    propagate_with_sensitivities_to_times(f, dfdy, dfdp, y0, t0, ts, u, p, model;
                                          solver = Tsit5(),
                                          reltol = 1e-12, abstol = 1e-12)
        -> (ys, Φys, Φps)

Integrate a state, its state transition matrix Φ(t, t0), and its sensitivity to the parameters
∂y/∂p once from `t0`, and sample all three at each time in `ts`. Solves ẏ = f,
Φ̇ = A Φ with Φ(t0) = I, and Φ̇p = A Φp + B with Φp(t0) = 0.

# Arguments
- `f`, `dfdy`, `y0`, `t0`, `u`, `model`: as for `propagate_with_stm`.
- `dfdp`: The Jacobian B = ∂f/∂p, n × length(`p`), called as `dfdp(y, u, p, t, model)`.
- `ts`: Sample times, each at or after `t0`.
- `p`: The parameters, a `NamedTuple`; `(;)` for none.

# Returns
One state, one n × n Φ(tₖ, t0) and one n × length(`p`) ∂y/∂p per entry of `ts`, in the order
given. With no parameters each ∂y/∂p is n × 0. Throws `ArgumentError` when a time in `ts`
precedes `t0`.
"""
function propagate_with_sensitivities_to_times(f, dfdy, dfdp,
                                               y0::AbstractVector,
                                               t0::Real,
                                               ts::AbstractVector{<:Real},
                                               u,
                                               p::NamedTuple,
                                               model;
                                               solver = Tsit5(),
                                               reltol = 1e-12,
                                               abstol = 1e-12)
    n   = length(y0)
    n_p = length(p)

    Y0 = zeros(Float64, n + n * n + n * n_p)
    @inbounds for i in 1:n
        Y0[i] = Float64(y0[i])
    end
    @inbounds for i in 1:n
        Y0[n + i + (i - 1) * n] = 1.0
    end

    t0f = Float64(t0)
    tsf = collect(Float64.(ts))
    # The solve runs forward from t0, so an earlier sample would be read off the interpolant
    # outside its span, which extrapolates rather than failing.
    all(>=(t0f), tsf) || throw(ArgumentError(
        "propagate_with_sensitivities_to_times: every time in ts must be at or after " *
        "t0 = $(t0f); got $(minimum(tsf))"))
    tend = max(t0f, maximum(tsf))

    ctx  = (f = f, A_of = dfdy, B_of = dfdp, n = n, n_p = n_p, u = u, p = p, model = model)
    prob = ODEProblem(_augmented_rhs_sens!, Y0, (t0f, tend), ctx)
    sol  = solve(prob, solver; reltol = reltol, abstol = abstol)

    ys  = Vector{Vector{Float64}}(undef, length(tsf))
    Φys = Vector{Matrix{Float64}}(undef, length(tsf))
    Φps = Vector{Matrix{Float64}}(undef, length(tsf))
    @inbounds for (k, t) in enumerate(tsf)
        Y = t == t0f ? Y0 : sol(t)
        ys[k]  = Y[1:n]
        Φys[k] = reshape(Y[n+1:n+n*n], n, n)
        Φps[k] = reshape(Y[n+n*n+1:n+n*n+n*n_p], n, n_p)
    end
    return ys, Φys, Φps
end
