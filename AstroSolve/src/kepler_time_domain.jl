# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Analytic two-body propagation over a given time, with the partials the
# shooting transcriptions need. Included from AstroSolve.jl in order; not
# standalone, it assumes what the files before it define.

"""
    KeplerPropagationResult

The propagated state, derivatives, coefficients, and convergence status from time-domain Kepler propagation.

# Fields
- `state::Vector{T}`: Final Cartesian state `[r; v]`.
- `stm::Union{Nothing,Matrix{T}}`: State transition matrix `∂x_f/∂x_0`, when requested.
- `dstate_dt::Union{Nothing,Vector{T}}`: Time derivative `∂x_f/∂t`, when requested.
- `F::T`, `G::T`, `Ft::T`, `Gt::T`: Lagrange coefficients.
- `iterations::Int`: Total universal-variable solver iterations.
- `converged::Bool`: Whether the universal-variable solve met tolerance.

# Example
```julia
using AstroSolve
result = AstroSolve.kepler_propagate_time_domain(
    [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0], 398600.4418, 60.0)
```
"""
struct KeplerPropagationResult{T<:Real}
    state::Vector{T}
    stm::Union{Nothing, Matrix{T}}
    dstate_dt::Union{Nothing, Vector{T}}
    F::T
    G::T
    Ft::T
    Gt::T
    iterations::Int
    converged::Bool
end

"""
    kepler_propagate_time_domain(state0, mu, dt; need_stm=true, need_time_partials=true,
                                 alphatol=1e-12, Xtol=1e-12,
                                 maximum_order=10, maximum_iterations_per_order=10)

Propagate a Cartesian two-body state through a specified elapsed time, with a
universal-variable Kepler solver following the EMTG `KeplerPropagatorTimeDomain` approach.

# Arguments
- `state0::AbstractVector{<:Real}`: Initial Cartesian state `[x,y,z,vx,vy,vz]`.
- `mu::Real`: Positive gravitational parameter.
- `dt::Real`: Elapsed time; positive propagates forward and negative propagates backward.
- `need_stm::Bool`: Whether to return the state transition matrix.
- `need_time_partials::Bool`: Whether to return the state derivative with respect to elapsed time.
- `alphatol::Float64`: Tolerance used to distinguish elliptic, parabolic, and hyperbolic motion.
- `Xtol::Float64`: Universal-anomaly convergence tolerance.
- `maximum_order::Int`: Highest universal-variable series order attempted.
- `maximum_iterations_per_order::Int`: Iteration limit at each series order.

# Notes
Position, velocity, `mu`, and `dt` must use consistent units. The state is
returned in the same frame as `state0`. Throws an `ArgumentError` when `state0`
does not contain six elements, `mu` is not positive, or the initial radius is zero.

# Returns
`KeplerPropagationResult` containing propagated state and requested partial derivatives.

# Example
```julia
using AstroSolve
result = AstroSolve.kepler_propagate_time_domain(
    [7000.0, 0.0, 0.0, 0.0, 7.546, 0.0], 398600.4418, 60.0)
result.state
```
"""
function kepler_propagate_time_domain(
    state0::AbstractVector{<:Real},
    mu::Real,
    dt::Real;
    need_stm::Bool = true,
    need_time_partials::Bool = true,
    alphatol::Float64 = 1.0e-12,
    Xtol::Float64 = 1.0e-12,
    maximum_order::Int = 10,
    maximum_iterations_per_order::Int = 10,
)
    length(state0) == 6 || throw(ArgumentError(
        "state0 must be a length-6 Cartesian state [r; v]; got length $(length(state0))"))
    mu > 0 || throw(ArgumentError(
        "gravitational parameter mu must be positive; got mu = $(mu)"))

    T = promote_type(eltype(state0), typeof(mu), typeof(dt))
    x0 = convert(Vector{T}, state0)
    r0vec = @view x0[1:3]
    v0vec = @view x0[4:6]
    r0 = norm(r0vec)
    v0 = norm(v0vec)
    r0 > 0 || throw(ArgumentError(
        "the initial position must have non-zero magnitude, since the orbit is undefined at " *
        "the centre; got |r| = $(r0)"))

    sqmu = sqrt(T(mu))
    alpha = 2 / r0 - v0 * v0 / T(mu)
    sigma0 = dot(r0vec, v0vec) / sqmu

    sqalpha  = alpha > alphatol ? sqrt(alpha)   : zero(T)
    sqmalpha = alpha < -alphatol ? sqrt(-alpha) : zero(T)

    X_new = alpha > alphatol ? alpha * sqmu * T(dt) : T(0.1) * sqmu * T(dt) / r0
    X  = T(1.0e100)
    dX = T(1.0e100)
    r  = r0
    sigma = zero(T)
    U0 = zero(T)
    U1 = zero(T)
    U2 = zero(T)
    U3 = zero(T)

    N = 2
    iteration_this_N = 0
    total_iterations = 0
    converged = false

    # Laguerre-Conway-Der iteration on Kepler's universal equation (Battin §4.5)
    while abs(X - X_new) > Xtol && N < maximum_order
        iteration_this_N += 1
        total_iterations += 1
        if iteration_this_N >= maximum_iterations_per_order
            N += 1
            iteration_this_N = 0
        end

        X = X_new

        if alpha > alphatol
            y = alpha * X * X
            sqrt_y = sqrt(y)
            C = (1.0 - cos(sqrt_y)) / y
            S = (sqrt_y - sin(sqrt_y)) / sqrt(y * y * y)
            U1 = X * (1.0 - y * S)
            U2 = X * X * C
            U3 = X * X * X * S
            U0 = 1.0 - alpha * U2
        elseif alpha < -alphatol
            sqmalphaX = sqmalpha * X
            if sqmalphaX > 30.0 || sqmalphaX < -30.0
                if N < maximum_order
                    N += 1
                    iteration_this_N = 0
                    continue
                else
                    break
                end
            end
            U0 = cosh(sqmalphaX)
            U1 = sinh(sqmalphaX) / sqmalpha
            U2 = (1.0 - U0) / alpha
            U3 = (X - U1) / alpha
        else  # parabola: Taylor limits of Stumpff functions as α→0
            U0 = 1.0
            U1 = X
            U2 = U1 * X / 2.0
            U3 = U2 * X / 3.0
        end

        r = r0 * U0 + sigma0 * U1 + U2
        sigma = sigma0 * U0 + (1.0 - alpha * r0) * U1

        FX = r0 * U1 + sigma0 * U2 + U3 - sqmu * T(dt)
        dFX = r
        ddFX = sigma

        sgn = dFX >= 0 ? 1.0 : -1.0
        denom = abs((N - 1) * (N - 1) * dFX * dFX - N * (N - 1) * FX * ddFX)
        if denom > 0.0
            dX = N * FX / (dFX + sgn * sqrt(denom))
        else
            dX = FX / dFX
        end

        X_new -= dX
    end

    converged = abs(X - X_new) <= Xtol

    # Lagrange coefficients F, G, Ḟ, Ġ (Battin §4.4)
    F  = 1.0 - U2 / r0
    G  = (r0 * U1 + sigma0 * U2) / sqmu
    Ft = -sqmu / (r0 * r) * U1
    Gt = 1.0 - U2 / r

    xf = Vector{T}(undef, 6)
    xf[1] = F * x0[1] + G * x0[4]
    xf[2] = F * x0[2] + G * x0[5]
    xf[3] = F * x0[3] + G * x0[6]
    xf[4] = Ft * x0[1] + Gt * x0[4]
    xf[5] = Ft * x0[2] + Gt * x0[5]
    xf[6] = Ft * x0[3] + Gt * x0[6]

    stm = nothing
    dstate_dt = nothing

    if need_stm || need_time_partials
        fr = r
        fr0 = r0

        # Parabolic limit: U_n(Ψ=0) = X^n / n!  (Stumpff c_n(0) = 1/n!)
        # The a*(X²/2 - U2) form gives 1e30*0 = 0 when alpha=0, which is wrong.
        if abs(alpha) < alphatol
            U4 = X * X * X * X / T(24)
            U5 = X * X * X * X * X / T(120)
        else
            a = one(T) / alpha
            U4 = a * (X * X / T(2) - U2)
            U5 = a * (X * X * X / T(6) - U3)
        end

        dXdt = sqmu / fr
        U0dot = -alpha * U1 * dXdt
        U1dot = U0 * dXdt
        U2dot = U1 * dXdt

        Cb = (1 / sqmu) * (3 * U5 - X * U4) - T(dt) * U2

        x00, y00, z00 = x0[1], x0[2], x0[3]
        xdot0, ydot0, zdot0 = x0[4], x0[5], x0[6]
        x, y, z = xf[1], xf[2], xf[3]
        xdot, ydot, zdot = xf[4], xf[5], xf[6]

        fr2 = fr * fr
        fr02 = fr0 * fr0
        fr3 = fr2 * fr
        fr03 = fr02 * fr0
        mu2 = T(mu) * T(mu)

        if need_stm
            stm = Matrix{T}(undef, 6, 6)

            stm[1, 1] = fr / mu * (xdot - xdot0)^2 + (fr0 * (1.0 - F) * (x * x00) + Cb * (x00 * xdot)) / fr03 + F
            stm[1, 2] = fr / mu * (xdot - xdot0) * (ydot - ydot0) + (fr0 * (1.0 - F) * (x * y00) + Cb * (y00 * xdot)) / fr03
            stm[1, 3] = fr / mu * (xdot - xdot0) * (zdot - zdot0) + (fr0 * (1.0 - F) * (x * z00) + Cb * (z00 * xdot)) / fr03
            stm[2, 1] = fr / mu * (ydot - ydot0) * (xdot - xdot0) + (fr0 * (1.0 - F) * (y * x00) + Cb * (x00 * ydot)) / fr03
            stm[2, 2] = fr / mu * (ydot - ydot0)^2 + (fr0 * (1.0 - F) * (y * y00) + Cb * (y00 * ydot)) / fr03 + F
            stm[2, 3] = fr / mu * (ydot - ydot0) * (zdot - zdot0) + (fr0 * (1.0 - F) * (y * z00) + Cb * (z00 * ydot)) / fr03
            stm[3, 1] = fr / mu * (zdot - zdot0) * (xdot - xdot0) + (fr0 * (1.0 - F) * (z * x00) + Cb * (x00 * zdot)) / fr03
            stm[3, 2] = fr / mu * (zdot - zdot0) * (ydot - ydot0) + (fr0 * (1.0 - F) * (z * y00) + Cb * (y00 * zdot)) / fr03
            stm[3, 3] = fr / mu * (zdot - zdot0)^2 + (fr0 * (1.0 - F) * (z * z00) + Cb * (z00 * zdot)) / fr03 + F

            stm[1, 4] = fr0 / mu * (1.0 - F) * (xdot0 * (x - x00) - x00 * (xdot - xdot0)) + Cb / mu * (xdot * xdot0) + G
            stm[1, 5] = fr0 / mu * (1.0 - F) * (ydot0 * (x - x00) - y00 * (xdot - xdot0)) + Cb / mu * (xdot * ydot0)
            stm[1, 6] = fr0 / mu * (1.0 - F) * (zdot0 * (x - x00) - z00 * (xdot - xdot0)) + Cb / mu * (xdot * zdot0)
            stm[2, 4] = fr0 / mu * (1.0 - F) * (xdot0 * (y - y00) - x00 * (ydot - ydot0)) + Cb / mu * (ydot * xdot0)
            stm[2, 5] = fr0 / mu * (1.0 - F) * (ydot0 * (y - y00) - y00 * (ydot - ydot0)) + Cb / mu * (ydot * ydot0) + G
            stm[2, 6] = fr0 / mu * (1.0 - F) * (zdot0 * (y - y00) - z00 * (ydot - ydot0)) + Cb / mu * (ydot * zdot0)
            stm[3, 4] = fr0 / mu * (1.0 - F) * (xdot0 * (z - z00) - x00 * (zdot - zdot0)) + Cb / mu * (zdot * xdot0)
            stm[3, 5] = fr0 / mu * (1.0 - F) * (ydot0 * (z - z00) - y00 * (zdot - zdot0)) + Cb / mu * (zdot * ydot0)
            stm[3, 6] = fr0 / mu * (1.0 - F) * (zdot0 * (z - z00) - z00 * (zdot - zdot0)) + Cb / mu * (zdot * zdot0) + G

            stm[4, 1] = (-Cb * mu2 * x * x00 + Ft * fr * fr03 * (mu * fr2 - mu * x * x + fr * (xdot0 - xdot) * (y * (xdot * y - ydot * x) + z * (xdot * z - zdot * x))) + mu * fr3 * fr0 * x00 * (xdot0 - xdot) + mu * fr * fr03 * x * (xdot0 - xdot)) / (mu * fr3 * fr03)
            stm[4, 2] = (-Cb * mu2 * x * y00 + Ft * fr * fr03 * (-mu * x * y + fr * (ydot0 - ydot) * (y * (xdot * y - ydot * x) + z * (xdot * z - zdot * x))) + mu * fr3 * fr0 * y00 * (xdot0 - xdot) + mu * fr * fr03 * x * (ydot0 - ydot)) / (mu * fr3 * fr03)
            stm[4, 3] = (-Cb * mu2 * x * z00 + Ft * fr * fr03 * (-mu * x * z + fr * (zdot0 - zdot) * (y * (xdot * y - ydot * x) + z * (xdot * z - zdot * x))) + mu * fr3 * fr0 * z00 * (xdot0 - xdot) + mu * fr * fr03 * x * (zdot0 - zdot)) / (mu * fr3 * fr03)
            stm[5, 1] = (-Cb * mu2 * x00 * y - Ft * fr * fr03 * (mu * x * y + fr * (xdot0 - xdot) * (x * (xdot * y - ydot * x) - z * (ydot * z - zdot * y))) + mu * fr3 * fr0 * x00 * (ydot0 - ydot) + mu * fr * fr03 * y * (xdot0 - xdot)) / (mu * fr3 * fr03)
            stm[5, 2] = (-Cb * mu2 * y * y00 - Ft * fr * fr03 * (-mu * fr2 + mu * y * y + fr * (ydot0 - ydot) * (x * (xdot * y - ydot * x) - z * (ydot * z - zdot * y))) + mu * fr3 * fr0 * y00 * (ydot0 - ydot) + mu * fr * fr03 * y * (ydot0 - ydot)) / (mu * fr3 * fr03)
            stm[5, 3] = (-Cb * mu2 * y * z00 - Ft * fr * fr03 * (mu * y * z + fr * (zdot0 - zdot) * (x * (xdot * y - ydot * x) - z * (ydot * z - zdot * y))) + mu * fr3 * fr0 * z00 * (ydot0 - ydot) + mu * fr * fr03 * y * (zdot0 - zdot)) / (mu * fr3 * fr03)
            stm[6, 1] = (-Cb * mu2 * x00 * z - Ft * fr * fr03 * (mu * x * z + fr * (xdot0 - xdot) * (x * (xdot * z - zdot * x) + y * (ydot * z - zdot * y))) + mu * fr3 * fr0 * x00 * (zdot0 - zdot) + mu * fr * fr03 * z * (xdot0 - xdot)) / (mu * fr3 * fr03)
            stm[6, 2] = (-Cb * mu2 * y00 * z - Ft * fr * fr03 * (mu * y * z + fr * (ydot0 - ydot) * (x * (xdot * z - zdot * x) + y * (ydot * z - zdot * y))) + mu * fr3 * fr0 * y00 * (zdot0 - zdot) + mu * fr * fr03 * z * (ydot0 - ydot)) / (mu * fr3 * fr03)
            stm[6, 3] = (-Cb * mu2 * z * z00 - Ft * fr * fr03 * (-mu * fr2 + mu * z * z + fr * (zdot0 - zdot) * (x * (xdot * z - zdot * x) + y * (ydot * z - zdot * y))) + mu * fr3 * fr0 * z00 * (zdot0 - zdot) + mu * fr * fr03 * z * (zdot0 - zdot)) / (mu * fr3 * fr03)

            stm[4, 4] = fr0 / mu * (xdot - xdot0)^2 + (fr0 * (1.0 - F) * (x * x00) - Cb * (x * xdot0)) / fr3 + Gt
            stm[4, 5] = fr0 / mu * (xdot - xdot0) * (ydot - ydot0) + (fr0 * (1.0 - F) * (x * y00) - Cb * (x * ydot0)) / fr3
            stm[4, 6] = fr0 / mu * (xdot - xdot0) * (zdot - zdot0) + (fr0 * (1.0 - F) * (x * z00) - Cb * (x * zdot0)) / fr3
            stm[5, 4] = fr0 / mu * (ydot - ydot0) * (xdot - xdot0) + (fr0 * (1.0 - F) * (y * x00) - Cb * (y * xdot0)) / fr3
            stm[5, 5] = fr0 / mu * (ydot - ydot0)^2 + (fr0 * (1.0 - F) * (y * y00) - Cb * (y * ydot0)) / fr3 + Gt
            stm[5, 6] = fr0 / mu * (ydot - ydot0) * (zdot - zdot0) + (fr0 * (1.0 - F) * (y * z00) - Cb * (y * zdot0)) / fr3
            stm[6, 4] = fr0 / mu * (zdot - zdot0) * (xdot - xdot0) + (fr0 * (1.0 - F) * (z * x00) - Cb * (z * xdot0)) / fr3
            stm[6, 5] = fr0 / mu * (zdot - zdot0) * (ydot - ydot0) + (fr0 * (1.0 - F) * (z * y00) - Cb * (z * ydot0)) / fr3
            stm[6, 6] = fr0 / mu * (zdot - zdot0)^2 + (fr0 * (1.0 - F) * (z * z00) - Cb * (z * zdot0)) / fr3 + Gt
        end

        if need_time_partials
            rdot = fr0 * U0dot + sigma0 * U1dot + U2dot
            r2 = fr * fr
            Ftt = -sqmu / fr0 * (U1dot / fr - U1 * rdot / r2)
            Gtt = -(U2dot / fr - U2 * rdot / r2)

            dstate_dt = Vector{T}(undef, 6)
            dstate_dt[1] = Ft * x0[1] + Gt * x0[4]
            dstate_dt[2] = Ft * x0[2] + Gt * x0[5]
            dstate_dt[3] = Ft * x0[3] + Gt * x0[6]
            dstate_dt[4] = Ftt * x0[1] + Gtt * x0[4]
            dstate_dt[5] = Ftt * x0[2] + Gtt * x0[5]
            dstate_dt[6] = Ftt * x0[3] + Gtt * x0[6]
        end
    end

    return KeplerPropagationResult(xf, stm, dstate_dt, F, G, Ft, Gt, total_iterations, converged)
end
