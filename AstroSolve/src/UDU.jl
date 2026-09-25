# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
# UDU.jl
#
# UDU' covariance factorization (Bierman / Thornton).
#
#   P = U · diag(D) · Uᵀ,    U unit upper triangular,  D diagonal (vector).
#
# Implements:
#   · `udu_from_P(P)`             — factor a SPD matrix into (U, D).
#   · `udu_to_P(U, D)`            — reconstruct P (for testing / reporting).
#   · `thornton_time_update!`     — TSB §5.7.2:  P⁺ = Φ P Φᵀ + G Q Gᵀ.
#   · `bierman_measurement_update!` — TSB §5.7.3: scalar-by-scalar update.
#
# Step-1 scope: nothing fancy.  No square-root form per se beyond UDU; no
# robust factorization tricks.  No process noise plumbing in the use case
# (Q is empty), but the Thornton routine accepts it so we don't have to
# revisit when DMC arrives.

module UDU

using LinearAlgebra

export udu_from_P, udu_to_P,
       thornton_time_update!, bierman_measurement_update!

"""
    udu_from_P(P) -> (U, D)

Factor a symmetric positive-definite matrix into unit upper triangular and diagonal parts.

# Arguments
- `P`: Symmetric positive-definite matrix, `n × n`. Only the upper triangle is read.

# Returns
`(U, D)` where `U` is `n × n` unit upper triangular and `D` is the length-`n` vector of
diagonal entries, satisfying `P == U * Diagonal(D) * U'`.

# Notes
Carrying a covariance as `U` and `D` rather than as `P` keeps it positive definite through
arithmetic that would otherwise lose that property, which is the reason the estimators use
this form throughout.

`P` is copied, not modified. Throws `DomainError` on the first non-positive pivot, naming
which one, since that means the input was not positive definite.

Bierman's outer-product algorithm; see Gelb, *Applied Optimal Estimation*, and Grewal &
Andrews, *Kalman Filtering* (3e) §6.4.

# Example
```julia
using AstroSolve.UDU: udu_from_P, udu_to_P

P = [4.0 2.0 0.5; 2.0 5.0 1.0; 0.5 1.0 3.0]
U, D = udu_from_P(P)
udu_to_P(U, D)        # returns P
```
"""
function udu_from_P(P::AbstractMatrix)
    n   = size(P, 1)
    A   = Matrix{Float64}(P)        # working copy (destroyed)
    U   = Matrix{Float64}(I, n, n)
    D   = zeros(Float64, n)

    @inbounds for j in n:-1:2
        D[j] = A[j, j]
        D[j] > 0 || throw(DomainError(D[j],
            "udu_from_P: P must be positive definite; pivot $(j) of $(n) is $(D[j]), " *
            "which is not positive"))
        @inbounds for k in 1:(j - 1)
            U[k, j] = A[k, j] / D[j]
        end
        @inbounds for k in 1:(j - 1)
            @inbounds for i in 1:k
                A[i, k] -= U[i, j] * D[j] * U[k, j]
            end
        end
    end
    D[1] = A[1, 1]
    D[1] > 0 || throw(DomainError(D[1],
        "udu_from_P: P must be positive definite; pivot 1 of $(n) is $(D[1]), " *
        "which is not positive"))
    return U, D
end

"""
    udu_to_P(U, D) -> Matrix

Reconstruct the covariance from its UDU factors.

# Arguments
- `U`: Unit upper triangular matrix, `n × n`.
- `D`: Diagonal entries, length `n`.

# Returns
The `n × n` matrix `U * Diagonal(D) * U'`.

# Notes
The inverse of [`udu_from_P`](@ref), and the round trip is how a factored covariance is
checked or reported. The estimators do not call it in their inner loops.

# Example
```julia
using AstroSolve.UDU: udu_from_P, udu_to_P

U, D = udu_from_P([4.0 2.0; 2.0 5.0])
udu_to_P(U, D)
```
"""
udu_to_P(U::AbstractMatrix, D::AbstractVector) = U * Diagonal(D) * U'

"""
    thornton_time_update!(U, D, Φ, G, Q) -> nothing

Propagate a factored covariance forward one step, in place.

# Arguments
- `U`: Unit upper triangular factor, `n × n`. Overwritten.
- `D`: Diagonal factor, length `n`. Overwritten.
- `Φ`: State transition matrix over the step, `n × n`.
- `G`: Process noise mapping, `n × nq`. May be `n × 0`.
- `Q`: Process noise variances, length `nq`. May be empty.

# Returns
`nothing`. `U` and `D` are updated in place to the factors of
`P⁺ = Φ (U D Uᵀ) Φᵀ + G diag(Q) Gᵀ`.

# Notes
The propagated covariance is never formed. Modified weighted Gram-Schmidt runs over the
augmented matrix `[Φ U | G]` with weights `[D; Q]`, which is what keeps the result positive
definite where forming `Φ P Φᵀ` and refactoring would not.

An empty `G` and `Q` mean deterministic dynamics and are the normal case until process noise
is modelled.

Throws `ArgumentError` if any argument's shape disagrees with `D`, and `DomainError` if the
updated covariance is not positive definite, which a non-positive entry of `D` or `Q` on
input will cause.

Thornton's algorithm; see Grewal & Andrews, *Kalman Filtering* (3e) §6.5.

# Example
```julia
using AstroSolve.UDU: udu_from_P, udu_to_P, thornton_time_update!

U, D = udu_from_P([4.0 0.2; 0.2 1.0])
Φ = [1.0 0.5; 0.0 1.0]
thornton_time_update!(U, D, Φ, zeros(2, 0), Float64[])
udu_to_P(U, D)        # Φ P Φ'
```
"""
function thornton_time_update!(U::AbstractMatrix, D::AbstractVector,
                                Φ::AbstractMatrix,
                                G::AbstractMatrix,
                                Q::AbstractVector)
    n  = length(D)
    nq = length(Q)
    size(U) == (n, n) || throw(ArgumentError(
        "thornton_time_update!: U must be ($(n), $(n)) to match D; got $(size(U))"))
    size(Φ) == (n, n) || throw(ArgumentError(
        "thornton_time_update!: Φ must be ($(n), $(n)) to match D; got $(size(Φ))"))
    size(G, 1) == n && size(G, 2) == nq ||
        throw(ArgumentError(
            "thornton_time_update!: G must be ($(n), $(nq)) to match D and Q; " *
            "got $(size(G))"))

    V    = Φ * U                              # n × n
    W    = hcat(V, Matrix{Float64}(G))        # n × (n + nq)
    Dbar = vcat(Float64.(D), Float64.(Q))     # n + nq

    Unew = Matrix{Float64}(I, n, n)
    Dnew = zeros(Float64, n)
    nw   = n + nq

    @inbounds for i in n:-1:1
        σ = 0.0
        @inbounds for k in 1:nw
            σ += W[i, k]^2 * Dbar[k]
        end
        σ > 0 || throw(DomainError(σ,
            "thornton_time_update!: the updated covariance must stay positive " *
            "definite; row $(i) of $(n) gives σ = $(σ). A non-positive D or Q entry " *
            "on input is the usual cause"))
        Dnew[i] = σ
        @inbounds for j in 1:(i - 1)
            s = 0.0
            @inbounds for k in 1:nw
                s += W[j, k] * Dbar[k] * W[i, k]
            end
            uji = s / σ
            Unew[j, i] = uji
            @inbounds for k in 1:nw
                W[j, k] -= uji * W[i, k]
            end
        end
    end

    @inbounds for j in 1:n, i in 1:n
        U[i, j] = Unew[i, j]
    end
    @inbounds for i in 1:n
        D[i] = Dnew[i]
    end
    return nothing
end

"""
    bierman_measurement_update!(U, D, x, h, R, ν) -> Vector

Apply one scalar measurement to a factored covariance and state, in place.

# Arguments
- `U`: Unit upper triangular factor, `n × n`. Overwritten.
- `D`: Diagonal factor, length `n`. Overwritten.
- `x`: State estimate, length `n`. Overwritten.
- `h`: Measurement partials `∂z/∂x`, length `n`.
- `R`: Measurement variance, in the square of the measurement's units. Must be positive.
- `ν`: Innovation `z − ẑ`, the observed measurement minus the one predicted at `x`.

# Returns
The Kalman gain, length `n`. `U`, `D` and `x` are updated in place.

# Notes
One measurement at a time, so a vector observation is applied by calling this once per
component. That is the point of the algorithm rather than a limitation: each scalar update
keeps the covariance factored and positive definite, where a batched matrix inverse need not.

The measurement model is `z = hᵀx + v` with `var(v) = R`. Correlated measurement noise must
be decorrelated before the call, since a scalar `R` cannot express it.

Throws `ArgumentError` if `x` or `h` disagrees in length with `D`, and `DomainError` for a
non-positive `R`.

Bierman's algorithm; see Bierman, *Factorization Methods for Discrete Sequential Estimation*
(1977), and Grewal & Andrews, *Kalman Filtering* (3e) §6.4.

# Example
```julia
using AstroSolve.UDU: udu_from_P, bierman_measurement_update!

U, D = udu_from_P([4.0 0.0; 0.0 1.0])
x = [0.0, 0.0]
K = bierman_measurement_update!(U, D, x, [1.0, 0.0], 0.25, 2.0)
```
"""
function bierman_measurement_update!(U::AbstractMatrix, D::AbstractVector,
                                      x::AbstractVector,
                                      h::AbstractVector, R::Real, ν::Real)
    n = length(D)
    length(x) == n || throw(ArgumentError(
        "bierman_measurement_update!: x must have length $(n) to match D; " *
        "got $(length(x))"))
    length(h) == n || throw(ArgumentError(
        "bierman_measurement_update!: h must have length $(n) to match D; " *
        "got $(length(h))"))
    R > 0 || throw(DomainError(R,
        "bierman_measurement_update!: measurement variance R must be positive; " *
        "got $(R)"))

    a = U' * h                # f vector,  length n
    b = D .* a                # v vector,  length n

    α = R + a[1] * b[1]
    γ = 1.0 / α
    D[1] = D[1] * R * γ

    K = zeros(Float64, n)
    K[1] = b[1]

    @inbounds for j in 2:n
        β = α
        α = β + a[j] * b[j]
        λ = -a[j] / β
        γ = 1.0 / α
        D[j] = D[j] * β * γ
        @inbounds for i in 1:(j - 1)
            τ        = U[i, j]
            U[i, j]  = τ + λ * K[i]
            K[i]     = K[i] + b[j] * τ
        end
        K[j] = b[j]
    end

    gain = K ./ α
    @inbounds for i in 1:n
        x[i] += gain[i] * ν
    end
    return gain
end

end # module UDU
