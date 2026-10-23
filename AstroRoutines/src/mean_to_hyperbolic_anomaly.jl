# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    mean_to_hyperbolic_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)

Convert mean anomaly to hyperbolic anomaly for a hyperbolic orbit by solving the hyperbolic Kepler
equation `M = e·sinh H − H` with Newton iteration.

# Arguments
- `M`: Mean anomaly.
- `e`: Eccentricity, `e > 1`.

# Keywords
- `tol`: Convergence tolerance on the Newton step, relative to `max(1, |H|)`.
- `maxiter`: Maximum number of Newton iterations.

# Returns
The hyperbolic anomaly `H`, with the sign of `M`.

# Notes
- Starts from `H₀ = sign(M)·ln(2|M|/e + 1.8)`, after Danby (1987), which stays near the root
  even close to a parabola, where the derivative `e·cosh H − 1` vanishes at `H = 0`. Newton
  converges within 20 iterations for `e` from `1 + 10⁻⁶` to `10⁴` and `|M|` up to `10⁵`.
- Differentiable with respect to `M` and `e` (ForwardDiff); convergence is judged on the value, so
  the derivatives are those of the converged root.
- Throws `ArgumentError` if `e ≤ 1`, and an error if the iteration does not converge.

# Example
```jldoctest
mean_to_hyperbolic_anomaly(1.0, 1.5)

# output

1.1616354445046073
```
"""
function mean_to_hyperbolic_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)
    _check_hyperbolic_eccentricity(e)

    M, e = promote(M, e)
    H = sign(M) * log(2 * abs(M) / e + 1.8)

    for _ in 1:maxiter
        ΔH = (e * sinh(H) - H - M) / (e * cosh(H) - 1)
        H -= ΔH
        if abs(ΔH) < tol * max(1, abs(H))
            return H
        end
    end

    error("mean_to_hyperbolic_anomaly: Newton iteration did not converge in $(maxiter) iterations " *
          "(M = $(M), e = $(e), tol = $(tol)).")
end
