# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    mean_to_true_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)

Convert mean anomaly to true anomaly for an elliptic orbit.

Solves Kepler's equation for `E`, then converts `E` to `ν`.

# Arguments
- `M`: Mean anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.
- `tol`: Convergence tolerance passed to the Kepler solve (default `1e-12`).
- `maxiter`: Maximum Newton iterations for the Kepler solve (default `50`).

# Returns
The true anomaly `ν` [rad] in the principal branch `(−π, π]`.

# Notes
- Differentiable with respect to `M` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ∉ [0, 1)`.
- Throws `ErrorException` if the Kepler solve does not converge within `maxiter`.

# Example
```jldoctest
mean_to_true_anomaly(deg2rad(120.0), 0.2)

# output

2.3975613842919206
```
"""
function mean_to_true_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)
    E = mean_to_eccentric_anomaly(M, e; tol=tol, maxiter=maxiter)
    return eccentric_to_true_anomaly(E, e)
end
