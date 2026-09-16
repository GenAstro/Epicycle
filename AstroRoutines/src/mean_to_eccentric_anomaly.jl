# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    mean_to_eccentric_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)

Convert mean anomaly to eccentric anomaly for an elliptic orbit.

Solves Kepler's equation `M = E - e·sin(E)` for `E`.

# Arguments
- `M`: Mean anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.
- `tol`: Convergence tolerance in radians on the Newton step `|ΔE|`; default `1e-12`.
- `maxiter`: Maximum Newton iterations (default `50`).

# Returns
The eccentric anomaly `E` [rad], as the smooth continuation of `M` (not wrapped to `[0, 2π)`).

# Notes
- Differentiable with respect to `M` and `e` (ForwardDiff): the analytic partials are
  `dE/dM = 1 / (1 - e·cos(E))` and `dE/de = sin(E) / (1 - e·cos(E))`.
- Throws `ArgumentError` if `e ∉ [0, 1)`.
- Throws `ErrorException` if the iteration does not converge within `maxiter`.

# Example
```jldoctest
mean_to_eccentric_anomaly(deg2rad(120.0), 0.2)

# output

2.250008654463779
```
"""
function mean_to_eccentric_anomaly(M::Real, e::Real; tol::Real=1e-12, maxiter::Integer=50)
    _check_elliptic_eccentricity(e)

    M, e = promote(M, e)
    # Danby (1987) starter: E₀ = M + sign(sin M)·0.85·e. Robust across 0 ≤ e < 1 — the naive
    # E₀ = M makes Newton overshoot at high e near periapsis (where 1 − e·cos E → 0) and diverge.
    # The starter only affects the iteration path, not the converged value or its AD partials.
    E = M + sign(sin(M)) * (0.85 * e)

    for _ in 1:maxiter
        ΔE = (E - e * sin(E) - M) / (1 - e * cos(E))
        E -= ΔE
        if abs(ΔE) < tol
            return E
        end
    end

    error("mean_to_eccentric_anomaly: Newton iteration did not converge in $(maxiter) iterations " *
          "(M = $(M), e = $(e), tol = $(tol)).")
end
