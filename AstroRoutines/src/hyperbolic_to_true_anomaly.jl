# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    hyperbolic_to_true_anomaly(H::Real, e::Real)

Convert hyperbolic anomaly to true anomaly for a hyperbolic orbit.

Evaluates the half-angle relation `ν = 2·atan(√((e+1)/(e−1))·tanh(H/2))`.

# Arguments
- `H`: Hyperbolic anomaly.
- `e`: Eccentricity, `e > 1`.

# Returns
The true anomaly `ν` [rad], in `(−ν∞, ν∞)` with `cos ν∞ = −1/e`; negative before periapsis.

# Notes
- Differentiable with respect to `H` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ≤ 1`.
- Exact analytic inverse of `true_to_hyperbolic_anomaly`.

# Example
```jldoctest
hyperbolic_to_true_anomaly(0.5283553629664819, 1.5)

# output

1.0471975511965976
```
"""
function hyperbolic_to_true_anomaly(H::Real, e::Real)
    _check_hyperbolic_eccentricity(e)
    return 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(H / 2))
end
