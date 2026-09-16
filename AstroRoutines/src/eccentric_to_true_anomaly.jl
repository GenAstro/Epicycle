# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    eccentric_to_true_anomaly(E::Real, e::Real)

Convert eccentric anomaly to true anomaly for an elliptic orbit.

Evaluates the half-angle relation `ν = 2·atan(√(1+e)·sin(E/2), √(1−e)·cos(E/2))`.

# Arguments
- `E`: Eccentric anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.

# Returns
The true anomaly `ν` [rad] in the principal branch `(−π, π]`.

# Notes
- Differentiable with respect to `E` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ∉ [0, 1)`.

# Example
```jldoctest
eccentric_to_true_anomaly(deg2rad(125.0), 0.2)

# output

2.337781537541991
```
"""
function eccentric_to_true_anomaly(E::Real, e::Real)
    _check_elliptic_eccentricity(e)
    return 2 * atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
end
