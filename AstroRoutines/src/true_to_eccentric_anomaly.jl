# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    true_to_eccentric_anomaly(ν::Real, e::Real)

Convert true anomaly to eccentric anomaly for an elliptic orbit.

Evaluates the half-angle relation `E = 2·atan(√(1−e)·sin(ν/2), √(1+e)·cos(ν/2))`.

# Arguments
- `ν`: True anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.

# Returns
The eccentric anomaly `E` [rad] in the principal branch `(−π, π]`.

# Notes
- Differentiable with respect to `ν` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ∉ [0, 1)`.
- Exact analytic inverse of `eccentric_to_true_anomaly` on a single revolution.

# Example
```jldoctest
true_to_eccentric_anomaly(deg2rad(130.0), 0.2)

# output

2.1037839463996946
```
"""
function true_to_eccentric_anomaly(ν::Real, e::Real)
    _check_elliptic_eccentricity(e)
    return 2 * atan(sqrt(1 - e) * sin(ν / 2), sqrt(1 + e) * cos(ν / 2))
end
