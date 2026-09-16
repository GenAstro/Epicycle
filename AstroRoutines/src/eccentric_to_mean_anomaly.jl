# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    eccentric_to_mean_anomaly(E::Real, e::Real)

Convert eccentric anomaly to mean anomaly for an elliptic orbit.

Evaluates Kepler's equation `M = E - e·sin(E)` (closed form).

# Arguments
- `E`: Eccentric anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.

# Returns
The mean anomaly `M` [rad], as the smooth continuation of `E` (not wrapped to `[0, 2π)`).

# Notes
- Differentiable with respect to `E` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ∉ [0, 1)`.
- This is the exact analytic inverse of `mean_to_eccentric_anomaly` and serves as its truth source.

# Example
```jldoctest
eccentric_to_mean_anomaly(deg2rad(125.0), 0.2)

# output

2.017831156135114
```
"""
function eccentric_to_mean_anomaly(E::Real, e::Real)
    _check_elliptic_eccentricity(e)
    return E - e * sin(E)
end
