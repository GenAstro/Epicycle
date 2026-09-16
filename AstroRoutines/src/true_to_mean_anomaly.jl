# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    true_to_mean_anomaly(ν::Real, e::Real)

Convert true anomaly to mean anomaly for an elliptic orbit.

Converts `ν` to `E`, then evaluates Kepler's equation `M = E - e·sin(E)`.

# Arguments
- `ν`: True anomaly [rad].
- `e`: Eccentricity, `0 ≤ e < 1`.

# Returns
The mean anomaly `M` [rad], as the smooth continuation of `ν` (not wrapped to `[0, 2π)`).

# Notes
- Differentiable with respect to `ν` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ∉ [0, 1)`.
- Closed form, so there is no tolerance to set and the cost does not depend on `e`.

# Example
```jldoctest
true_to_mean_anomaly(deg2rad(130.0), 0.2)

# output

1.9315253702414643
```
"""
function true_to_mean_anomaly(ν::Real, e::Real)
    E = true_to_eccentric_anomaly(ν, e)
    return eccentric_to_mean_anomaly(E, e)
end
