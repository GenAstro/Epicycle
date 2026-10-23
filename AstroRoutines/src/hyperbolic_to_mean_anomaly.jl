# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    hyperbolic_to_mean_anomaly(H::Real, e::Real)

Convert hyperbolic anomaly to mean anomaly for a hyperbolic orbit.

Evaluates the hyperbolic Kepler equation `M = e·sinh H − H` in closed form.

# Arguments
- `H`: Hyperbolic anomaly.
- `e`: Eccentricity, `e > 1`.

# Returns
The mean anomaly `M`, negative before periapsis and positive after. It is not an angle and does
not wrap.

# Notes
- Differentiable with respect to `H` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ≤ 1`.
- Exact; the truth the iterative `mean_to_hyperbolic_anomaly` is checked against.

# Example
```jldoctest
hyperbolic_to_mean_anomaly(1.0, 1.5)

# output

0.7628017904657021
```
"""
function hyperbolic_to_mean_anomaly(H::Real, e::Real)
    _check_hyperbolic_eccentricity(e)
    return e * sinh(H) - H
end
