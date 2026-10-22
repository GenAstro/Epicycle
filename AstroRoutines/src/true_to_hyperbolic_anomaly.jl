# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    true_to_hyperbolic_anomaly(ν::Real, e::Real)

Convert true anomaly to hyperbolic anomaly for a hyperbolic orbit.

Evaluates `H = asinh(√(e²−1)·sin ν / (1 + e·cos ν))`, the form of `tanh(H/2) = √((e−1)/(e+1))·tan(ν/2)`
that stays finite across ν = ±π.

# Arguments
- `ν`: True anomaly [rad], between the asymptotes, `1 + e·cos ν > 0`. A value and the same value
  less 2π give the same `H`.
- `e`: Eccentricity, `e > 1`.

# Returns
The hyperbolic anomaly `H`, negative before periapsis and positive after. It is not an angle and
does not wrap.

# Notes
- Differentiable with respect to `ν` and `e` (ForwardDiff).
- Throws `ArgumentError` if `e ≤ 1`, or if `ν` lies outside the asymptotes.
- Exact analytic inverse of `hyperbolic_to_true_anomaly`.

# Example
```jldoctest
true_to_hyperbolic_anomaly(deg2rad(60.0), 1.5)

# output

0.5283553629664819
```
"""
function true_to_hyperbolic_anomaly(ν::Real, e::Real)
    _check_hyperbolic_eccentricity(e)
    denominator = 1 + e * cos(ν)
    denominator > 0 || throw(ArgumentError(
        "true anomaly ν = $(ν) rad lies outside the asymptotes of a hyperbola with eccentricity " *
        "e = $(e); it must satisfy 1 + e·cos ν > 0."))
    return asinh(sqrt(e^2 - 1) * sin(ν) / denominator)
end
