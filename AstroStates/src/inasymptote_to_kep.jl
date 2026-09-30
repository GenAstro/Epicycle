# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    inasymptote_to_kep(inasym::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)

Convert incoming asymptote elements to Keplerian elements.

# Arguments
- `inasym`: incoming asymptote elements:
    - `rₚ`  : periapsis radius [length], > 0
    - `C₃`  : characteristic energy [length²/time²]
    - `λₐ` : right ascension of the asymptote [rad]
    - `δₐ` : declination of the asymptote [rad]
    - `θᵦ` : B-plane angle [rad]
    - `ν`  : true anomaly [rad]

- `μ`: Gravitational parameter [length³/time²]
- `tol`: Singularity tolerance (default = 1e-12)

# Returns
- Keplerian state vector `[a, e, i, Ω, ω, ν]`, angles in the ranges of [`cart_to_kep`](@ref)

# Notes
- Returns `fill(NaN, 6)`, with a warning, when the elements describe no orbit: C₃ ≈ 0
  (parabolic), rₚ ≤ 0, an elliptic C₃ with rₚ beyond the semi-major axis, a circular orbit, or an
  asymptote along the z-axis.
- Angles in radians. Units consistent with `μ`.

# Examples
```julia
inasym = [6778.0, 5.0, 0.0, π/4, π/2, π/2]
kep = inasymptote_to_kep(inasym, 398600.4418)
```
"""
inasymptote_to_kep(inasym::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12) =
    _asymptote_to_kep(inasym, μ, -1, tol)
