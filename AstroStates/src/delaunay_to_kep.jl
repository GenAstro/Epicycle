# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    delaunay_to_kep(d::AbstractVector{<:Real}, μ::Real)

Convert a Delaunay state to a classical Keplerian state.

# Arguments
- `d::AbstractVector{<:Real}`: Delaunay state vector `[l, g, h, L, G, H]`
- `μ::Real`: gravitational parameter

# Returns
Keplerian state `[a, e, i, Ω, ω, ν]`, true anomaly in [0, 2π), or `fill(NaN, 6)` if invalid.

# Notes
- Elements that describe no orbit (L ≤ 0, G outside (0, L], or |H| > G) or a non-positive μ
  log a warning and return `NaN`s.

# Examples
```julia
d = kep_to_delaunay([7000.0, 0.01, 0.5, 1.0, 2.0, 0.3], 398600.4418)
kep = delaunay_to_kep(d, 398600.4418)
```
"""
function delaunay_to_kep(d::AbstractVector{<:Real}, μ::Real)
    if length(d) != 6
        error("Input vector must contain six Delaunay elements: [l, g, h, L, G, H]")
    end
    T = float(promote_type(eltype(d), typeof(μ)))

    l, g, h, L, G, H = d
    if μ <= 0 || L <= 0 || !(0 < G <= L) || abs(H) > G
        @warn "Conversion failed: Delaunay elements describe no orbit; got L = $L, G = $G, " *
              "H = $H, μ = $μ."
        return fill(T(NaN), 6)
    end

    a = L^2 / μ
    e = sqrt(max(1 - (G / L)^2, zero(T)))
    i = acos(clamp(H / G, -one(T), one(T)))
    E = _eccentric_anomaly(l, e)
    ν = atan(sqrt(1 - e^2) * sin(E), cos(E) - e)
    return T[a, e, i, h, g, _wrap_2pi(ν)]
end
