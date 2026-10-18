# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    kep_to_delaunay(kep::AbstractVector{<:Real}, μ::Real)

Convert a classical Keplerian state to a Delaunay state.

# Arguments
- `kep::AbstractVector{<:Real}`: Keplerian state vector `[a, e, i, Ω, ω, ν]`
- `μ::Real`: gravitational parameter

# Returns
Delaunay state `[l, g, h, L, G, H]`, angles in radians and momenta in km²/s, or `fill(NaN, 6)`
if invalid.

# Notes
- Defined for closed orbits and positive μ; otherwise logs a warning and returns `NaN`s.

# Examples
```julia
kep = [7000.0, 0.01, 0.5, 1.0, 2.0, 0.3]
d = kep_to_delaunay(kep, 398600.4418)
```
"""
function kep_to_delaunay(kep::AbstractVector{<:Real}, μ::Real)
    if length(kep) != 6
        error("Input vector must contain six elements: [a, e, i, Ω, ω, ν]")
    end
    T = float(promote_type(eltype(kep), typeof(μ)))

    a, e, i, Ω, ω, ν = kep
    if !(0 <= e < 1) || a <= 0 || μ <= 0
        @warn "Conversion failed: Delaunay elements are defined for a closed orbit and positive μ; " *
              "got e = $e, a = $a, μ = $μ."
        return fill(T(NaN), 6)
    end

    L = sqrt(μ * a)
    G = L * sqrt(1 - e^2)
    return T[_wrap_2pi(_true_to_mean_anom(ν, e)), ω, Ω, L, G, G * cos(i)]
end
