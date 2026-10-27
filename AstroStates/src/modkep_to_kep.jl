# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    modkep_to_kep(modkep::AbstractVector{<:Real}; tol::Real = 1e-12)

Convert a Modified Keplerian state to a classical Keplerian state.

# Arguments
- `modkep::AbstractVector{<:Real}`: Modified Keplerian state vector `[rₚ, rₐ, i, Ω, ω, ν]`
    - `rₚ` : radius of periapsis
    - `rₐ` : radius of apoapsis 
    - `i`  : inclination 
    - `Ω`  : right ascension of ascending node 
    - `ω`  : argument of periapsis 
    - `ν`  : true anomaly

- `tol::Real`: tolerance for singularity and consistency checks (default = 1e-12)

# Returns
A classical Keplerian state `[a, e, i, Ω, ω, ν]` or `fill(NaN, 6)` if invalid.

# Notes
- Returns `NaN` vector if input describes a parabolic, undefined, or inconsistent orbit.
- Assumes input angles are in radians and distances use consistent units.

# Examples
```julia
modkep = [6778.0, 42164.0, π/6, 0.0, 0.0, 0.0]
kep = modkep_to_kep(modkep)
```
"""
function modkep_to_kep(modkep::AbstractVector{<:Real}; tol::Real = 1e-12)
    if length(modkep) != 6
        error("Input vector must contain six elements: [rₚ, rₐ, i, Ω, ω, ν]; got $(length(modkep)) elements.")
    end
    _has_nan(modkep) && return fill(float(eltype(modkep))(NaN), 6)

    rₚ, rₐ, i, Ω, ω, ν = modkep

    if _inclination_out_of_range(i)
        @warn "Conversion failed: Inclination ($i) is outside [0, π]."
        return fill(float(eltype(modkep))(NaN), 6)
    end

    # Check for singular orbits
    if abs(rₐ) < tol
        @warn "Conversion failed: Radius of apoapsis must not be zero."
        return fill(float(eltype(modkep))(NaN), 6)
    end
    # An infinite apoapsis is a parabola, which Keplerian elements cannot hold; a hyperbola has a
    # finite, negative rₐ
    if isinf(rₐ)
        @warn "Conversion failed: Radius of apoapsis is infinite, a parabolic orbit."
        return fill(float(eltype(modkep))(NaN), 6)
    end
    if rₐ < rₚ && rₐ > 0.0
        @warn "Conversion failed: Inconsistent Modified Keplerian state. If rₐ < rₚ, then rₐ must be negative."
        return fill(float(eltype(modkep))(NaN), 6)
    end
    if rₚ <= tol || abs(rₐ) <= tol
        @warn "Conversion failed: Singular conic section."
        return fill(float(eltype(modkep))(NaN), 6)
    end
    # A hyperbola has rₐ = a(1 + e) < 0 with |rₐ| = |a|(e + 1) > rₚ = |a|(e - 1). Equal radii give
    # e = ∞, and |rₐ| < rₚ gives e < 0.
    if rₐ < 0 && -rₐ <= rₚ + tol
        @warn "Conversion failed: Inconsistent Modified Keplerian state. For a hyperbola " *
              "(rₐ < 0), |rₐ| must exceed rₚ; got rₚ = $(rₚ), rₐ = $(rₐ)."
        return fill(float(eltype(modkep))(NaN), 6)
    end

    # Compute eccentricity and semi-major axis
    rp_by_ra = rₚ / rₐ
    e = (1 - rp_by_ra) / (1 + rp_by_ra)
    a = rₚ / (1 - e)
    # On a hyperbola the true anomaly stops short of the asymptote, where 1 + e cos ν = 0, as in
    # kep_to_cart
    if e > 1 && 1 + e * cos(ν) <= tol
        @warn "Conversion failed: True anomaly $(ν) is at or beyond the asymptote of a hyperbola " *
              "with eccentricity $(e)."
        return fill(float(eltype(modkep))(NaN), 6)
    end

    return [a, e, i, Ω, ω, ν]
end
