# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    kep_to_modkep(kep::AbstractVector{<:Real}; tol::Real = 1e-12)

Convert a classical Keplerian state to a Modified Keplerian state.

# Arguments
- `kep::AbstractVector{<:Real}`: Keplerian state vector `[a, e, i, Ω, ω, ν]`
    - `a`  : semi-major axis 
    - `e`  : eccentricity
    - `i`  : inclination 
    - `Ω`  : right ascension of ascending node 
    - `ω`  : argument of periapsis 
    - `ν`  : true anomaly 

- `tol::Real`: tolerance for singularity and consistency checks (default = 1e-12)

# Returns
Modified Keplerian state `[rₚ, rₐ, i, Ω, ω, ν]` or `fill(NaN, 6)` if invalid.

# Notes
- Parabolic orbits (`e ≈ 1`) and singular conics are not supported.
- Units must be consistent. Angles in radians.

# Examples
```julia
kep = [7000.0, 0.01, π/4, 0.0, 0.0, π/3]
modkep = kep_to_modkep(kep)
```
"""
function kep_to_modkep(kep::AbstractVector{<:Real}; tol::Real = 1e-12)
    if length(kep) != 6
        error("Input vector must contain six elements: [a, e, i, Ω, ω, ν]; got $(length(kep)) elements.")
    end
    _has_nan(kep) && return fill(float(eltype(kep))(NaN), 6)

    a, e, i, Ω, ω, ν = kep

    if _inclination_out_of_range(i)
        @warn "Conversion failed: Inclination ($i) is outside [0, π]."
        return fill(float(eltype(kep))(NaN), 6)
    end

    # Check for parabolic or undefined orbit
    if abs(1.0 - e) < tol
        @warn "Conversion failed: Orbit is parabolic."
        return fill(float(eltype(kep))(NaN), 6)
    end
    if e < 0
        @warn "Conversion failed: Eccentricity ($e) cannot be negative."
        return fill(float(eltype(kep))(NaN), 6)
    end
    if abs(a) < tol
        @warn "Conversion failed: Semimajor axis ($a) is zero, a collapsed orbit."
        return fill(float(eltype(kep))(NaN), 6)
    end

    # Ensure a and e are compatible
    if a > 0.0 && e > 1.0
        @warn "Conversion failed: Semimajor axis ($a) cannot be positive if eccentricity ($e) > 1."
        return fill(float(eltype(kep))(NaN), 6)
    elseif a < 0.0 && e < 1.0
        @warn "Conversion failed: Semimajor axis ($a) cannot be negative if eccentricity ($e) < 1."
        return fill(float(eltype(kep))(NaN), 6)
    end

    # On a hyperbola the true anomaly stops short of the asymptote, where 1 + e cos ν = 0, as in
    # kep_to_cart
    if e > 1 && 1 + e * cos(ν) <= tol
        @warn "Conversion failed: True anomaly $(ν) is at or beyond the asymptote of a hyperbola " *
              "with eccentricity $(e)."
        return fill(float(eltype(kep))(NaN), 6)
    end

    rₚ = a * (1 - e)
    rₐ = a * (1 + e)

    return float(eltype(kep))[rₚ, rₐ, i, Ω, ω, ν]
end
