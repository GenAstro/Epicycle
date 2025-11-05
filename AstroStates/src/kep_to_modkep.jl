# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    kep_to_modkep(kep::Vector{<:Real}; tol::Float64 = 1e-12)

Convert a classical Keplerian state to a Modified Keplerian state.

# Arguments
- `kep::Vector{<:Real}`: Keplerian state vector `[a, e, i, Ω, ω, ν]`
    - `a`  : semi-major axis 
    - `e`  : eccentricity
    - `i`  : inclination 
    - `Ω`  : right ascension of ascending node 
    - `ω`  : argument of periapsis 
    - `ν`  : true anomaly 

- `tol::Float64`: tolerance for singularity and consistency checks (default = 1e-12)

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
function kep_to_modkep(kep::Vector{<:Real}; tol::Float64 = 1e-12)
    if length(kep) != 6
        error("Input vector must contain six elements: [a, e, i, Ω, ω, ν]")
    end

    a, e, i, Ω, ω, ν = kep

    # Check for parabolic or undefined orbit
    if abs(1.0 - e) < tol
        @warn "Conversion failed: Orbit is parabolic."
        return fill(NaN, 6)
    end
    if e < 0
        @warn "Conversion failed: Eccentricity ($e) cannot be negative."
        return fill(NaN, 6)
    end

    # Ensure a and e are compatible
    if a > 0.0 && e > 1.0
        @warn "Conversion failed: Semimajor axis ($a) cannot be positive if eccentricity ($e) > 1."
        return fill(NaN, 6)
    elseif a < 0.0 && e < 1.0
        @warn "Conversion failed: Semimajor axis ($a) cannot be negative if eccentricity ($e) < 1."
        return fill(NaN, 6)
    end

    rₚ = a * (1 - e)
    rₐ = a * (1 + e)

    # Reject nearly singular conics or tiny periapsis
    if abs(rₚ) < tol
        @warn "Conversion failed: Singular conic section."
        return fill(NaN, 6)
    end

    return [rₚ, rₐ, i, Ω, ω, ν]
end
