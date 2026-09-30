# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    kep_to_cart(state::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)

Convert a Keplerian state vector to a Cartesian state vector.

# Arguments
- `state`: Keplerian elements `[a, e, i, Ω, ω, ν]`
  * `a`: semi-major axis; negative for a hyperbola
  * `e`: eccentricity
  * `i`: inclination
  * `Ω`: right ascension of ascending node
  * `ω`: argument of periapsis
  * `ν`: true anomaly
- `μ`: Gravitational parameter
- `tol`: Tolerance for singularities like p ≈ 0 (default: 1e-12)

# Returns
A 6-element vector `[x, y, z, vx, vy, vz]` representing Cartesian position and velocity. The
element type follows the inputs.

# Examples
```julia
kep = [7000.0, 0.01, π/4, 0.0, 0.0, π/3]
cart = kep_to_cart(kep, 398600.4418)
```

# Notes
- Angles must be in radians.
- Dimensional quantities must be consistent units with μ.
- Returns a vector of `NaN`s, with a warning, when the conversion is undefined: a parabolic or
  collapsed orbit, or a hyperbolic true anomaly at or beyond the asymptote
  (`1 + e cos ν ≤ 0`), which would give a negative radius.
"""
function kep_to_cart(state::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)
    if length(state) != 6
        error("Input vector must have exactly six elements: a, e, i, Ω, ω, ν.")
    end
    T = float(promote_type(eltype(state), typeof(μ)))

    if μ < tol
        @warn "Conversion Failed: μ < tolerance."
        return fill(T(NaN), 6)
    end

    # Unpack the elements
    a, e, i, Ω, ω, ν = state

    # Semi-latus rectum: p = a * (1 - e²)
    p = a * (1 - e^2)

    # Degenerate orbit (parabolic or collapsed)
    if p < tol || abs(1 - e) < tol
        @warn "Conversion Failed: Orbit is parabolic or singular."
        return fill(T(NaN), 6)
    end

    # Radial distance: r = p / (1 + e cos ν). On a hyperbola the denominator reaches zero at the
    # asymptote; beyond it r would be negative and the state a mirror image of no real point.
    denom = 1 + e * cos(ν)
    if denom <= tol
        @warn "Conversion Failed: True anomaly $(ν) is at or beyond the asymptote of a hyperbola " *
              "with eccentricity $(e)."
        return fill(T(NaN), 6)
    end
    r = p / denom

    # Position and velocity in the perifocal frame
    factor = sqrt(μ / p)
    sν, cν = sincos(ν)
    r̄ₚ = SVector{3,T}(r * cν, r * sν, 0)
    v̄ₚ = SVector{3,T}(-factor * sν, factor * (e + cν), 0)

    # Rotation from perifocal to inertial
    sΩ, cΩ = sincos(Ω)
    sω, cω = sincos(ω)
    si, ci = sincos(i)
    R = SMatrix{3,3,T}(cω * cΩ - sω * ci * sΩ,  cω * sΩ + sω * ci * cΩ,  sω * si,
                       -sω * cΩ - cω * ci * sΩ, -sω * sΩ + cω * ci * cΩ, cω * si,
                       si * sΩ,                 -si * cΩ,                ci)      # column-major

    pos = R * r̄ₚ
    vel = R * v̄ₚ
    return T[pos[1], pos[2], pos[3], vel[1], vel[2], vel[3]]
end
