# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    outasymptote_to_kep(outasym::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12)

Convert outgoing asymptote elements to Keplerian elements.

# Arguments
- `outasym`: outgoing asymptote elements:
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
outasym = [6778.0, 5.0, 0.0, π/4, π/2, π/2]
kep = outasymptote_to_kep(outasym, 398600.4418)
```
"""
outasymptote_to_kep(outasym::AbstractVector{<:Real}, μ::Real; tol::Real=1e-12) =
    _asymptote_to_kep(outasym, μ, 1, tol)

# The outgoing (dir = 1) and incoming (dir = -1) asymptotes differ only in which side of the
# asymptote the eccentricity vector lies.
function _asymptote_to_kep(asym::AbstractVector{<:Real}, μ::Real, dir::Int, tol::Real)
    if length(asym) != 6
        error("Input must be a 6-element vector: [rₚ, C₃, λₐ, δₐ, θᵦ, ν]")
    end
    T = float(promote_type(eltype(asym), typeof(μ)))

    rₚ, c₃, λₐ, δₐ, θᵦ, ν = asym

    # Parabolic orbits cannot be represented by asymptote parameters
    if abs(c₃) < tol
        @warn "Conversion failed: Orbit is nearly parabolic."
        return fill(T(NaN), 6)
    end
    if rₚ <= 0
        @warn "Conversion failed: Periapsis radius $(rₚ) must be positive."
        return fill(T(NaN), 6)
    end

    # Semi-major axis from energy, eccentricity from periapsis radius
    a = -μ / c₃
    e = 1 - rₚ / a
    if c₃ < 0 && rₚ > a
        @warn "Conversion failed: Periapsis radius $(rₚ) exceeds the semi-major axis $(a) of the " *
              "elliptic orbit that C₃ = $(c₃) gives."
        return fill(T(NaN), 6)
    end
    if e < tol
        @warn "Conversion failed: Orbit is nearly circular."
        return fill(T(NaN), 6)
    end

    # Asymptote direction unit vector
    sδ, cδ = sincos(δₐ)
    sλ, cλ = sincos(λₐ)
    ŝ = SVector{3,T}(cδ * cλ, cδ * sλ, sδ)

    # The B-plane axes need an asymptote off the z-axis
    sxy = hypot(ŝ[1], ŝ[2])
    if sxy < tol
        @warn "Conversion failed: Asymptote vector is aligned with z-axis."
        return fill(T(NaN), 6)
    end

    # B-plane axes, and the angular momentum direction from the B-plane angle
    Ê = SVector{3,T}(-ŝ[2], ŝ[1], 0) / sxy          # ẑ × ŝ, normalised
    N̂ = cross(ŝ, Ê)
    sθ, cθ = sincos(θᵦ)
    ĥ = cθ * Ê + sθ * N̂                             # sin(π/2 - θᵦ) Ê + cos(π/2 - θᵦ) N̂

    # Eccentricity direction
    if c₃ <= -tol
        # Elliptic: the "asymptote" is the apoapsis direction
        ê = -ŝ
    else
        # Hyperbolic: periapsis is the asymptote turned back through the true anomaly of the
        # asymptote, νₘ = acos(-1/e), written with atan to keep its derivative.
        νₘ = atan(sqrt(e^2 - 1), -one(T))
        ô = cross(ĥ, ŝ)
        ê = -dir * sin(νₘ) * ô + cos(νₘ) * ŝ
    end

    i, Ω, ω, _, _ = _orbit_orientation(ĥ, e * ê, e, tol)
    return T[a, e, i, Ω, ω, mod(ν, 2 * T(π))]
end
