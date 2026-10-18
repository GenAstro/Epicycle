# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    equinoctial_to_nonsingular(eq::AbstractVector{<:Real})

Convert an equinoctial state to a nonsingular Keplerian state.

# Arguments
- `eq::AbstractVector{<:Real}`: equinoctial state vector `[a, h, k, p, q, λ]`
    - `a` : semi-major axis
    - `h`, `k` : e sin(Ω + ω), e cos(Ω + ω)
    - `p`, `q` : tan(i/2) sin Ω, tan(i/2) cos Ω
    - `λ` : mean longitude Ω + ω + M

# Returns
Nonsingular Keplerian state `[a, e1, e2, e3, e4, e5]`, or `fill(NaN, 6)` if invalid.

# Notes
- The eccentricity pair carries over unchanged, (e1, e2) = (k, h), and the inclination pair
  scales by cos(i/2) = 1/√(1 + p² + q²). Neither step divides by e or sin i, so the conversion
  and its derivative are defined at e = 0 and i = 0.
- Defined for closed orbits; `a ≤ 0` or `e ≥ 1` logs a warning and returns `NaN`s.

# Examples
```julia
eq = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
ns = equinoctial_to_nonsingular(eq)
```
"""
function equinoctial_to_nonsingular(eq::AbstractVector{<:Real})
    if length(eq) != 6
        error("Input vector must contain six equinoctial elements: [a, h, k, p, q, λ]")
    end
    T = float(eltype(eq))

    # A failed conversion upstream has already warned; pass its NaNs on without a second warning
    any(isnan, eq) && return fill(T(NaN), 6)

    a, h, k, p, q, λ = eq
    if a <= 0 || !(h^2 + k^2 < 1)
        @warn "Conversion failed: Nonsingular Keplerian elements are defined for a closed orbit; " *
              "got a = $a, e = $(sqrt(h^2 + k^2))."
        return fill(T(NaN), 6)
    end

    cos_half_i = 1 / sqrt(1 + p^2 + q^2)
    return T[a, k, h, p * cos_half_i, q * cos_half_i, _wrap_2pi(mod(λ, 2 * T(π)))]
end
