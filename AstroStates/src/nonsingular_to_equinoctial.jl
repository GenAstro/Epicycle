# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    nonsingular_to_equinoctial(ns::AbstractVector{<:Real})

Convert a nonsingular Keplerian state to an equinoctial state.

# Arguments
- `ns::AbstractVector{<:Real}`: nonsingular Keplerian state vector `[a, e1, e2, e3, e4, e5]`
    - `a`  : semi-major axis
    - `e1`, `e2` : e cos(Ω + ω), e sin(Ω + ω)
    - `e3`, `e4` : sin(i/2) sin Ω, sin(i/2) cos Ω
    - `e5` : mean longitude Ω + ω + M

# Returns
Equinoctial state `[a, h, k, p, q, λ]`, or `fill(NaN, 6)` if invalid.

# Notes
- The eccentricity pair carries over unchanged, (h, k) = (e2, e1), and the inclination pair
  scales by 1/cos(i/2) = 1/√(1 - e3² - e4²). Neither step divides by e or sin i, so the
  conversion and its derivative are defined at e = 0 and i = 0.
- Defined for closed orbits away from i = π, where cos(i/2) is zero; `a ≤ 0`, `e ≥ 1` or
  `sin(i/2) ≥ 1` logs a warning and returns `NaN`s.

# Examples
```julia
ns = [7000.0, 0.0, 0.01, 0.0995, 0.0, π/4]
eq = nonsingular_to_equinoctial(ns)
```
"""
function nonsingular_to_equinoctial(ns::AbstractVector{<:Real})
    if length(ns) != 6
        error("Input vector must contain six nonsingular Keplerian elements: [a, e1, e2, e3, e4, e5]")
    end
    T = float(eltype(ns))

    # A failed conversion upstream has already warned; pass its NaNs on without a second warning
    any(isnan, ns) && return fill(T(NaN), 6)

    a, e1, e2, e3, e4, e5 = ns
    if a <= 0 || !(e1^2 + e2^2 < 1) || !(e3^2 + e4^2 < 1)
        @warn "Conversion failed: Nonsingular Keplerian elements describe no closed orbit away " *
              "from i = π; got a = $a, e = $(sqrt(e1^2 + e2^2)), sin(i/2) = $(sqrt(e3^2 + e4^2))."
        return fill(T(NaN), 6)
    end

    sec_half_i = 1 / sqrt(1 - e3^2 - e4^2)
    return T[a, e2, e1, e3 * sec_half_i, e4 * sec_half_i, _wrap_2pi(mod(e5, 2 * T(π)))]
end
