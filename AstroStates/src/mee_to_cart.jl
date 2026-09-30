# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    mee_to_cart(mod_equinoct::AbstractVector{<:Real}, μ::Real; j::Real = 1.0, tol::Real = 1e-12)

Convert Modified Equinoctial Elements to Cartesian state.

# Arguments
- `mod_equinoct`: the Modified Equinoctial Elements `[p, f, g, h, k, L]`
- `μ::Real`: Gravitational parameter
- `j::Real=1.0`: retrograde factor, 1 for the prograde set and -1 for the retrograde set
- `tol::Real`: tolerance for singularity checking

# Returns
- A 6-element vector `[x, y, z, vx, vy, vz]` representing the Cartesian position and velocity.

The elements describe no state, and the result is `NaN`s with a warning, when μ or `p` is not
positive, or when `1 + f cos L + g sin L ≤ 0`: a hyperbolic true longitude at or beyond the
asymptote, which would give a negative radius.

# Examples
```julia
mee = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
cart = mee_to_cart(mee, 398600.4418)
```
"""
function mee_to_cart(mod_equinoct::AbstractVector{<:Real}, μ::Real; j::Real = 1.0, tol::Real = 1e-12)
    if length(mod_equinoct) != 6
        error("Input vector must have exactly six elements: [p, f, g, h, k, L].")
    end
    if !(j == 1 || j == -1)
        error("Invalid value for j: must be 1.0 or -1.0")
    end
    T = float(promote_type(eltype(mod_equinoct), typeof(μ)))

    p, f, g, h, k, L = mod_equinoct

    if μ < tol
        @warn "Conversion failed: μ < tolerance."
        return fill(T(NaN), 6)
    end
    if p < tol
        @warn "Conversion failed: Semi-latus rectum p = $(p) must be positive."
        return fill(T(NaN), 6)
    end

    # Radius
    sL, cL = sincos(L)
    w = 1 + f * cL + g * sL
    if w <= tol
        @warn "Conversion failed: True longitude $(L) is at or beyond the asymptote of the orbit."
        return fill(T(NaN), 6)
    end
    r = p / w

    # Position and velocity in the orbital plane
    X1 = r * cL
    Y1 = r * sL
    dotX1 = -sqrt(μ / p) * (g + sL)
    dotY1 =  sqrt(μ / p) * (f + cL)

    # Equinoctial frame
    α2 = h^2 - k^2
    s2 = 1 + h^2 + k^2
    f̂ = SVector{3,T}(1 + α2, 2k * h, -2k * j) / s2
    ĝ = SVector{3,T}(2k * h * j, (1 - α2) * j, 2h) / s2

    # Cartesian position and velocity
    reci = X1 * f̂ + Y1 * ĝ
    veci = dotX1 * f̂ + dotY1 * ĝ
    return T[reci[1], reci[2], reci[3], veci[1], veci[2], veci[3]]
end
