# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier:

"""
    mee_to_cart(mod_equinoct::Vector{<:Real}, μ::Real; j::Real = 1.0)

Convert Modified Equinoctial Elements to Cartesian state.

# Arguments
- `mod_equinoct::Vector{<:Real}`: Vector containing the Modified Equinoctial Elements [p, f, g, h, k, L]
- `μ::Real`: Gravitational parameter
- `j::Real=1.0`: Optional constant (1 for prograde, -1 for retrograde), defaults to `1.0`

# Returns
- A 6-element vector `[x, y, z, vx, vy, vz]` representing the Cartesian position and velocity.

# Examples
```julia
mee = [7000.0, 0.01, 0.0, 0.1, 0.0, π/4]
cart = mee_to_cart(mee, 398600.4418)
```

"""
function mee_to_cart(mod_equinoct::Vector{<:Real}, μ::Real; j::Real = 1.0)
    if length(mod_equinoct) != 6
        error("Input vector must have exactly six elements: [p, f, g, h, k, L].")
    end

    # Extract Modified Equinoctial Elements
    p, f, g, h, k, L = mod_equinoct

    # Validate j
    if j ∉ (-1.0, 1.0)
        error("Invalid value for j: must be 1.0 or -1.0")
    end

    # Validate semi-latus rectum
    if p < 0
        error("Semi-latus rectum must be greater than 0")
    end

    # Radius computation
    r = p / (1 + f * cos(L) + g * sin(L))

    # Position in orbital plane
    X1 = r * cos(L)
    Y1 = r * sin(L)

    # Velocity in orbital plane
    if p == 0
        dotX1 = 0.0
        dotY1 = 0.0
    else
        dotX1 = -sqrt(μ / p) * (g + sin(L))
        dotY1 =  sqrt(μ / p) * (f + cos(L))
    end

    # Direction cosine transformation
    α2 = h^2 - k^2
    s2 = 1 + h^2 + k^2

    f̂ = [(1 + α2), 2k*h, -2k*j] ./ s2
    ĝ = [2k*h*j, (1 - α2)*j, 2h] ./ s2

    # Cartesian position and velocity
    reci = X1 * f̂ + Y1 * ĝ
    veci = dotX1 * f̂ + dotY1 * ĝ

    return vcat(reci, veci)
end
