# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using LinearAlgebra

"""
    cart_to_sphazfpa(cart::AbstractVector{<:Real}; tol::Real = 1e-12)

Convert a Cartesian state to Spherical AZ-FPA representation.

# Arguments
- `cart`: Cartesian state `[x, y, z, vx, vy, vz]`
- `tol::Real`: Numerical tolerance for singularity checks (default: `1e-12`)

# Returns
A 6-element Spherical AZ-FPA state `[r, λ, δ, v, αₚ, ψ]`:
- `r`   : radial distance [length]
- `λ`   : right ascension [rad], in [0, 2π)
- `δ`   : declination [rad], in [-π/2, π/2]
- `v`   : velocity magnitude [length/time]
- `αₚ`  : flight path azimuth [rad], east of north in the local horizontal plane, in [0, 2π)
- `ψ`   : flight path angle [rad], measured from the radial direction, in [0, π]; π/2 is
          horizontal flight, less than π/2 is climbing. This is GMAT's convention.

# Notes
- Returns `fill(NaN, 6)`, with a warning, if `r` or `v` is near zero.
- All angles are in radians.
- For purely radial motion the azimuth is undefined; the value returned comes from rounding.

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
sphazfpa = cart_to_sphazfpa(cart)
```
"""
function cart_to_sphazfpa(cart::AbstractVector{<:Real}; tol::Real = 1e-12)
    if length(cart) != 6
        error("Input vector must have six elements: [x, y, z, vx, vy, vz]")
    end
    T = float(eltype(cart))

    r̄ = SVector{3,T}(cart[1], cart[2], cart[3])
    v̄ = SVector{3,T}(cart[4], cart[5], cart[6])

    r = norm(r̄)
    if r < tol
        @warn "Conversion failed: Position magnitude r = $r is below tolerance."
        return fill(T(NaN), 6)
    end

    v = norm(v̄)
    if v < tol
        @warn "Conversion failed: Velocity magnitude v = $v is below tolerance."
        return fill(T(NaN), 6)
    end

    rxy = hypot(r̄[1], r̄[2])
    λ = _wrap_2pi(atan(r̄[2], r̄[1]))
    δ = atan(r̄[3], rxy)

    # Flight path angle from the radial direction. atan keeps its precision near 0 and π, where
    # acos loses it; exactly radial motion is 0 or π with no derivative to give.
    c̄ = cross(r̄, v̄)
    ψ = iszero(c̄) ? (dot(r̄, v̄) >= 0 ? zero(T) : T(π)) : atan(norm(c̄), dot(r̄, v̄))

    # Azimuth from north toward east, in the local horizontal plane
    sδ, cδ = sincos(δ)
    sλ, cλ = sincos(λ)
    east  = SVector{3,T}(-sλ, cλ, 0)
    north = SVector{3,T}(-sδ * cλ, -sδ * sλ, cδ)
    αₚ = _wrap_2pi(atan(dot(v̄, east), dot(v̄, north)))

    return T[r, λ, δ, v, αₚ, ψ]
end
