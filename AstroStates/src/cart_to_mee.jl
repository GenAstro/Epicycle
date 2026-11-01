# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

"""
    cart_to_mee(cart::AbstractVector{<:Real}, μ::Real; j::Real = 1.0, tol::Real = 1e-12)

Convert Cartesian state to Modified Equinoctial Elements (MEE).

# Arguments
- `cart`: 6-element vector `[x, y, z, vx, vy, vz]`
- `μ::Real`: Gravitational parameter
- `j::Real=1.0`: retrograde factor, 1 for the prograde set and -1 for the retrograde set, which
  is singular at i = 0 instead of i = π
- `tol::Real`: tolerance for singularity checking

# Returns
- A 6-element vector `[p, f, g, h, k, L]` representing the modified equinoctial elements, with
  the true longitude `L` in [0, 2π).

A singular state logs a warning and returns `NaN`s: μ below `tol`, a zero position, velocity or
angular momentum, or an orbit at the singularity of the chosen set (i = π for `j = 1`, i = 0 for
`j = -1`).

# Examples
```julia
cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
mee = cart_to_mee(cart, 398600.4418)
```
"""
function cart_to_mee(cart::AbstractVector{<:Real}, μ::Real; j::Real = 1.0, tol::Real = 1e-12)
    if length(cart) != 6
        error("Input vector must have exactly six elements: [x, y, z, vx, vy, vz]; got $(length(cart)) elements.")
    end
    if !(j == 1 || j == -1)
        error("Invalid value for j: must be 1.0 or -1.0; got j = $j.")
    end
    T = float(promote_type(eltype(cart), typeof(μ)))
    _has_nan(cart, μ) && return fill(T(NaN), 6)

    if μ < tol
        @warn "Conversion failed: μ < tolerance."
        return fill(T(NaN), 6)
    end

    r̄ = SVector{3,T}(cart[1], cart[2], cart[3])
    v̄ = SVector{3,T}(cart[4], cart[5], cart[6])
    r = norm(r̄)
    # Each check is !(x >= tol) so that a NaN magnitude counts as singular: under ForwardDiff, the
    # norm of an exactly zero vector is NaN, and x < tol would let it through.
    if !(r >= tol && norm(v̄) >= tol)
        @warn "Conversion failed: Orbit is singular due to degenerate position or velocity vector."
        return fill(T(NaN), 6)
    end

    # Angular momentum vector and magnitude
    h̄ = cross(r̄, v̄)
    h = norm(h̄)
    if !(h >= tol)
        @warn "Conversion failed: Orbit is singular due to degenerate angular momentum."
        return fill(T(NaN), 6)
    end
    ĥ = h̄ / h
    r̂ = r̄ / r

    # Eccentricity vector and semi-latus rectum
    ē = cross(v̄, h̄) / μ - r̂
    p = h^2 / μ

    # The set is singular where its inclination vector is infinite
    denom = 1 + ĥ[3] * j
    if abs(denom) < tol
        @warn "Conversion failed: modified equinoctial elements are singular at i = π for j = 1 and at i = 0 for j = -1."
        return fill(T(NaN), 6)
    end
    # Near the singularity the elements still convert, but h and k grow as 1 / (1 + j cos i) and
    # precision goes with them. The threshold is where the equinoctial set refuses outright,
    # within about 0.008° of the singular inclination.
    if abs(denom) < _EQUINOCTIAL_RETROGRADE_TOL
        singular, other = j == 1 ? ("π", -1) : ("0", 1)
        @warn "Modified equinoctial elements with j = $(Int(j)) are near their singularity at " *
              "i = $singular (1 + j cos i = $(denom)) and lose precision there; j = $other " *
              "suits this orbit."
    end

    # Equinoctial frame
    f̂ = SVector{3,T}(1 - ĥ[1]^2 / denom, -ĥ[1] * ĥ[2] / denom, -ĥ[1] * j)
    ĝ = cross(ĥ, f̂)

    # Modified equinoctial elements
    f  = dot(ē, f̂)
    g  = dot(ē, ĝ)
    hh = -ĥ[2] / denom
    k  =  ĥ[1] / denom

    # True longitude: the direction of the position in the equinoctial frame
    L = _plane_angle(r̂, f̂, ĝ)

    return T[p, f, g, hh, k, L]
end
