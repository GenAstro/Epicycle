# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

using LinearAlgebra

"""
    cart_to_outasymptote(cart::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12) -> Vector{Float64}

Convert a Cartesian state to outgoing asymptote parameters.

# Arguments
- `cart::Vector{<:Real}`: Cartesian state `[x, y, z, vx, vy, vz]`
- `μ::Real`: gravitational parameter
- `tol::Float64`: tolerance for detecting singularities

# Returns
A 6-element vector: `[rₚ, C₃, λₐ, δₐ, θᵦ, ν]`  
- rₚ : radius of periapsis
- C₃ : characteristic energy
- λₐ : right ascension of outgoing asymptote
- δₐ : declination of outgoing asymptote
- θᵦ : B-plane angle
- ν : true anomaly

# Notes
- Angles are in radians.
- Dimensional quantities are consistent with units of μ 
"""
function cart_to_outasymptote(cart::Vector{<:Real}, μ::Real; tol::Float64 = 1e-12)
    if length(cart) != 6
        error("Input must be a 6-element Cartesian vector: [x, y, z, vx, vy, vz]")
    end

    # Unpack position and velocity 
    r̄ = cart[1:3]
    v̄ = cart[4:6]

    # Compute position, velocity, angular momentum magnitudes
    r = norm(r̄)
    v = norm(v̄)
    h̄ = cross(r̄, v̄)
    h = norm(h̄)

    # Check for degenerate cases
    if r < tol || v < tol 
        @warn "Conversion failed: Orbit is singular due to degenerate position or velocity vector."
        return fill(NaN, 6)
    end
    if h < tol 
        @warn "Conversion failed: Orbit is singular with zero angular momentum."
        return fill(NaN, 6)
    end

    # Compute eccentricity vector and magnitude
    ē = (cross(v̄, h̄) ./ μ) .- (r̄ ./ r)  # TODO
    e = norm(ē)
    if e < tol 
        @warn "Conversion failed: Orbit is circular."
        return fill(NaN, 6)
    end

    # Compute characteristic energy and check for parabolic orbit
    C₃ = v^2 - 2μ / r
    if isapprox(C₃, 0.0; atol=tol)
        @warn "Conversion failed: Orbit is parabolic."
        return fill(NaN, 6)
    end

    # Compute radius of periapsis and check for singular orbit
    a = -μ / C₃
    rₚ = a * (1 - e)
    if rₚ < tol
        @warn "Conversion failed: Orbit is singular due to near-zero periapsis radius."
        return fill(NaN, 6)
    end

    # Build incoming asymptote unit vector ŝ
    fac = 1 / (1 + (C₃ * h^2) / μ^2)
    if C₃ > tol
        ŝ = fac * (sqrt(C₃)/μ * cross(h̄, ē) - ē )
    else
        ŝ = -ē / e
    end

    # Check if asymptote is aligned with z-vector
    ẑ = [0.0, 0.0, 1.0]
    angle_with_z = acos(clamp(abs(dot(ŝ,ẑ)), -1.0, 1.0))
    if angle_with_z < tol
        @warn "Conversion failed: Asymptote vector is aligned with the z-axis."
        return fill(NaN, 6)
    end

    # Construct B-plane coordinates and the B-plane angle
    Ê = normalize(cross(ẑ, ŝ))
    N̂ = cross(ŝ, Ê)
    b̄ = cross(h̄, ŝ)
    sinθᵦ = dot(b̄, Ê) / h
    cosθᵦ = dot(b̄, -N̂) / h
    θᵦ = atan(sinθᵦ, cosθᵦ)
    θᵦ = θᵦ < 0 ? θᵦ + 2π : θᵦ

    # Asymptote orientation angles
    δₐ = asin(ŝ[3])
    λₐ = atan(ŝ[2], ŝ[1])
    λₐ = λₐ < 0 ? λₐ + 2π : λₐ

    # True anomaly
    ν = acos(clamp(dot(ē, r̄) / (e * r), -1.0, 1.0))
    ν = dot(v̄, r̄) < 0 ? 2π - ν : ν
 
    return [rₚ, C₃, λₐ, δₐ, θᵦ, ν]
end
