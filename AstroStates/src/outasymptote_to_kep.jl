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

"""
    outasymptote_to_kep(outasym::Vector{<:Real}, μ::Real; tol::Float64=1e-12)

Convert outgoing asymptote elements to Keplerian elements.

# Arguments
- `outasym`: Vector{<:Real} of outgoing asymptote elements:
    - `rₚ`  : periapsis radius [length]
    - `C₃`  : characteristic energy [length²/time²]
    - `λₐ` : right ascension of the asymptote [rad]
    - `δₐ` : declination of the asymptote [rad]
    - `θᵦ` : B-plane angle [rad]
    - `ν`  : true anomaly [rad]

- `μ`: Gravitational parameter [length³/time²]
- `tol`: Singularity tolerance (default = 1e-12)

# Returns
- Keplerian state vector `[a, e, i, Ω, ω, ν]`

# Notes
- Returns `fill(NaN, 6)` if singularity is detected.
- Angles in radians. Units consistent with `μ`.
"""
function outasymptote_to_kep(outasym::Vector{<:Real}, μ::Real; tol::Float64=1e-12)
    if length(outasym) != 6
        error("Input must be a 6-element vector: [a, e, i, Ω, ω, ν]")
    end

    # Unpack outgoing asymptote elements
    rₚ, c₃, λₐ, δₐ, θᵦ, ν = outasym

    # Compute semi-major axis from energy
    a = -μ / c₃

    # Compute eccentricity from periapsis radius
    e = 1 - rₚ / a

    # Parabolic or circular orbits cannot be represented by asymptote parameters
    if abs(c₃) < tol
        @warn "Conversion failed: Orbit is nearly parabolic."
        return fill(NaN, 6)
    end
    if e < tol
        @warn "Conversion failed: Orbit is nearly circular."
        return fill(NaN, 6)
    end

    # Asymptote direction unit vector
    ŝ = [cos(δₐ) * cos(λₐ), cos(δₐ) * sin(λₐ), sin(δₐ)]

    # Define inertial unit vectors
    ẑ = [0.0, 0.0, 1.0]  # z-axis
    x̂ = [1.0, 0.0, 0.0]  # x-axis

    # Ensure asymptote vector is not aligned with z-axis
    if acos(abs(dot(ŝ, ẑ))) < tol
        @warn "Conversion failed: Asymptote vector is aligned with z-axis."
        return fill(NaN, 6)
    end

    # Build B-plane coordinate frame
    Ê = cross(ẑ, ŝ) / norm(cross(ẑ, ŝ))
    N̂ = cross(ŝ, Ê)

    # Angular momentum direction in inertial frame
    ami = π/2 - θᵦ
    ĥ = sin(ami) * Ê + cos(ami) * N̂

    # Inclination is angle between angular momentum and z-axis
    i = acos(clamp(dot(ẑ, ĥ), -1.0, 1.0))

    # Node vector points along line of nodes (intersection of orbit plane with equator)
    nodevec = cross(ẑ, ĥ)
    n = norm(nodevec)

    # Determine eccentricity direction unit vector
    if c₃ <= -tol
        # Elliptical orbit: Eccentricity vector opposite to asymptote
        ê = -ŝ
    else
        # Hyperbolic orbit: Compute ê from turning angle
        νₘ = acos(clamp(-1 / e, -1.0, 1.0))
        ô = cross(ĥ, ŝ)  # Orbit-normal vector
        ê = -sin(νₘ) * ô + cos(νₘ) * ŝ
    end

    # Compute Ω and ω based on inclination and eccentricity direction
    if e >= tol && i >= tol
        # General inclined case
        if n < tol
            @warn "Conversion failed: Line-of-nodes is undefined."
            return fill(NaN, 6)
        end

        # Ω from node vector projection on x-axis
        Ω = acos(clamp(dot(x̂, nodevec) / n, -1.0, 1.0))
        Ω = nodevec[2] < 0 ? 2π - Ω : Ω

        # ω from projection of eccentricity direction into orbital plane
        ω = acos(clamp(dot(nodevec / n, ê), -1.0, 1.0))
        ω = ê[3] < 0 ? 2π - ω : ω

    elseif e >= tol && i < tol
        # Equatorial prograde orbit
        Ω = 0.0
        ω = acos(clamp(ê[1], -1.0, 1.0))
        ω = ê[2] < 0 ? 2π - ω : ω

    elseif e >= tol && i ≥ π - tol
        # Equatorial retrograde orbit
        Ω = 0.0
        ω = -acos(clamp(ê[1], -1.0, 1.0))
        ω = ê[2] < 0 ? 2π - ω : ω

    end

    # Return Keplerian state
    return [a, e, i, Ω, ω, ν]
end
