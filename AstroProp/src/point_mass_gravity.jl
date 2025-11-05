# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    PointMassGravity

Combined force model for central body gravity and N-body point-mass perturbations.

# Fields
- `central_body::CelestialBody`: The central body about which dynamics are referenced.
- `perturbers::Tuple{Vararg{CelestialBody}}`: Other celestial bodies treated as point-mass perturbers.
- `dependencies::Vector{Type{<:AbstractVar}}`: Vector of variable dependencies (e.g., `PosVel`)
- `num_funs::Int`: Number of functions in the ODE
"""
struct PointMassGravity <: OrbitODE
    central_body::CelestialBody
    pert_bodies::Tuple{Vararg{CelestialBody}}
    dependencies::Vector{Type{<:AbstractVar}}
    num_funs::Int

    function PointMassGravity(central_body::CelestialBody, perturbers::Tuple{Vararg{CelestialBody}})
        all_bodies = (central_body, perturbers...)
        check_duplicates(all_bodies)
        return new(central_body, perturbers, [PosVel], 6)
    end
end

"""
    nbody_perts(t::Time, center::CelestialBody, pert_bodies::Tuple{Vararg{CelestialBody}}; jac::Bool=false)

Compute the gravitational acceleration on a central body due to a tuple of perturbing bodies using Newtonian point-mass gravity.

# Arguments
- `t::Time`: State epoch
- `posvel::Vector`  Orbit state (position and velocity) 
- `center::CelestialBody`: The central body of propagation
- `pert_bodies::Tuple{Vararg{CelestialBody}}`: Tuple of perturbing celestial bodies (e.g., Moon, Sun)
- `jac::Bool`: Optional keyword (default = `false`) to return the Jacobian of the perturbing acceleration with respect to the central body's position
- `tol::Real`: Optional tolerance on singularity testing (default 1e-12)
# Returns
- If `jac == false`: `a_pert::Vector{Float64}` — total perturbing acceleration
- If `jac == true`: `(a_pert::Vector{Float64}, jacobian::Matrix{Float64})`

# Notes
- Requires `AstroUniverse.translate(from::CelestialBody, to::CelestialBody, jd_tdb::Float64)` to return the vector from `from` to `to` in inertial coordinates.
- Units must be consistent with the gravitational parameters (`mu`) of the celestial bodies.
"""
function compute_point_mass_gravity!(
    t::Time,
    posvel::AbstractVector{T},
    x̄̇::AbstractVector{T},
    center::CelestialBody,
    pert_bodies::Tuple{Vararg{CelestialBody}};
    jac::Dict = Dict(),
    include_center::Bool = true,
    tol::Real = 1e-12,
) where T
    # TODO: Trap jacobian calls that don't exist
    # Time and state preparations (ephemeris is computed in TDB)
    t_tdb = t.tdb
    jd_tdb = t_tdb.jd
    r̄ = posvel[1:3]
    r = norm(r̄)

    if r < tol
        error("Computation of acceleration failed: Position is less than tol and approaching singularity.")
    end

    # Initialize Jacobian for PosVel if requested
    compute_jac_posvel = haskey(jac, PosVel)
    if compute_jac_posvel
        I3 = Matrix{T}(I, 3, 3)
        jac[PosVel][1:3, 4:6] .= I3
        ∂r̄̈∂r̄ = include_center ? center.mu * (3 * (r̄ * r̄') / r^5 - I3 / r^3) : zeros(T, 3, 3)
    end

    # Compute the acceleration of central body if requested
    acc = include_center ? -center.mu * r̄ / r^3 : zeros(T, 3)

    # Compute the acceleration and jacobian of perturbing bodies
    for pert in pert_bodies
        r̄ₖ = translate(center, pert, jd_tdb)
        r̄ᵣ = r̄ₖ - r̄
        rᵣ = norm(r̄ᵣ)
        if rᵣ < tol
            error("Computation of acceleration failed: Perturbing body vector is less than tol and approaching singularity.")
        end

        acc += pert.mu * (r̄ᵣ / rᵣ^3 - r̄ₖ / norm(r̄ₖ)^3)

        if compute_jac_posvel
            ∂r̄̈∂r̄ += pert.mu * (-I3 / rᵣ^3 + 3 * (r̄ᵣ * r̄ᵣ') / rᵣ^5)
        end
    end

    x̄̇[1:3] = posvel[4:6]
    x̄̇[4:6] = acc

    if compute_jac_posvel
        jac[PosVel][4:6, 1:3] .= ∂r̄̈∂r̄
    end

    return nothing
end

""" 
    accel_eval!(model::PointMassGravity, t::Time, x̄::Vector, 
                 x̄̇::Vector, sc::Spacecraft, params; jac::Dict = Dict())

Evaluate the acceleration due to point-mass gravity from central and perturbing bodies.
"""
function accel_eval!(model::PointMassGravity, t::Time, x̄::Vector, 
                        x̄̇::Vector, sc::Spacecraft, params; jac::Dict = Dict())
    compute_point_mass_gravity!(t,x̄, x̄̇, model.central_body, model.pert_bodies; jac = jac) 
    return x̄̇ 
end

"""
    function check_duplicates(bodies::Tuple{Vararg{CelestialBody}})
    
Validate that all names in force model are unique
"""
function check_duplicates(bodies::Tuple{Vararg{CelestialBody}})
    seen = Dict{String, Int}()
    for b in bodies
        name = b.name
        seen[name] = get(seen, name, 0) + 1
        if seen[name] > 1
            error("The CelestialBody $name is included in force model multiple times.")
        end
    end
end
