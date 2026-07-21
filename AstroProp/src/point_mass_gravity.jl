# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    PointMassGravity

Combined force model for central body gravity and N-body point-mass perturbations.

# Fields
- `central_body::CelestialBody`: The central body about which dynamics are referenced.
- `perturbers::Tuple{Vararg{CelestialBody}}`: Other celestial bodies treated as point-mass perturbers.
- `dependencies::Vector{Type{<:AbstractVar}}`: Legacy descriptor of variable dependencies; not read by any caller and slated for removal.
- `num_funs::Int`: Legacy state-dim hint; not read by any caller and slated for removal.
"""
struct PointMassGravity <: AbstractGravityForce
    central_body::CelestialBody
    pert_bodies::Tuple{Vararg{CelestialBody}}
    dependencies::Vector{DataType}
    num_funs::Int

    function PointMassGravity(central_body::CelestialBody, perturbers::Tuple{Vararg{CelestialBody}})
        all_bodies = (central_body, perturbers...)
        check_duplicates(all_bodies)
        return new(central_body, perturbers, DataType[PosVel], 6)
    end
end

# ---------------------------------------------------------------------------
# Kernel — generic in mu for ForwardDiff
# ---------------------------------------------------------------------------

"""
    _gravity_accel(mu, r) -> Vector

Point-mass gravitational acceleration kernel.  Generic in `mu::T` so that
`ForwardDiff.derivative` can seed mu with a Dual number for exact ∂f/∂μ.

a = -μ · r / ‖r‖³
"""
@inline function _gravity_accel(mu::T, r̄::AbstractVector) where T
    return -mu .* r̄ ./ norm(r̄)^3
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
    include_center::Bool = true,
    tol::Real = 1e-12,
) where T
    t_tdb  = t.tdb
    jd_tdb = t_tdb.jd
    r̄ = posvel[1:3]
    r = norm(r̄)

    if r < tol
        error("Computation of acceleration failed: Position is less than tol and approaching singularity.")
    end

    acc = include_center ? _gravity_accel(center.mu, r̄) : zeros(T, 3)

    for pert in pert_bodies
        r̄ₖ = translate(center, pert, jd_tdb)
        r̄ᵣ = r̄ₖ - r̄
        rᵣ = norm(r̄ᵣ)
        if rᵣ < tol
            error("Computation of acceleration failed: Perturbing body vector is less than tol and approaching singularity.")
        end
        acc += pert.mu * (r̄ᵣ / rᵣ^3 - r̄ₖ / norm(r̄ₖ)^3)
    end

    x̄̇[1:3] = posvel[4:6]
    x̄̇[4:6] = acc
    return nothing
end

"""
    accel_eval!(model::PointMassGravity, t::Time, x̄::AbstractVector,
                 x̄̇::AbstractVector, sc::Spacecraft, params)

Evaluate the acceleration due to point-mass gravity from central and perturbing bodies.
"""
function accel_eval!(model::PointMassGravity, t::Time, x̄::AbstractVector,
                        x̄̇::AbstractVector, sc::Spacecraft, params)
    compute_point_mass_gravity!(t, x̄, x̄̇, model.central_body, model.pert_bodies)
    return x̄̇
end

# ---------------------------------------------------------------------------
# state_jac! — A = ∂f/∂y  (analytic registration via dispatch)
# ---------------------------------------------------------------------------

function state_jac!(out::AbstractMatrix, m::PointMassGravity, t::Time,
                     y::AbstractVector, sc::Spacecraft)
    r̄ = y[1:3]
    r = norm(r̄)
    I3 = Matrix{Float64}(I, 3, 3)

    # ∂ṙ/∂v = I₃
    out[1:3, 4:6] .+= I3

    # ∂v̇/∂r = μ·(3r̂r̂ᵀ − I)/r³  for the central body
    ∂v̇∂r = m.central_body.mu * (3 * (r̄ * r̄') / r^5 - I3 / r^3)

    # Perturber contributions
    t_tdb  = t.tdb
    jd_tdb = t_tdb.jd
    for pert in m.pert_bodies
        r̄ₖ = translate(m.central_body, pert, jd_tdb)
        r̄ᵣ = r̄ₖ .- r̄
        rᵣ = norm(r̄ᵣ)
        ∂v̇∂r .+= pert.mu * (-I3 / rᵣ^3 + 3 * (r̄ᵣ * r̄ᵣ') / rᵣ^5)
    end

    out[4:6, 1:3] .+= ∂v̇∂r
    return nothing
end

# ---------------------------------------------------------------------------
# param_jac! — B = ∂f/∂p  (analytic registration via dispatch)
# ---------------------------------------------------------------------------

function param_jac!(out::AbstractVector, m::PointMassGravity, ::Mu, t::Time,
                     y::AbstractVector, sc::Spacecraft)
    r̄ = y[1:3]
    out[4:6] .+= ForwardDiff.derivative(
        μ -> _gravity_accel(μ, r̄), m.central_body.mu)
    return nothing
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
