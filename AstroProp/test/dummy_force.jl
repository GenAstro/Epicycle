

import AstroProp: accel_eval!

"Symbolic marker for thrust input (3D control)"
struct DummyThrust <: AbstractControl
    numvars::Int
    DummyThrust() = new(3)
end

"Symbolic marker for scalar parameter"
struct DummyParam <: AbstractParam
    numvars::Int
    DummyParam() = new(1)
end

"Symbolic marker for scalar mass state"
struct MassVal <: AbstractState
    numvars::Int
    MassVal() = new(1)
end

"""
    TheForceIsWithYou

Fake force model for testing STM superposition and Jacobian block types.
Includes:
- Constant acceleration normal to orbital plane.
- Constant uniform acceleration in all directions.
- Dependencies on state, control, and param types.

# Dependencies
- `PosVel()`: position and velocity
- `MassVal()`: scalar mass
- `DummyThrust()`: 3-element control
- `DummyParam()`: scalar parameter
"""
struct TheForceIsWithYou <: OrbitODE
    dependencies::Tuple{AbstractVar,AbstractVar,AbstractVar,AbstractVar}
    num_funs::Int
    params::Dict{AbstractVar, Any}  # Add a params field

    function TheForceIsWithYou(params=Dict(
        DummyParam() => 3e-8,
        DummyThrust() => -1e-8 .* ones(3)
    ))
        deps = (PosVel(), MassVal(), DummyThrust(), DummyParam())
        return new(deps, 6, params)
    end
end

function compute_the_force!(
    t::Time,
    posvel::AbstractVector{T},
    x̄̇::AbstractVector{T},
    params;  
    jac::Dict = Dict(),
) where T

    r̄ = posvel[1:3]
    v̄ = posvel[4:6]
    h̄ = cross(r̄, v̄)
    h = norm(h̄)
    ĥ = h̄ / h

    # Extract param and control
    p1 = params[DummyParam()]
    p2 = params[DummyThrust()]

    a_total = p1 * ĥ + p2

    x̄̇[1:3] = v̄
    x̄̇[4:6] = a_total

    # Jacobians
    if !isempty(jac)
        I3 = Matrix{T}(I, 3, 3)
        for k in keys(jac)
            if k === PosVel
                ∂a_∂r = zeros(T, 3, 3)
                ∂a_∂v = zeros(T, 3, 3)
                

                # h̄ = cross(r, v)
                # ∂h̄/∂r = [0  v3 -v2; -v3 0 v1; v2 -v1 0]
                # ∂h̄/∂v = [0 -r3 r2; r3 0 -r1; -r2 r1 0]
                ∂h̄_∂r = [
                    0        v̄[3]    -v̄[2];
                   -v̄[3]    0         v̄[1];
                    v̄[2]   -v̄[1]     0
                ]
                ∂h̄_∂v = [
                    0       -r̄[3]     r̄[2];
                    r̄[3]    0        -r̄[1];
                   -r̄[2]    r̄[1]     0
                ]

                # ∂(ĥ) = (I - ĥĥᵀ)/h ⋅ ∂h̄
                proj = I - ĥ * ĥ'
                ∂ĥ_∂r = (proj / h) * ∂h̄_∂r
                ∂ĥ_∂v = (proj / h) * ∂h̄_∂v

                ∂a_∂r .= p1 .* ∂ĥ_∂r
                ∂a_∂v .= p1 .* ∂ĥ_∂v

                jac[k][1:3, 4:6] .= I3  # ∂ẋ/∂v = I
                jac[k][4:6, 1:3] .= ∂a_∂r
                jac[k][4:6, 4:6] .= ∂a_∂v

            elseif k === DummyThrust
                jac[k][4:6,1:3] .= I3
            elseif k === DummyParam
                jac[k][4:6, 1] .= ĥ
            else
                error("Jacobian requested for unknown variable type: $k")
            end
        end
    end

    return nothing
end

function accel_eval!(force::TheForceIsWithYou, t::Time, x̄::Vector{Float64},
    x̄̇::Vector, sc::Spacecraft, params; jac::Dict = Dict()) 
    compute_the_force!(t, x̄, x̄̇, force.params; jac)
    return x̄̇, jac
end

