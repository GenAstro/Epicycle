# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

"""
    STMConfig

Declares which integrated sensitivity quantities to propagate alongside the
nominal trajectory.

# Fields
- `Φ::Bool` — propagate the 6×6 state transition matrix Φ = ∂y(t_f)/∂y(t_0)
- `S_p::Vector{ModelVariable}` — propagate ∂y(t_f)/∂p_i for each listed variable

# Augmented ODE
The propagated state is `z = [y(6); vec(Φ)(36); S_p_1(6); ...; S_p_n(6)]`
with augmented dynamics:
    ẏ       = f(y, t)
    Φ̇       = A·Φ           where A = ∂f/∂y
    Ṡ_p_i   = A·S_p_i + B_i  where B_i = ∂f/∂p_i

A and B_i are evaluated at each ODE step using the `JacobianConfig` machinery.
"""
struct STMConfig
    Φ::Bool
    S_p::Vector{ModelVariable}

    function STMConfig(; Φ::Bool = false, S_p = [])
        return new(Φ, collect(ModelVariable, S_p))
    end
end

"""
    OrbitODEProblem

Fully specified orbital propagation problem, optionally with STM / parameter
sensitivity propagation.

# Fields
- `prop::OrbitPropagator`
- `sc::Spacecraft`
- `duration_s::Float64`
- `stm::Union{STMConfig, Nothing}`
- `dense::Bool` — retain the integrator's dense interpolant in
  `PropagationResult.sol` (`save_everystep = dense`,
  `dense = dense`).  Default `false` for backward compatibility.
  Used by callers (e.g. measurement simulators) that need to query
  the trajectory at arbitrary epochs after the solve.

# Example
```julia
prob   = OrbitODEProblem(prop, sc, duration_s=5400.0, stm=stm_cfg)
result = solve(prob)
```
"""
struct OrbitODEProblem
    prop::OrbitPropagator
    sc::Spacecraft
    duration_s::Float64
    stm::Union{STMConfig, Nothing}
    dense::Bool

    function OrbitODEProblem(prop::OrbitPropagator, sc::Spacecraft;
                              duration_s::Real, stm = nothing,
                              dense::Bool = false)
        return new(prop, sc, Float64(duration_s), stm, dense)
    end
end

"""
    PropagationResult

Results of a `solve(OrbitODEProblem(...))` call.

# Fields
- `y_final::Vector{Float64}` — state at t_f
- `Φ::Union{Matrix{Float64}, Nothing}` — STM (nothing if not requested)
- `S_p::Dict{ModelVariable, Vector{Float64}}` — parameter sensitivities
- `sol::Any` — raw `ODESolution` if `prob.dense == true`, otherwise
  `nothing`.  Typed `Any` to avoid forcing every caller to depend on
  the OrdinaryDiffEq concrete solution type.
"""
struct PropagationResult{T<:Real}
    y_final::Vector{T}
    Φ::Union{Matrix{T}, Nothing}
    S_p::Dict{ModelVariable, Vector{T}}
    sol::Any
end

# Backward-compatible 3-arg constructor (no dense solution).
PropagationResult(y_final::Vector{T},
                  Φ::Union{Matrix{T}, Nothing},
                  S_p::Dict{ModelVariable, Vector{T}}) where {T<:Real} =
    PropagationResult{T}(y_final, Φ, S_p, nothing)

# ---------------------------------------------------------------------------
# solve  (extends OrdinaryDiffEq.solve via the import in AstroProp.jl)
# ---------------------------------------------------------------------------

"""
    solve(prob::OrbitODEProblem) -> PropagationResult

Integrate the orbit ODE and, if `prob.stm` is set, the augmented sensitivity
equations alongside it.

The plain-propagation path (no STM) is a thin wrapper around OrdinaryDiffEq.
The augmented path propagates `z = [y; vec(Φ); S_p_1; ...; S_p_n]` with a
single ODE solve; A and B_i blocks are recomputed at each step via
`eval_jacobian!` using the pre-allocated `JacobianResult`.
"""
function solve(prob::OrbitODEProblem)
    sc          = prob.sc
    prop        = prob.prop
    y0          = collect(sc.state.state)
    tspan       = (0.0, prob.duration_s)
    start_epoch = sc.time
    integ       = prop.integ

    if prob.stm === nothing
        # ------------------------------------------------------------------
        # Plain propagation
        # ------------------------------------------------------------------
        function plain_rhs!(dy, y, _p, t_rel)
            t = start_epoch + t_rel / 86400.0
            fill!(dy, 0.0)
            for force in prop.forces.forces
                accel_eval!(force, t, y, dy, sc, [])
            end
        end

        ode = ODEProblem(plain_rhs!, y0, tspan)
        sol = OrdinaryDiffEq.solve(ode, integ.integrator;
                                    reltol = integ.reltol,
                                    abstol = integ.abstol,
                                    save_everystep = prob.dense,
                                    dense          = prob.dense)
        u = sol.u[end]
        return PropagationResult{eltype(u)}(u, nothing,
                                            Dict{ModelVariable, Vector{eltype(u)}}(),
                                            prob.dense ? sol : nothing)
    else
        return _augmented_solve(prob, y0, tspan, start_epoch)
    end
end

function _augmented_solve(prob::OrbitODEProblem, y0, tspan, start_epoch)
    stm    = prob.stm
    prop   = prob.prop
    sc     = prob.sc
    forces = prop.forces
    integ  = prop.integ
    n_y    = length(y0)          # 6
    n_p    = length(stm.S_p)
    n_Φ    = stm.Φ ? n_y^2 : 0
    n_aug  = n_y + n_Φ + n_y * n_p

    # Pre-allocate Jacobian result — reused at every RHS call, zero allocation
    jac_cfg    = JacobianConfig(partial_y = stm.Φ || n_p > 0,
                                 partial_p = stm.S_p)
    jac_result = JacobianResult(jac_cfg; n_state = n_y,
                                 forces = forces, sc = sc,
                                 t_example = start_epoch)

    # Augmented initial condition: z = [y0; vec(I₆); zeros(6*n_p)]
    z0 = zeros(n_aug)
    z0[1:n_y] .= y0
    if stm.Φ
        for i in 1:n_y
            z0[n_y + (i-1)*n_y + i] = 1.0   # Φ(t0) = I, column-major
        end
    end

    function augmented_rhs!(dz, z, _p, t_rel)
        t = start_epoch + t_rel / 86400.0
        y = @view z[1:n_y]

        # Nominal dynamics
        dy = @view dz[1:n_y]
        fill!(dy, 0.0)
        for force in forces.forces
            accel_eval!(force, t, y, dy, sc, [])
        end

        # Jacobian blocks A = ∂f/∂y  and  B_i = ∂f/∂p_i
        eval_jacobian!(jac_result, forces, y, sc, t)
        A = jac_result.partial_y   # 6×6

        # Φ̇ = A·Φ
        if stm.Φ
            Φ_view  = reshape(@view(z[(n_y+1):(n_y+n_Φ)]),  n_y, n_y)
            dΦ_view = reshape(@view(dz[(n_y+1):(n_y+n_Φ)]), n_y, n_y)
            mul!(dΦ_view, A, Φ_view)
        end

        # Ṡ_p_i = A·S_p_i + B_i
        base = n_y + n_Φ
        for (i, mv) in enumerate(stm.S_p)
            idx  = (base + (i-1)*n_y + 1):(base + i*n_y)
            S_i  = @view  z[idx]
            dS_i = @view dz[idx]
            B_i  = jac_result.partial_p[mv]
            mul!(dS_i, A, S_i)
            dS_i .+= B_i
        end

        return nothing
    end

    ode = ODEProblem(augmented_rhs!, z0, tspan)
    sol = OrdinaryDiffEq.solve(ode, integ.integrator;
                                reltol = integ.reltol,
                                abstol = integ.abstol,
                                save_everystep = prob.dense,
                                dense          = prob.dense)

    z_f     = sol.u[end]
    T       = eltype(z_f)
    y_final = z_f[1:n_y]

    Φ_result = stm.Φ ?
        reshape(z_f[(n_y+1):(n_y+n_Φ)], n_y, n_y) :
        nothing

    S_p_result = Dict{ModelVariable, Vector{T}}()
    base = n_y + n_Φ
    for (i, mv) in enumerate(stm.S_p)
        idx = (base + (i-1)*n_y + 1):(base + i*n_y)
        S_p_result[mv] = z_f[idx]
    end

    return PropagationResult{T}(y_final, Φ_result, S_p_result,
                                prob.dense ? sol : nothing)
end
