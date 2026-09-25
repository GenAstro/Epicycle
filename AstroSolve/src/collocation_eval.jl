# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Guess interpolation and the collocation evaluation path every transcription uses.
#

# Linear interpolation of a guess onto the mesh. Generic, not LGL's.
function _interp1(xs::AbstractVector, ys::AbstractVector, x::Float64)
    x <= xs[1]   && return ys[1]
    x >= xs[end] && return ys[end]
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i+1] - xs[i])
    return ys[i] * (1 - t) + ys[i+1] * t
end

# =============================================================================
# Collocation phase evaluation
#
# These two lived with LGL and were written untyped on the phase, which
# made them look like LGL's while being the machinery every collocation phase
# uses. They call `build_mesh` on whatever transcription the phase holds, so
# they were always transcription-agnostic — only their location said otherwise.
#
# That mattered: dropping LGL took the collocation evaluation path with
# it, and Hermite-Simpson could not solve. A boundary cannot pass through a file
# that holds both a transcription and the manager.
# =============================================================================

#
# Returns the flat function vector the NLP solver evaluates:
#   [defect residuals; boundary constraint residuals]
#
# Order aligns with get_constraint_bounds(phase).
# Defects are zero at a feasible trajectory; boundary residuals are zero when
# their Constraint lower_bounds ≤ value ≤ upper_bounds.
#
# Requires: phase._Y, phase._U, phase._t0, phase._tf populated via
#           set_decision_vector! before calling.
# ─────────────────────────────────────────────────────────────────────────────

function get_functions(phase::CollocationPhase)
    mesh = build_mesh(phase.transcription)

    defects = if size(phase._Y) == (phase._n_states, mesh.N)
        compute_defects(phase.dynamics, phase._Y, phase._U,
                        phase._params, phase.model, mesh, phase._t0, phase._tf)
    else
        zeros(n_defect_rows(phase.transcription, phase._n_states))
    end

    # Only evaluate boundary/path closures if state has been initialized
    bc_vals = Float64[]
    pc_vals = Float64[]
    if !isempty(phase._y0)
        for c in phase.constraints
            if c isa BoundaryConstraint
                ctx = BoundaryContext(_state_named(phase, phase._Y[:, 1]),
                                      _state_named(phase, phase._Y[:, end]),
                                      phase._t0, phase._tf, phase._params)
                # Support both new ctx-arg API and legacy zero-arg closures
                vals = applicable(c.calc.fn, ctx) ? c.calc.fn(ctx) : c.calc.fn()
                append!(bc_vals, vals)
            end
        end
        N = mesh.N
        t_nodes = node_times(mesh, phase._t0, phase._tf)
        for pc in phase.path_constraints
            if pc isa PathConstraint
                g_k = zeros(pc.n_pc)
                for k in 1:N
                    fill!(g_k, 0.0)
                    ctx = EvalContext(phase.model, phase._U[:, k], phase._params)
                    pc.fn(g_k, phase._Y[:, k], ctx, t_nodes[k])
                    append!(pc_vals, g_k)
                end
            end
        end
    end

    return [defects; bc_vals; pc_vals]
end

#   phase._Y / _U        — dense decision vector at LGL nodes
#
# Interpolates waypoint arrays linearly onto the LGL node grid, then
# populates _t0, _tf, _Y, _U, _y0, _yf.  After this call,
# get_decision_vector, get_functions, and all accessors work correctly.
#
# Conventions assumed about SolverVariable.var:
#   time variables    — .value::Float64
#   state variable    — .data::Matrix{Float64}  (n_states  × n_waypoints)
#   control variable  — .data::Matrix{Float64}  (n_controls × n_waypoints)
# ─────────────────────────────────────────────────────────────────────────────

function initialize!(phase::CollocationPhase)
    phase._t0 = phase.t0_var.var.value
    phase._tf = phase.tf_var.var.value

    mesh = build_mesh(phase.transcription)
    τ    = mesh.τ_global   # normalized ∈ [0,1], length N
    N    = mesh.N

    # Interpolate state guess (n_states × n_waypoints) → (n_states × N)
    Y_wp  = phase.state_var.var.data
    ns    = size(Y_wp, 1)
    n_wp  = size(Y_wp, 2)
    τ_wp  = collect(LinRange(0.0, 1.0, n_wp))
    phase._Y = Matrix{Float64}(undef, ns, N)
    for s in 1:ns
        phase._Y[s, :] = [_interp1(τ_wp, Y_wp[s, :], τk) for τk in τ]
    end

    # Interpolate control guess (n_controls × n_waypoints) → (n_controls × N)
    if phase.control_var !== nothing
        U_wp  = phase.control_var.var.data
        nc    = size(U_wp, 1)
        n_uwp = size(U_wp, 2)
        τ_uwp = collect(LinRange(0.0, 1.0, n_uwp))
        phase._U = Matrix{Float64}(undef, nc, N)
        for c in 1:nc
            phase._U[c, :] = [_interp1(τ_uwp, U_wp[c, :], τk) for τk in τ]
        end
    else
        phase._U = Matrix{Float64}(undef, 0, N)
    end

    # Initialize parameters (if present)
    if phase.param_var !== nothing
        phase._params = copy(vec(phase.param_var.var.data))
    end

    phase._y0 = phase._Y[:, 1]
    phase._yf = phase._Y[:, end]
    return nothing
end
