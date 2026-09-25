# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Hermite-Simpson (Lobatto IIIA, 4th order, S=3) collocation transcription.
#
# The open tier's reference implementation of the transcription extension
# interface in transcription.jl. It implements the four methods on the
# transcription type and the seven on the mesh, and nothing else about it is
# known to the framework that drives it.
#
# Reference: Betts, "Practical Methods for Optimal Control", Table (S=3).
#            CSALT/GMAT: LobattoIIIA_4HSOrder, NLPFuncUtil_ImplicitRungeKutta.
#
# Grid layout (n_steps steps, N = 2*n_steps+1 total points):
#   k=1: [t_0,  t_{0.5},  t_1]
#   k=2: [t_1,  t_{1.5},  t_2]
#   ...
#   k=n: [t_{n-1}, t_{n-0.5}, t_n]
#
# Per-step defect equations (2 per step, from paramDepArray / funcConstArray):
#
#   r1 = (-½)y_L + y_M + (-½)y_R  +  h·[(-⅛)f_L + 0·f_M + (⅛)f_R]  = 0
#   r2 = (-1)y_L + 0·y_M + y_R    +  h·[(-⅙)f_L + (-⅔)f_M + (-⅙)f_R] = 0
#
# where h = (tf-t0)/n_steps (physical step size).
#
# r1 enforces Hermite midpoint interpolation: y_M = (y_L+y_R)/2 + h/8·(f_L-f_R)
# r2 enforces Simpson quadrature:             y_R - y_L = h/6·(f_L + 4f_M + f_R)
#
# Jacobians (analytic):
#   ∂r_i/∂y_p  =  HS_A[i,p]·I  +  h·HS_B[i,p]·∂f_p/∂y_p
#   ∂r_i/∂u_p  =  h·HS_B[i,p]·∂f_p/∂u_p
#   ∂r_i/∂tf   =  Σ_p HS_B[i,p]·f_p / n_steps  (autonomous)
#   ∂r_i/∂t0   = -Σ_p HS_B[i,p]·f_p / n_steps  (autonomous)
#

# ─────────────────────────────────────────────────────────────────────────────
# The transcription type, and the four methods the interface asks for
# ─────────────────────────────────────────────────────────────────────────────

"""
    HermiteSimpson(; n_steps = 10)

Hermite-Simpson collocation over `n_steps` mesh steps.

Each step carries three points, its two ends and a midpoint, and contributes two
defect equations per state: one enforcing Hermite midpoint interpolation and one
enforcing Simpson quadrature. Adjacent steps share an end point, so a phase of
`n` steps has `2n + 1` distinct points.

# Fields
- `n_steps::Int`: Number of mesh steps.

# Notes
A transcription with `n_steps` contains `2n_steps + 1` nodes.

# Examples
```julia
t    = HermiteSimpson(n_steps = 20)
mesh = build_mesh(t)
ts   = node_times(mesh, 0.0, 1.0)     # 41 points, midpoints included
```
"""
struct HermiteSimpson
    n_steps::Int
end

HermiteSimpson(; n_steps::Int = 10) = HermiteSimpson(n_steps)

n_intervals(t::HermiteSimpson)    = t.n_steps
n_unique_nodes(t::HermiteSimpson) = 2 * t.n_steps + 1
n_defect_rows(t::HermiteSimpson, ns::Int) = 2 * ns * t.n_steps
build_mesh(t::HermiteSimpson)     = HermiteSimpsonMesh(t)

has_defect_blocks(::HermiteSimpson) = true

"""Step k occupies `2·ns` rows and reads nodes 2k-1, 2k and 2k+1. See `defect_row_groups`."""
function defect_row_groups(t::HermiteSimpson, ns::Int)
    groups = Vector{Tuple{UnitRange{Int}, Vector{Int}}}(undef, t.n_steps)
    for k in 1:t.n_steps
        groups[k] = ((2 * (k - 1)) * ns + 1 : 2k * ns, [2k - 1, 2k, 2k + 1])
    end
    return groups
end

# ─────────────────────────────────────────────────────────────────────────────
# Butcher table constants (fixed; same for all HS problems)
# ─────────────────────────────────────────────────────────────────────────────

# paramDepArray:  A[i, p_local] — state coefficient for defect i, point p
const HS_A = [ -1/2   1.0  -1/2 ;   # defect 1: Hermite midpoint interpolation
               -1.0   0.0   1.0  ]   # defect 2: Simpson quadrature

# funcConstArray: B[i, p_local] — dynamics coefficient for defect i, point p
const HS_B = [ -1/8   0.0   1/8 ;   # defect 1
               -1/6  -2/3  -1/6  ]   # defect 2

# betaVec: per-step quadrature weights for Lagrange integrals
const HS_β = [1/6, 2/3, 1/6]

# ─────────────────────────────────────────────────────────────────────────────
# HermiteSimpsonMesh
# ─────────────────────────────────────────────────────────────────────────────

"""
    HermiteSimpsonMesh(transcription::HermiteSimpson)

The mesh a `HermiteSimpson` transcription works on, built once per phase.

Hermite-Simpson places a midpoint in every step and shares each step boundary
with its neighbour, so a mesh of `n` steps carries `2n + 1` nodes. The node
positions are held in normalised time and scaled to the phase's own span when
[`node_times`](@ref) is called, which is what lets `t0` and `tf` be solver
variables without rebuilding the mesh.

# Fields
- `n_steps::Int`: Number of mesh steps.
- `N::Int`: Total nodes, `2 * n_steps + 1`.
- `τ_global::Vector{Float64}`: Normalized node positions from zero to one.
- `h::Float64`: Normalized step size, `1 / n_steps`.

# Examples
```jldoctest
using AstroSolve: HermiteSimpson, HermiteSimpsonMesh, n_unique_nodes
mesh = HermiteSimpsonMesh(HermiteSimpson(n_steps = 3))
(mesh.N, mesh.N == n_unique_nodes(HermiteSimpson(n_steps = 3)))

# output
(7, true)
```
"""
struct HermiteSimpsonMesh
    n_steps  ::Int
    N        ::Int                   # total collocation points = 2*n_steps + 1
    τ_global ::Vector{Float64}       # normalized times τ ∈ [0,1], length N
    h        ::Float64               # normalized step size = 1/n_steps
end

function HermiteSimpsonMesh(hs::HermiteSimpson)
    n = hs.n_steps
    N = 2 * n + 1
    h = 1.0 / n
    τ = [(i - 1) / (2 * n) for i in 1:N]   # 0, 1/(2n), 2/(2n), ..., 1
    return HermiteSimpsonMesh(n, N, τ, h)
end

node_times(m::HermiteSimpsonMesh, t0::Float64, tf::Float64) =
    t0 .+ m.τ_global .* (tf - t0)

# quadrature_weights: composite Simpson rule weights for Lagrange integrals.
# Shared mesh endpoints (right of step k = left of step k+1) accumulate
# contributions from both adjacent steps.
function quadrature_weights(m::HermiteSimpsonMesh, t0::Float64, tf::Float64)
    h_phys = m.h * (tf - t0)
    q = zeros(m.N)
    for k in 1:m.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        q[L] += h_phys * HS_β[1]
        q[M] += h_phys * HS_β[2]
        q[R] += h_phys * HS_β[3]
    end
    return q
end

Base.show(io::IO, m::HermiteSimpsonMesh) =
    print(io, "HermiteSimpsonMesh(n_steps=$(m.n_steps), N=$(m.N))")

# ─────────────────────────────────────────────────────────────────────────────
# compute_defects
# ─────────────────────────────────────────────────────────────────────────────

function compute_defects(dynamics!::Function,
                         Y::Matrix{Float64},
                         U::Matrix{Float64},
                         params::Vector{Float64},
                         model,
                         mesh::HermiteSimpsonMesh,
                         t0::Float64,
                         tf::Float64)
    ns      = size(Y, 1)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dy      = zeros(ns)

    # Evaluate dynamics at all N collocation points
    F = zeros(ns, mesh.N)
    for k in 1:mesh.N
        ctx = EvalContext(model, U[:, k], params)
        dynamics!(dy, Y[:, k], ctx, t_nodes[k])
        F[:, k] .= dy
    end

    # Assemble defect residuals: 2 equations per step
    res = zeros(2 * ns * mesh.n_steps)
    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for i in 1:2
            rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
            res[rows] .= (HS_A[i,1] .* Y[:, L] .+ HS_A[i,2] .* Y[:, M] .+ HS_A[i,3] .* Y[:, R]
                       .+ h_phys .* (HS_B[i,1] .* F[:, L] .+ HS_B[i,2] .* F[:, M] .+ HS_B[i,3] .* F[:, R]))
        end
    end
    return res
end

# ─────────────────────────────────────────────────────────────────────────────
# Analytic defect Jacobians
# ─────────────────────────────────────────────────────────────────────────────

"""
    defect_jacobian_Y(jac_fn!, Y, U, params, model, mesh::HermiteSimpsonMesh, t0, tf)

Analytic Jacobian ∂(defects)/∂(vec(Y)).
Shape: (2·ns·n_steps) × (ns·N).

For each step k and each of the 3 points p ∈ {L,M,R}:
  block[d_i, p] += HS_A[i,p]·I  +  h·HS_B[i,p]·∂f_p/∂y_p
"""
function defect_jacobian_Y(dynamics_jac_y!::Function,
                            Y::Matrix{Float64}, U::Matrix{Float64},
                            params::Vector{Float64}, model,
                            mesh::HermiteSimpsonMesh,
                            t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    n_rows  = 2 * ns * mesh.n_steps
    n_cols  = ns * mesh.N
    J       = zeros(n_rows, n_cols)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dFdy    = zeros(ns, ns)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            fill!(dFdy, 0.0)
            ctx = EvalContext(model, U[:, gi], params)
            dynamics_jac_y!(dFdy, Y[:, gi], ctx, t_nodes[gi])

            y_col = (gi - 1) * ns + 1 : gi * ns   # column block for state at point gi
            for i in 1:2
                d_rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
                a = HS_A[i, p_local]
                b = HS_B[i, p_local]
                blk = view(J, d_rows, y_col)
                b != 0.0 && (blk .+= h_phys * b .* dFdy)
                if a != 0.0
                    for s in 1:ns
                        blk[s, s] += a
                    end
                end
            end
        end
    end
    return J
end

"""
    defect_blocks_Y(sink, jac_fn!, Y, U, params, model, mesh::HermiteSimpsonMesh, t0, tf)

Emit the nonzero sub-blocks of the state defect Jacobian instead of assembling it.

`defect_jacobian_Y` returns a `(2·ns·n_steps)` by `(ns·N)` matrix, which is the whole trajectory
squared and cannot exist for a large problem. The blocks it fills are `ns` by `ns` and there are
`6·n_steps` of them, so emitting them costs memory proportional to the mesh rather than its square.

`sink(block, row_offset, col_offset)` is called once per block, with offsets counted from the top
left of the chunk the assembler placed. **The block is a buffer reused between calls**, so a sink
reads it before returning and does not keep it.

# Returns
`nothing`. Everything is delivered through `sink`.
"""
function defect_blocks_Y(sink, dynamics_jac_y!::Function,
                         Y::Matrix{Float64}, U::Matrix{Float64},
                         params::Vector{Float64}, model,
                         mesh::HermiteSimpsonMesh, t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dFdy    = zeros(ns, ns)
    blk     = zeros(ns, ns)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            fill!(dFdy, 0.0)
            ctx = EvalContext(model, U[:, gi], params)
            dynamics_jac_y!(dFdy, Y[:, gi], ctx, t_nodes[gi])

            for i in 1:2
                a = HS_A[i, p_local]
                b = HS_B[i, p_local]
                (a == 0.0 && b == 0.0) && continue
                fill!(blk, 0.0)
                b != 0.0 && (blk .+= h_phys * b .* dFdy)
                if a != 0.0
                    for s in 1:ns
                        blk[s, s] += a
                    end
                end
                sink(blk, (2(k - 1) + i - 1) * ns, (gi - 1) * ns)
            end
        end
    end
    return nothing
end

"""
    defect_blocks_U(sink, jac_fn!, Y, U, params, model, mesh::HermiteSimpsonMesh, t0, tf)

Emit the nonzero sub-blocks of the control defect Jacobian instead of assembling it.

The counterpart to [`defect_blocks_Y`](@ref), with the same contract: `ns` by `nc` blocks through
`sink`, offsets from the top left of the chunk, and a buffer that is reused between calls. The first
defect contributes nothing at the midpoint, since `HS_B[1,2]` is zero, so that block is not emitted
at all rather than emitted as zeros.

# Returns
`nothing`. Everything is delivered through `sink`.
"""
function defect_blocks_U(sink, dynamics_jac_u!::Function,
                         Y::Matrix{Float64}, U::Matrix{Float64},
                         params::Vector{Float64}, model,
                         mesh::HermiteSimpsonMesh, t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    nc      = size(U, 1)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dFdu    = zeros(ns, nc)
    blk     = zeros(ns, nc)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            fill!(dFdu, 0.0)
            ctx = EvalContext(model, U[:, gi], params)
            dynamics_jac_u!(dFdu, Y[:, gi], ctx, t_nodes[gi])

            for i in 1:2
                b = HS_B[i, p_local]
                b == 0.0 && continue
                blk .= h_phys .* b .* dFdu
                sink(blk, (2(k - 1) + i - 1) * ns, (gi - 1) * nc)
            end
        end
    end
    return nothing
end

"""
    defect_jacobian_U(jac_fn!, Y, U, params, model, mesh::HermiteSimpsonMesh, t0, tf)

Analytic Jacobian ∂(defects)/∂(vec(U)).
Shape: (2·ns·n_steps) × (nc·N).

  block[d_i, p] += h·HS_B[i,p]·∂f_p/∂u_p
"""
function defect_jacobian_U(dynamics_jac_u!::Function,
                            Y::Matrix{Float64}, U::Matrix{Float64},
                            params::Vector{Float64}, model,
                            mesh::HermiteSimpsonMesh,
                            t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    nc      = size(U, 1)
    n_rows  = 2 * ns * mesh.n_steps
    n_cols  = nc * mesh.N
    J       = zeros(n_rows, n_cols)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dFdu    = zeros(ns, nc)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            fill!(dFdu, 0.0)
            ctx = EvalContext(model, U[:, gi], params)
            dynamics_jac_u!(dFdu, Y[:, gi], ctx, t_nodes[gi])

            u_col = (gi - 1) * nc + 1 : gi * nc
            for i in 1:2
                b = HS_B[i, p_local]
                b == 0.0 && continue   # skip zero contributions (e.g. defect 1 / midpoint)
                d_rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
                view(J, d_rows, u_col) .+= h_phys * b .* dFdu
            end
        end
    end
    return J
end

"""
    defect_jacobian_t0(Y, U, params, model, dynamics!, mesh::HermiteSimpsonMesh, t0, tf)

Jacobian ∂(defects)/∂t0 for HS transcription.  Shape: (2·ns·n_steps, 1).

For autonomous dynamics (∂F/∂t = 0):
  ∂r_i/∂t0 = -Σ_p HS_B[i,p]·f_p / n_steps

The non-autonomous correction is applied separately by _add_nonauto_time_correction!.
"""
function defect_jacobian_t0(Y::Matrix{Float64}, U::Matrix{Float64},
                             params::Vector{Float64}, model, dynamics!,
                             mesh::HermiteSimpsonMesh,
                             t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    n_rows  = 2 * ns * mesh.n_steps
    J       = zeros(n_rows, 1)
    t_nodes = node_times(mesh, t0, tf)
    dy      = zeros(ns)

    # Evaluate F at all points (needed for the h-scaling term)
    F = zeros(ns, mesh.N)
    for k in 1:mesh.N
        ctx = EvalContext(model, U[:, k], params)
        dynamics!(dy, Y[:, k], ctx, t_nodes[k])
        F[:, k] .= dy
    end

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for i in 1:2
            d_rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
            # ∂h/∂t0 = -1/n_steps → contribution = -B*f/n_steps
            contrib = zeros(ns)
            for (p_local, gi) in enumerate((L, M, R))
                b = HS_B[i, p_local]
                b == 0.0 && continue
                contrib .+= b .* F[:, gi]
            end
            J[d_rows, 1] .= .-contrib ./ mesh.n_steps
        end
    end
    return J
end

"""
    defect_jacobian_tf(Y, U, params, model, dynamics!, mesh::HermiteSimpsonMesh, t0, tf)

Jacobian ∂(defects)/∂tf for HS transcription.  Shape: (2·ns·n_steps, 1).
Negation of defect_jacobian_t0 for autonomous dynamics.
"""
function defect_jacobian_tf(Y::Matrix{Float64}, U::Matrix{Float64},
                             params::Vector{Float64}, model, dynamics!,
                             mesh::HermiteSimpsonMesh,
                             t0::Float64, tf::Float64)
    J = defect_jacobian_t0(Y, U, params, model, dynamics!, mesh, t0, tf)
    J .= .-J
    return J
end

"""
    defect_jacobian_P(jac_fn!, Y, U, params, model, mesh::HermiteSimpsonMesh, t0, tf)

Analytic Jacobian ∂(defects)/∂(params).  Shape: (2·ns·n_steps) × np.

  block[d_i, :] -= h·HS_B[i,p]·∂f_p/∂params  (accumulated over all points)
"""
function defect_jacobian_P(dynamics_jac_p!::Function,
                            Y::Matrix{Float64}, U::Matrix{Float64},
                            params::Vector{Float64}, model,
                            mesh::HermiteSimpsonMesh,
                            t0::Float64, tf::Float64)
    ns      = size(Y, 1)
    np      = length(params)
    n_rows  = 2 * ns * mesh.n_steps
    J       = zeros(n_rows, np)
    h_phys  = mesh.h * (tf - t0)
    t_nodes = node_times(mesh, t0, tf)
    dFdp    = zeros(ns, np)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            fill!(dFdp, 0.0)
            ctx = EvalContext(model, U[:, gi], params)
            dynamics_jac_p!(dFdp, Y[:, gi], ctx, t_nodes[gi])

            for i in 1:2
                b = HS_B[i, p_local]
                b == 0.0 && continue
                d_rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
                J[d_rows, :] .+= h_phys * b .* dFdp
            end
        end
    end
    return J
end

# 🔴 The non-autonomous time correction for Hermite-Simpson lives outside this package.
# defect_jacobian_t0 and defect_jacobian_tf below return the autonomous answer,
# and the framework patches in the ∂F/∂t term afterwards. That split is a leak in
# the interface: an interface method returns an incomplete result and the caller
# is expected to know it. A third-party transcription could not know it. Fixing
# it means every transcription returning its complete Jacobian, which touches LGL
# as well, so it is recorded rather than done here.
