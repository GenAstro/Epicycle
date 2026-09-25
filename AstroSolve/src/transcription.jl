# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A transcription discretizes a continuous arc for a numerical solver. This file declares what an
# implementation must provide and nothing else: twelve names, all public, so an implementation can
# live in any package on either side of the open/Enterprise boundary.
#
# That the interface says nothing about which tier implements it is the point. Hermite-Simpson is
# the open reference implementation and LGL is commercial, and moving either across that line is a
# packaging change rather than an interface change. See Architecture 11.1.1.

"""
    EvalContext(model, u, params)

The model, control, and phase parameters supplied to the dynamics at one node.

A transcription builds these internally, one per node, rather than being handed
one. It exists so a right-hand side written as `f!(dy, y, p, t)` can reach the
control and the model through `p` without the transcription and the dynamics
agreeing on a calling convention.

# Fields
- `model::M`: Model supplied when the phase was built.
- `u::U`: Control at this node.
- `params::P`: Phase parameters, empty when there are none.

# Example
```jldoctest
julia> using AstroSolve: EvalContext

julia> ctx = EvalContext(nothing, [0.3, -0.1], [3.986e5]);

julia> ctx.u                      # the control at this node
2-element Vector{Float64}:
  0.3
 -0.1

julia> EvalContext(nothing, [0.0]).params   # no parameters is the common case
Float64[]
```
"""
struct EvalContext{M, U<:AbstractVector, P<:AbstractVector}
    model::M
    u::U
    params::P
end

"""
    defect_row_groups(transcription, n_states) -> Vector{Tuple{UnitRange{Int}, Vector{Int}}}

Which nodes each group of defect rows reads, for declaring the constraint Jacobian's sparsity.

One entry per group: the rows it occupies inside the defect block, counted from one, and the node
indices whose state and control those rows depend on. A transcription that couples only neighbouring
nodes returns many small groups; one whose defects couple every node in an interval returns one group
per interval, which is the honest answer for a scheme built on an integration matrix.

The default is an empty vector, and a transcription that does not declare its stencil has its defect
rows declared dense. That is correct and costs a dense factorisation, so a transcription can be added
without this and gain it later.

Like the other extension points in this file, this one lets a transcription sit on either side of
the open/Enterprise boundary without an interface change.

# Arguments
- `transcription`: Transcription whose defect stencil is being declared.
- `n_states::Int`: Number of state components at each node.

# Returns
A vector of `(rows, nodes)` pairs, empty when the transcription declares no stencil.

# Example
```julia
defect_row_groups(HermiteSimpson(n_steps = 3), 2)
```
"""
defect_row_groups(::Any, ::Int) = Tuple{UnitRange{Int}, Vector{Int}}[]

"""
    has_defect_blocks(transcription) -> Bool

Whether this transcription can emit its defect Jacobian as sub-blocks rather than one matrix.

The companion to `defect_row_groups`: that one declares where the nonzeros are, this one says the
values can be produced without assembling the chunk. A transcription that returns `true` provides
`defect_blocks_Y` and `defect_blocks_U` methods on its own mesh type.

False by default, and the assembler then calls `jacobian_chunk`, which is correct and allocates a
matrix the size of the trajectory squared. So the two halves are gained independently: a transcription
can declare its stencil first and emit blocks later.

# Arguments
- `transcription`: Transcription to query for blockwise Jacobian support.

# Returns
`true` when the blockwise producers exist for this transcription.

# Example
```julia
has_defect_blocks(HermiteSimpson(n_steps = 3))
```
"""
has_defect_blocks(::Any) = false


# A phase with no parameters still needs a context, so `params` may be omitted.
# The three-argument form is the struct's own constructor.
EvalContext(model, u::AbstractVector) = EvalContext(model, u, Float64[])

# --- what a transcription implements -----------------------------------------

"""
    n_intervals(transcription) -> Int

How many mesh intervals the transcription divides an arc into.

# Arguments
- `transcription`: Transcription whose interval count is requested.

# Returns
The interval count, as an `Int`. One number for the whole arc; the mesh is not
built to answer this.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, n_intervals

julia> n_intervals(HermiteSimpson(n_steps = 20))
20
```
"""
function n_intervals end

"""
    n_unique_nodes(transcription) -> Int

Total distinct nodes across the mesh, midpoints included and shared interval
boundaries counted once.

# Arguments
- `transcription`: Transcription whose node count is requested.

# Returns
The node count, as an `Int`. This is the column count of the state matrix `Y`
the framework will hand back, so it fixes the size of the NLP.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, n_unique_nodes

julia> n_unique_nodes(HermiteSimpson(n_steps = 20))   # a midpoint per step, ends shared
41
```
"""
function n_unique_nodes end

"""
    n_defect_rows(transcription, n_states) -> Int

How many defect constraint rows the transcription contributes for a phase with
`n_states` states. LGL gives one equation per node; Hermite-Simpson gives two
per step.

# Arguments
- `transcription`: Transcription whose defect row count is requested.
- `n_states::Int`: Number of state components at each node.

# Returns
The row count, as an `Int`. The framework reserves exactly this many rows, so a
transcription that returns a count disagreeing with what `compute_defects`
produces corrupts every row after its own.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, n_defect_rows

julia> n_defect_rows(HermiteSimpson(n_steps = 20), 6)   # two per step, six states
240
```
"""
function n_defect_rows end

"""
    build_mesh(transcription) -> mesh

The mesh object the seven functions below are called on. Its type is the
transcription's own and the framework does not inspect it.

Called once when a phase is set up, so a transcription may put anything
expensive here — quadrature weights, a differentiation matrix — rather than
recomputing it each iteration.

# Arguments
- `transcription`: Transcription for which to build a mesh.

# Returns
Whatever the transcription uses to answer the seven mesh functions. The
framework only passes it back.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, build_mesh

julia> mesh = build_mesh(HermiteSimpson(n_steps = 4))
HermiteSimpsonMesh(n_steps=4, N=9)
```
"""
function build_mesh end

# --- what the mesh answers ---------------------------------------------------

"""
    node_times(mesh, t0, tf) -> Vector

Epoch at each node, for a phase spanning `t0` to `tf`.

# Arguments
- `mesh`: Mesh returned by [`build_mesh`](@ref).
- `t0`: Phase start epoch.
- `tf`: Phase end epoch.

# Returns
A vector of `n_unique_nodes` times, ascending, starting at `t0` and ending at
`tf`. These are the times the dynamics are evaluated at, so a path constraint
and a Lagrange integrand see the same ones.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, build_mesh, node_times

julia> node_times(build_mesh(HermiteSimpson(n_steps = 2)), 0.0, 1.0)
5-element Vector{Float64}:
 0.0
 0.25
 0.5
 0.75
 1.0
```
"""
function node_times end

"""
    quadrature_weights(mesh, t0, tf) -> Vector

Weights approximating an integral over the arc, for a Lagrange cost term.

The weights scale with the arc, so they carry the `tf - t0` factor and a caller
sums `w .* g` without rescaling.

# Arguments
- `mesh`: Mesh returned by [`build_mesh`](@ref).
- `t0`: Phase start epoch.
- `tf`: Phase end epoch.

# Returns
A vector of `n_unique_nodes` weights summing to `tf - t0`, which is the integral
of one over the arc and the cheapest check that a transcription got them right.

# Example
```jldoctest
julia> using AstroSolve: HermiteSimpson, build_mesh, quadrature_weights

julia> w = quadrature_weights(build_mesh(HermiteSimpson(n_steps = 4)), 0.0, 2.0);

julia> sum(w) ≈ 2.0               # the integral of 1 over the arc
true
```
"""
function quadrature_weights end

"""
    compute_defects(dynamics!, Y, U, params, model, mesh, t0, tf) -> Matrix

The defect residuals: how far the discretised trajectory is from satisfying the
equations of motion. Zero defects mean the discretisation is consistent.

`Y` holds the state at each node, one column per node, and `U` the control the
same way. `dynamics!` is called as `dynamics!(dy, y, ctx, t)` with an
[`EvalContext`](@ref) carrying the model, the control and the parameters.

# Arguments
- `dynamics!`: In-place dynamics function called as `dynamics!(dy, y, ctx, t)`.
- `Y`: State matrix with one node per column.
- `U`: Control matrix with one node per column.
- `params`: Phase parameters.
- `model`: Model passed to the dynamics through [`EvalContext`](@ref).
- `mesh`: Mesh returned by [`build_mesh`](@ref).
- `t0`: Phase start epoch.
- `tf`: Phase end epoch.

# Returns
A matrix the framework flattens into `n_defect_rows` constraint rows. Its shape
is the transcription's own; what matters is that the element count matches what
`n_defect_rows` promised.
"""
function compute_defects end

"""
    defect_jacobian_Y(dynamics_jac_y!, Y, U, params, model, mesh, t0, tf)

Derivative of the defects with respect to the node states. The remaining four
follow the same shape for the control, the two phase times, and the parameters.

`dynamics_jac_y!` supplies the derivative of the right-hand side, and the
transcription's job is to propagate that through its own discretisation.

# Arguments
- `dynamics_jac_y!`: In-place state Jacobian of the dynamics.
- `Y`: State matrix with one node per column.
- `U`: Control matrix with one node per column.
- `params`: Phase parameters.
- `model`: Model passed to the dynamics.
- `mesh`: Mesh returned by [`build_mesh`](@ref).
- `t0`: Phase start epoch.
- `tf`: Phase end epoch.

# Returns
A matrix with one row per defect row and one column per state variable, both in
the order the framework laid them out. A transcription that gets an entry wrong
is worse than one that omits a whole block: the framework merges columns when it
colours the Jacobian, so a wrong entry corrupts entries that were right.
"""
function defect_jacobian_Y end

"""
    defect_jacobian_U(args...)

Derivative of the defects with respect to the control. Same shape as
[`defect_jacobian_Y`](@ref).

# Arguments
- `args...`: Control-Jacobian function, phase data, mesh, and phase endpoints required by the
  transcription implementation.

# Returns
A matrix with one row per defect row and one column per component of the control.
"""
function defect_jacobian_U end

"""
    defect_jacobian_t0(args...)

Derivative of the defects with respect to the start time. Same shape as
[`defect_jacobian_Y`](@ref).

A non-autonomous right-hand side contributes a term here that an autonomous
one does not.

# Arguments
- `args...`: Dynamics partials, phase data, mesh, and phase endpoints required by the
  transcription implementation.

# Returns
A matrix with one row per defect row and one column per component of the start time.
"""
function defect_jacobian_t0 end

"""
    defect_jacobian_tf(args...)

Derivative of the defects with respect to the end time. Same shape as
[`defect_jacobian_Y`](@ref).

# Arguments
- `args...`: Dynamics partials, phase data, mesh, and phase endpoints required by the
  transcription implementation.

# Returns
A matrix with one row per defect row and one column per component of the end time.
"""
function defect_jacobian_tf end

"""
    defect_jacobian_P(args...)

Derivative of the defects with respect to the phase parameters. Same shape as
[`defect_jacobian_Y`](@ref).

# Arguments
- `args...`: Parameter-Jacobian function, phase data, mesh, and phase endpoints required by the
  transcription implementation.

# Returns
A matrix with one row per defect row and one column per component of the phase parameters.
"""
function defect_jacobian_P end
