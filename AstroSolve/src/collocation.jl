# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Registration, sizing, decision-vector packing, sparsity and the Jacobian chunks.
#

# ─────────────────────────────────────────────────────────────────────────────
# Registration functions
# ─────────────────────────────────────────────────────────────────────────────

# set_dynamics! — supports the domain-native mutating convention
# f!(dy, y_struct, u_struct, p, t, model) — for new API (model passed explicitly)
# f!(dy, y_vec, ctx, t)                   — for old API (ctx carries model)
const _raw_dynamics = Dict{UInt64, Any}()

function set_dynamics!(phase::CollocationPhase, f::Function; model=nothing)
    phase.model = model
    _raw_dynamics[objectid(phase)] = f
    ph = phase   # capture phase identity (not mutable ref) for registry lookup
    phase.dynamics = function _dynamics_wrapper!(dy, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            # `ctx.params` and not `nothing`. The p slot of a six-argument
            # dynamics function was always nothing, so a phase could declare
            # static parameters and its own dynamics could not read them.
            if reg.control_type !== nothing
                f(dy, y_s, u_s, ctx.params, t, ph.model)
            else
                f(dy, y_s, nothing, ctx.params, t, ph.model)
            end
        else
            f(dy, y_vec, ctx, t)
        end
    end

    nothing
end

"""
    set_dynamics_jacobian!(phase, var, fn)
    set_dynamics_jacobian!(fn, phase, var)   # do-block sugar

Register an analytic node-level Jacobian for one variable block of the dynamics.

Signature: `fn(dF, y, p, t)`
  - `dF :: Matrix{Float64}` — pre-zeroed, sized (ns × block_size) for this variable
  - `y, p, t`               — same arguments as `dynamics!(dy, y, p, t)`

Call once per variable block (state_var, control_var, param_var, ...).
Unregistered blocks fall back to AutoDiff (not yet implemented).
"""
set_dynamics_jacobian!(p::CollocationPhase,
                       var::DirectSolverVariable, fn::Function) =
    register_jacobian!(p, var, fn)

set_dynamics_jacobian!(fn::Function, p::CollocationPhase,
                       var::DirectSolverVariable) =
    set_dynamics_jacobian!(p, var, fn)   # do-block sugar

has_dynamics_jacobian(p::CollocationPhase, var::DirectSolverVariable) =
    has_registered_jacobian(p, var)

get_dynamics_jacobian(p::CollocationPhase, var::DirectSolverVariable) =
    registered_jacobian(p, var)

"""
    add_jacobian!(fn, phase, State())   do y, u, p, t, model ... end
    add_jacobian!(fn, phase, Control()) do y, u, p, t, model ... end

Register an analytic dynamics Jacobian using the domain-native 5-arg functional form.
`fn(y, u, p, t, model) -> Matrix` is wrapped into the internal mutating convention.
"""
function add_jacobian!(fn::Function, ph::CollocationPhase, ::State)
    ph_ref = ph
    wrapped! = function _jac_state_6arg!(dF, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph_ref), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            dF .= fn(y_s, u_s, nothing, t, ph_ref.model)
        else
            dF .= fn(y_vec, ctx.u, nothing, t, ph_ref.model)
        end
    end
    set_dynamics_jacobian!(ph, ph.state_var, wrapped!)
    nothing
end

function add_jacobian!(fn::Function, ph::CollocationPhase, ::Control)
    ph_ref = ph
    wrapped! = function _jac_ctrl_6arg!(dF, y_vec, ctx, t)
        reg = get(_phase_registry, objectid(ph_ref), nothing)
        if reg !== nothing && reg.state_type !== nothing
            Ty = promote_type(eltype(y_vec), eltype(ctx.u))
            y_s, u_s = _prepare_model!(reg, Ty, y_vec, ctx.u)
            dF .= fn(y_s, u_s, nothing, t, ph_ref.model)
        else
            dF .= fn(y_vec, ctx.u, nothing, t, ph_ref.model)
        end
    end
    set_dynamics_jacobian!(ph, ph.control_var, wrapped!)
    nothing
end

"""
    add_jacobian!(fn, phase, TimeTag())

Register an analytic ∂f/∂t for the dynamics.  Physics layer only.

Signature:  fn(dft::Vector, y, ctx, t) → nothing  (fills ns-Vector in-place)

The transcription assembles ∂defect/∂t0 and ∂defect/∂tf from this per-node
partial via the mesh chain rule.  Users never write that chain rule;
each transcription (LGL, H-S, ZOH) owns it exactly once.
"""
function add_jacobian!(fn::Function, ph::CollocationPhase, ::TimeTag)
    register_jacobian!(ph, TimeTag(), fn)
    nothing
end

has_dynamics_jac_t(ph::CollocationPhase) = has_registered_jacobian(ph, TimeTag())
get_dynamics_jac_t(ph::CollocationPhase) = registered_jacobian(ph, TimeTag())

"""
    add_jacobian!(bf::BoundaryFunction, var::DirectSolverVariable, fn::Function)
    add_jacobian!(fn, bf, var)   # do-block sugar

Register an analytic Jacobian closure for the (BoundaryFunction, variable) pair.

`fn` is a zero-arg closure (captures phase via enclosing scope) that returns
a `Matrix{Float64}` of shape `(n_bc_outputs, nlp_length(phase, var))`.
Only call this for variables where `sparsity_structure` returns true.
"""
add_jacobian!(bf::BoundaryFunction, var::DirectSolverVariable, fn::Function) =
    register_jacobian!(bf, var, fn)

add_jacobian!(fn::Function, bf::BoundaryFunction, var::DirectSolverVariable) =
    add_jacobian!(bf, var, fn)   # do-block sugar

has_jacobian(bf::BoundaryFunction, var::DirectSolverVariable) =
    has_registered_jacobian(bf, var)

# Calls what was registered. See the 🔴 on the storage layer.
get_jacobian(bf::BoundaryFunction, var::DirectSolverVariable) =
    registered_jacobian(bf, var)()

function set_state!(p::CollocationPhase, v)
    p.state_var   = v
    p._n_states   = length(v.lower_bounds)
    nothing
end

function set_control!(p::CollocationPhase, v)
    p.control_var  = v
    p._n_controls  = length(v.lower_bounds)
    nothing
end

function set_parameter!(p::CollocationPhase, v)
    p.param_var = v
    p._n_params = length(v.lower_bounds)
    nothing
end

"""
    get_param_value(phase) -> Vector{Float64}

A phase's static parameters. Before `solve!` these are the guess given to `Vary(parameter, ...)`;
after it, the values found.

# Returns
The parameter vector, in whatever units the dynamics read it as.

# Example
```julia
get_param_value(phase)
```
"""
get_param_value(p::CollocationPhase) = p._params

set_initial_time!(p::CollocationPhase, v) = (p.t0_var = v; nothing)
set_final_time!(p::CollocationPhase, v)   = (p.tf_var = v; nothing)

function add_constraint!(p::CollocationPhase, con)
    if con isa BoundaryConstraint
        phases = con.calc.phases
        if length(phases) != 1 || phases[1] !== p
            throw(ArgumentError(
                "Cross-phase BoundaryFunctions cannot be registered on a phase. " *
                "Use add_constraint!(seq, con) instead."
            ))
        end
    end
    push!(p.constraints, con)
    return nothing
end

add_path_constraint!(p::CollocationPhase, pc::PathConstraint) =
    (push!(p.path_constraints, pc); nothing)

# Legacy stub overload (keyword form) — kept for compatibility
add_path_constraint!(p::CollocationPhase, f; lower_bounds, upper_bounds) =
    (push!(p.path_constraints, (f=f, lb=lower_bounds, ub=upper_bounds)); nothing)

set_objective!(p::CollocationPhase, obj) = (p.objective = obj; nothing)

"""
    add_objective_jacobian!(obj, var, fn)
    add_objective_jacobian!(fn, obj, var)   # do-block sugar

Register an analytic gradient closure for one variable block of the Mayer objective.

`fn` is a zero-arg closure returning a `Vector{Float64}` of length
`nlp_length(phase, var)` — i.e. the partial ∂J/∂(block_j).
Only register variables for which the objective actually depends on that block.

SequenceManager calls `objective_gradient_chunk(phase, var)` for each variable
block and places the result at the correct global column offset — the same
pattern as `jacobian_chunk` for constraint blocks.

For :Max objectives, `objective_gradient_chunk` negates the partial automatically.
"""
add_objective_jacobian!(obj::MayerObjective,
                        var::DirectSolverVariable, fn::Function) =
    register_jacobian!(obj, var, fn)

add_objective_jacobian!(fn::Function, obj::MayerObjective,
                        var::DirectSolverVariable) =
    add_objective_jacobian!(obj, var, fn)   # do-block sugar

has_objective_jacobian(obj::MayerObjective, var::DirectSolverVariable) =
    has_registered_jacobian(obj, var)

# Calls what was registered. See the 🔴 on the storage layer.
get_objective_jacobian(obj::MayerObjective, var::DirectSolverVariable) =
    registered_jacobian(obj, var)()

# ── BolzaObjective uses the same jac_fns pattern ─────────────────────────────

add_objective_jacobian!(obj::BolzaObjective,
                        var::DirectSolverVariable, fn::Function) =
    register_jacobian!(obj, var, fn)

add_objective_jacobian!(fn::Function, obj::BolzaObjective,
                        var::DirectSolverVariable) =
    add_objective_jacobian!(obj, var, fn)

has_objective_jacobian(obj::BolzaObjective, var::DirectSolverVariable) =
    has_registered_jacobian(obj, var)

# Calls what was registered. See the 🔴 on the storage layer.
get_objective_jacobian(obj::BolzaObjective, var::DirectSolverVariable) =
    registered_jacobian(obj, var)()

# ─────────────────────────────────────────────────────────────────────────────
# Accessor functions — return stored values (populated by framework at eval time)
# ─────────────────────────────────────────────────────────────────────────────

"""
    get_initial_state(phase) -> state

The state at the start of the phase, as the phase's state type when the phase declares one and as a
plain vector otherwise. Before `solve!` it is the guess; after it, the solution.

Components come back in the units of the dynamics, in the order the state type declares its fields.

# Returns
The phase's state type holding the initial values, or the raw `Vector{Float64}` when the phase
declares no state type or the values have not been populated.

# Example
```julia
y0 = get_initial_state(phase)
```
"""
function get_initial_state(p::CollocationPhase)
    reg = get(_phase_registry, objectid(p), nothing)
    reg !== nothing && reg.state_type !== nothing && !isempty(p._y0) &&
        return reg.state_type{Float64}(p._y0...)
    return p._y0
end

"""
    get_final_state(phase) -> state

The state at the end of the phase, as the phase's state type when the phase declares one and as a
plain vector otherwise. Before `solve!` it is the guess; after it, the solution.

Components come back in the units of the dynamics, in the order the state type declares its fields,
so a field is read by name: `get_final_state(phase).m` is the delivered mass of a phase whose state
carries one.

# Returns
The phase's state type holding the final values, or the raw `Vector{Float64}` when the phase declares
no state type or the values have not been populated.

# Example
```julia
yf = get_final_state(phase)
```
"""
function get_final_state(p::CollocationPhase)
    reg = get(_phase_registry, objectid(p), nothing)
    reg !== nothing && reg.state_type !== nothing && !isempty(p._yf) &&
        return reg.state_type{Float64}(p._yf...)
    return p._yf
end
"""
    get_initial_time(phase::CollocationPhase) -> Float64

The phase's initial time, in the time units of its dynamics. Before `solve!` it is the start of
`tspan`, or the guess given to `Vary(initial_time, ...)`; after it, the solution.

# Returns
The start time as a `Float64`, in the time units of the dynamics.

# Example
```julia
get_initial_time(phase)
```
"""
get_initial_time(p::CollocationPhase)  = p._t0

"""
    get_final_time(phase::CollocationPhase) -> Float64

The phase's final time, in the time units of its dynamics. Before `solve!` it is the end of
`tspan`, or the guess given to `Vary(final_time, ...)`; after it, the solution.
# Returns
The end time as a `Float64`, in the time units of the dynamics.

# Example
```julia
get_final_time(phase)
```
"""
get_final_time(p::CollocationPhase)    = p._tf

# ─────────────────────────────────────────────────────────────────────────────
# get_control
#   get_control(ctx, t)        — EvalContext / SystemContext (SciML dynamics interface)
#   get_control(phase, k)      — node-index access into stored U matrix
#   get_control_matrix(phase)  — full n_controls × N control matrix
# ─────────────────────────────────────────────────────────────────────────────

get_control(p, t) = zeros(1)                                    # fallback / stub
get_control(p::CollocationPhase, k::Int) = p._U[:, k]           # node-index
get_control_matrix(p::CollocationPhase)  = p._U

# ─────────────────────────────────────────────────────────────────────────────
# Phase sizing
# ─────────────────────────────────────────────────────────────────────────────

"""Total unique LGL nodes for the phase (same for state and control grids)."""
n_control_nodes(p::CollocationPhase) = n_unique_nodes(p.transcription)

"""
    get_node_times(phase) → Vector{Float64}

Return the physical time at each collocation node.  This is the only
post-processing helper users need; it hides build_mesh and the τ → t
mapping, which are transcription internals.
# Returns
The times of the phase's mesh nodes, in the time units of the dynamics, spanning the phase from its
initial time to its final time. A Hermite-Simpson phase of `n` steps has `2n + 1` of them, since the
midpoints are nodes too.

# Example
```julia
get_node_times(phase)
```
"""
get_node_times(p::CollocationPhase) =
    node_times(build_mesh(p.transcription), p._t0, p._tf)

"""Length of the flat decision vector: vec(Y) ++ [vec(U)] ++ [params] ++ [t0, tf]."""
function n_decisions(p::CollocationPhase)
    N  = n_control_nodes(p)
    nc = p.control_var !== nothing ? p._n_controls : 0
    np = p.param_var   !== nothing ? p._n_params   : 0
    return p._n_states * N + nc * N + np + 2
end

"""
Total number of constraint function outputs:
  - Defect equations: n_states × (nodes per interval, shared nodes counted once each side)
  - Boundary constraints: sum of lower_bounds lengths over BoundaryConstraints
"""
function n_constraints(p::CollocationPhase)
    n_defect = n_defect_rows(p.transcription, p._n_states)
    n_bc = sum(length(c.lower_bounds) for c in p.constraints
               if c isa BoundaryConstraint; init=0)
    N = n_control_nodes(p)
    n_path = sum(pc.n_pc * N for pc in p.path_constraints
                 if pc isa PathConstraint; init=0)
    return n_defect + n_bc + n_path
end

# ─────────────────────────────────────────────────────────────────────────────
# Decision vector packing  x = [vec(Y); vec(U); t0; tf]
# Y is n_states × N, Julia col-major → x[1:ns*N] = [y_1; y_2; ...; y_N]
# ─────────────────────────────────────────────────────────────────────────────

function get_decision_vector(p::CollocationPhase)
    N  = n_control_nodes(p)
    ns = p._n_states
    nc = p.control_var !== nothing ? p._n_controls : 0
    np = p.param_var   !== nothing ? p._n_params   : 0
    x  = zeros(n_decisions(p))
    if size(p._Y) == (ns, N)
        x[1 : ns*N] .= vec(p._Y)
    end
    offset = ns*N
    if p.control_var !== nothing && size(p._U) == (nc, N)
        x[offset+1 : offset+nc*N] .= vec(p._U)
        offset += nc*N
    end
    if p.param_var !== nothing && length(p._params) == np
        x[offset+1 : offset+np] .= p._params
        offset += np
    end
    x[offset+1] = p._t0
    x[offset+2] = p._tf
    return x
end

function set_decision_vector!(p::CollocationPhase, x::Vector{Float64})
    N  = n_control_nodes(p)
    ns = p._n_states
    nc = p.control_var !== nothing ? p._n_controls : 0
    np = p.param_var   !== nothing ? p._n_params   : 0
    p._Y  = reshape(x[1 : ns*N], ns, N)
    offset = ns*N
    if p.control_var !== nothing
        p._U = reshape(x[offset+1 : offset+nc*N], nc, N)
        offset += nc*N
    else
        p._U = Matrix{Float64}(undef, 0, N)
    end
    if p.param_var !== nothing
        p._params = x[offset+1 : offset+np]
        offset += np
    end
    p._t0 = x[offset+1]
    p._tf = x[offset+2]
    p._y0 = p._Y[:, 1]
    p._yf = p._Y[:, end]
    nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Constraint bounds  — returns (lb, ub) aligned with get_functions output
# Defects are equality constraints (lb = ub = 0).
# ─────────────────────────────────────────────────────────────────────────────

function get_constraint_bounds(p::CollocationPhase)
    n_defect = n_defect_rows(p.transcription, p._n_states)
    lb = zeros(n_defect)
    ub = zeros(n_defect)
    for c in p.constraints
        if c isa BoundaryConstraint
            append!(lb, c.lower_bounds)
            append!(ub, c.upper_bounds)
        end
    end
    N = n_control_nodes(p)
    for pc in p.path_constraints
        if pc isa PathConstraint
            append!(lb, repeat(pc.lower_bounds, N))
            append!(ub, repeat(pc.upper_bounds, N))
        end
    end
    return lb, ub
end

# ─────────────────────────────────────────────────────────────────────────────
# PhaseFunction — descriptor returned by function_list(phase)
#
# Carries the source object (DefectBlock or BoundaryConstraint), row count,
# and bounds.  Equality constraints have lb == ub.
#
# SequenceManager iterates function_list to:
#   - assign global row offsets
#   - build NLP constraint bounds
#   - request Jacobian chunks per function block
# ─────────────────────────────────────────────────────────────────────────────

struct DefectBlock end   # sentinel: phase computes LGL defects internally

struct PhaseFunction
    source :: Any                 # DefectBlock() or BoundaryConstraint
    n_nlp  :: Int                 # residual row count
    lb     :: Vector{Float64}     # constraint lower bounds, length n_nlp
    ub     :: Vector{Float64}     # constraint upper bounds, length n_nlp
    name   :: String              # display label; "defects" for DefectBlock, bf.name for BoundaryConstraint
end

# ─────────────────────────────────────────────────────────────────────────────
# variable_list / function_list
#
# Primary interface for SequenceManager at setup time.
# Phase declares its NLP structure in the order it wants in the global vector.
# SequenceManager assigns offsets; phase never knows global indices.
#
# variable_list returns the four SolverVariables in NLP order.
# SequenceManager dispatches on typeof(v.var) to determine expansion:
#   AbstractStateArray   → ns × N scalars, repeat bounds N times
#   AbstractControlArray → nc × N scalars, repeat bounds N times
#   AbstractTime         → 1 scalar, bounds as-is
#
# nlp_length / nlp_bounds are two-argument: (phase, var) so the phase
# supplies N without storing it inside the variable descriptor.
# ─────────────────────────────────────────────────────────────────────────────

variable_list(p::CollocationPhase) = begin
    vars = Any[p.state_var]
    p.control_var !== nothing && push!(vars, p.control_var)
    p.param_var   !== nothing && push!(vars, p.param_var)
    push!(vars, p.t0_var)
    push!(vars, p.tf_var)
    vars
end

function function_list(p::CollocationPhase)
    n_defect = n_defect_rows(p.transcription, p._n_states)
    result = PhaseFunction[
        PhaseFunction(DefectBlock(), n_defect, zeros(n_defect), zeros(n_defect), "defects")
    ]
    for c in p.constraints
        if c isa BoundaryConstraint
            n = length(c.lower_bounds)
            push!(result, PhaseFunction(c, n, c.lower_bounds, c.upper_bounds, c.calc.name))
        end
    end
    N = n_control_nodes(p)
    for pc in p.path_constraints
        if pc isa PathConstraint
            n_total = pc.n_pc * N
            lb_tiled = repeat(pc.lower_bounds, N)
            ub_tiled = repeat(pc.upper_bounds, N)
            push!(result, PhaseFunction(PathConstraintBlock(pc), n_total,
                                        lb_tiled, ub_tiled, pc.name))
        end
    end
    return result
end

# nlp_length / nlp_bounds — dispatch on var.var type, phase supplies N
nlp_length(p::CollocationPhase) = sum(nlp_length(p, v) for v in variable_list(p))

function nlp_length(p::CollocationPhase, v::DirectSolverVariable)
    N = n_control_nodes(p)
    if v.var isa AbstractStateArray   || v.var isa State;         return p._n_states   * N end
    if v.var isa AbstractControlArray || v.var isa Control;       return p._n_controls * N end
    if v.var isa AbstractParameter;                               return p._n_params       end
    if v.var isa AbstractTime || v.var isa InitialTime || v.var isa FinalTime; return 1 end
    throw(ArgumentError(
        "collocation variable must be a state, control, parameter, initial-time or " *
        "final-time variable; got $(typeof(v.var))"))
end

function nlp_bounds(p::CollocationPhase, v::DirectSolverVariable)
    N = n_control_nodes(p)
    if v.var isa AbstractStateArray || v.var isa State
        return repeat(v.lower_bounds, N), repeat(v.upper_bounds, N)
    end
    if v.var isa AbstractControlArray || v.var isa Control
        return repeat(v.lower_bounds, N), repeat(v.upper_bounds, N)
    end
    if v.var isa AbstractParameter
        return v.lower_bounds, v.upper_bounds
    end
    if v.var isa AbstractTime || v.var isa InitialTime || v.var isa FinalTime
        return v.lower_bounds, v.upper_bounds
    end
    throw(ArgumentError(
        "collocation variable must be a state, control, parameter, initial-time or " *
        "final-time variable; got $(typeof(v.var))"))
end

# nlp_length / nlp_bounds for PhaseFunction — single-arg, no phase needed
nlp_length(pf::PhaseFunction)  = pf.n_nlp
nlp_bounds(pf::PhaseFunction)  = (pf.lb, pf.ub)

# ─────────────────────────────────────────────────────────────────────────────
# variable_ranges
#
# Returns a Vector{UnitRange{Int}} in the same order as variable_list(phase),
# giving each variable's column indices in the phase-local decision vector x.
#
# Layout matches get_decision_vector / set_decision_vector!:
#   state   → 1 : ns*N
#   control → ns*N+1 : ns*N+nc*N
#   t0      → ns*N+nc*N+1 : ns*N+nc*N+1
#   tf      → ns*N+nc*N+2 : ns*N+nc*N+2   (= n_decisions)
#
# SequenceManager uses phase-local ranges plus a global column offset to place
# each phase's Jacobian chunk into the full NLP matrix.
# ─────────────────────────────────────────────────────────────────────────────

function variable_ranges(p::CollocationPhase)
    ranges = UnitRange{Int}[]
    offset = 0
    for v in variable_list(p)
        len = nlp_length(p, v)
        push!(ranges, offset+1 : offset+len)
        offset += len
    end
    ranges
end

# ─────────────────────────────────────────────────────────────────────────────
# sparsity_structure
#
# Block-level boolean dependency: does function_i depend on variable_j?
#
# Two-argument form  (function-level):
#   sparsity_structure(phase, pf) → Vector{Bool}, length == length(variable_list)
#   true at position j means pf has a nonzero Jacobian block w.r.t. variable j.
#
# One-argument form  (phase-level):
#   sparsity_structure(phase) → Matrix{Bool}, size (n_functions, n_variables)
#   Assembled by calling the function-level form for every PhaseFunction.
#
# Derived from the evaluation context, not guessed from the variable's type.
#
# The framework builds the context each problem function is handed, so it knows
# what that function can read. A function depends on what it is given and on
# nothing else, which makes the set of context fields an upper bound that is
# exact at the block level and needs no declaration from the user.
#
#   DefectBlock         → Y, U, P, t0, tf
#   PathConstraintBlock → fn!(g, y, ctx, t), ctx = EvalContext(model, u, params)
#                         and t = t0 + tau * (tf - t0), so all five again
#   BoundaryConstraint  → BoundaryContext(y0, yf, t0, tf, params): all but control
#
# This replaces a guess that returned state only for a boundary constraint and
# state plus control for a path constraint. Both could drop a real dependence,
# and neither gave a user any way to say so: the `depends` keyword the old
# comment promised here was never built.
#
# The cost is a denser pattern where a function does not read everything it is
# handed. That is the right way round. An over-declared entry is computed and
# found to be zero; a dropped one never reaches the solver, and nothing
# downstream can tell a derivative that is zero from one never asked for.
# ─────────────────────────────────────────────────────────────────────────────

function sparsity_structure(p::CollocationPhase, pf::PhaseFunction)
    vlist = variable_list(p)
    n = length(vlist)
    if pf.source isa DefectBlock || pf.source isa PathConstraintBlock
        # Both are evaluated along the path, from the state and control at a
        # node, the params in the EvalContext, and a node time that is
        # t0 + tau * (tf - t0). That is every variable the phase has.
        return fill(true, n)
    elseif pf.source isa BoundaryConstraint
        # BoundaryContext(y0, yf, t0, tf, params) carries no control, and that
        # is the only block a boundary function cannot reach.
        return Bool[v !== p.control_var for v in vlist]
    end
    return fill(false, n)
end

function sparsity_structure(p::CollocationPhase)
    flist = function_list(p)
    vlist = variable_list(p)
    # n_functions × n_variables Bool matrix
    return Bool[sparsity_structure(p, f)[j]
                for f in flist, j in eachindex(vlist)]
end

# ─────────────────────────────────────────────────────────────────────────────
# print_sparsity  — human-readable block-level sparsity table for a phase
#
# Example output for brachistochrone:
#
#   Sparsity structure: CollocationPhase(:brachistochrone)
#   ──────────────────────────────────────────────────────────────
#              state  control     t0     tf
#   defects      ■      ■         ■      ■
#   start_bc     ■      ·         ·      ·
#   target_bc    ■      ·         ·      ·
# ─────────────────────────────────────────────────────────────────────────────

# Role labels for variable_list positions — derived from variable type
const _VARIABLE_ROLES = ("state", "control", "t0", "tf")   # kept for legacy reference

function _var_role(v::DirectSolverVariable)
    (v.var isa AbstractStateArray   || v.var isa State)   && return "state"
    (v.var isa AbstractControlArray || v.var isa Control) && return "control"
    v.var isa AbstractParameter && return "param"
    (v.var isa AbstractTime || v.var isa InitialTime || v.var isa FinalTime) &&
        return isempty(v.name) ? "time" : v.name
    return "var"
end

function print_sparsity(io::IO, p::CollocationPhase)
    S     = sparsity_structure(p)
    flist = function_list(p)
    vlist = variable_list(p)

    # Row label: use pf.name if set, else "fn_i" fallback
    row_labels = [isempty(pf.name) ? "fn_$i" : pf.name  for (i, pf) in enumerate(flist)]
    # Column label: use v.name if set, else derive from variable type
    col_labels = [isempty(v.name) ? _var_role(v) : v.name for v in vlist]

    row_w = max(maximum(length, row_labels), 3)
    col_w = [max(length(h), 7) for h in col_labels]

    sep = "─" ^ (row_w + 2 + sum(col_w) + length(col_w) * 2)
    println(io, "\nSparsity structure: $(p)")
    println(io, sep)

    # Header row
    print(io, " " ^ (row_w + 2))
    for (j, h) in enumerate(col_labels)
        print(io, lpad(h, col_w[j] + 2))
    end
    println(io)

    # Data rows
    for (i, rl) in enumerate(row_labels)
        print(io, rpad(rl, row_w + 2))
        for j in eachindex(col_labels)
            sym = S[i, j] ? "■" : "·"
            print(io, lpad(sym, col_w[j] + 2))
        end
        println(io)
    end
    println(io, sep)
end

print_sparsity(p::CollocationPhase) = print_sparsity(stdout, p)

# ─────────────────────────────────────────────────────────────────────────────
# objective_gradient_chunk
#
# Returns the partial ∂J/∂(block_j) as a Vector{Float64} of length
# nlp_length(phase, var) — one chunk per (objective, variable) pair.
#
# Mirrors jacobian_chunk exactly: caller (SequenceManager) is responsible for
# placing this at the correct global column offset in the NLP gradient.
#
# A block with no registered partial is differentiated with ForwardDiff. A
# Bolza objective sums its two halves, each registered or differentiated.
#
# Sign convention: negated for :Max (NLP solver always minimizes).
# ─────────────────────────────────────────────────────────────────────────────

function objective_gradient_chunk(p::CollocationPhase, var::DirectSolverVariable)
    obj = p.objective
    # Phase may have no objective in a multi-phase problem — return zeros
    obj isa Union{MayerObjective, BolzaObjective} || return zeros(nlp_length(p, var))
    chunk = obj isa MayerObjective ? _mayer_chunk(p, obj, var) :
            _lagrange_chunk(p, obj, var) .+ _bolza_mayer_chunk(p, var)
    return obj.sense === :Max ? .-chunk : chunk
end

# Each half of a Bolza objective takes its own partial when one was declared and is
# differentiated otherwise, so a declared running-cost partial does not stand in for the
# terminal cost's, or the reverse. The Mayer half is the phase's MayerObjective, a function of
# the boundary state; the Bolza table holds only the Lagrange half's partials.
_mayer_chunk(p, obj::MayerObjective, var) =
    has_objective_jacobian(obj, var) ? get_objective_jacobian(obj, var) :
                                       _mayer_gradient_ad_chunk(p, obj, var)

_lagrange_chunk(p, obj::BolzaObjective, var) =
    has_objective_jacobian(obj, var) ? get_objective_jacobian(obj, var) :
                                       _bolza_gradient_ad_chunk(p, obj, var)

function _bolza_mayer_chunk(p::CollocationPhase, var)
    m = p._mayer_internal
    m === nothing && return zeros(nlp_length(p, var))
    return _mayer_chunk(p, m, var)
end

# AD gradient of the Lagrange integral in a BolzaObjective w.r.t. one variable block.
# Returns a Vector{Float64} of length nlp_length(p, var).
# Differentiates lagrange_fn(y_k, u_k, t_k) at each node and assembles the
# block-diagonal result weighted by LGL quadrature weights.
# t0/tf gradient is not yet implemented (returns zeros) — for fixed-time problems
# this is exact; for free-time problems register an analytic Jacobian instead.
function _bolza_gradient_ad_chunk(p::CollocationPhase, obj::BolzaObjective,
                                   var::DirectSolverVariable)
    mesh    = build_mesh(p.transcription)
    N       = mesh.N
    qs      = quadrature_weights(mesh, p._t0, p._tf)
    t_nodes = node_times(mesh, p._t0, p._tf)
    if var.var isa AbstractStateArray
        ns = p._n_states
        g  = zeros(ns * N)
        for k in 1:N
            y_k = p._Y[:, k]; u_k = p._U[:, k]; t_k = t_nodes[k]
            gk  = ForwardDiff.gradient(yd -> obj.lagrange_fn(yd, u_k, t_k), y_k)
            g[(k-1)*ns+1 : k*ns] .= qs[k] .* gk
        end
        return g
    elseif var.var isa AbstractControlArray
        nc = p._n_controls
        g  = zeros(nc * N)
        for k in 1:N
            y_k = p._Y[:, k]; u_k = p._U[:, k]; t_k = t_nodes[k]
            gk  = ForwardDiff.gradient(ud -> obj.lagrange_fn(y_k, ud, t_k), u_k)
            g[(k-1)*nc+1 : k*nc] .= qs[k] .* gk
        end
        return g
    elseif var.var isa AbstractTime
        # Chain rule through both quadrature weights and node times.
        # q_k = Δτ_k * (tf-t0)/2 * w_k  →  ∂q_k/∂tf =  q_k/(tf-t0),  ∂q_k/∂t0 = -q_k/(tf-t0)
        # t_k = t0 + τ_k*(tf-t0)         →  ∂t_k/∂tf = τ_k,           ∂t_k/∂t0 = 1-τ_k
        dt     = p._tf - p._t0
        τ      = mesh.τ_global
        L_k    = [obj.lagrange_fn(p._Y[:,k], p._U[:,k], t_nodes[k]) for k in 1:N]
        dLdt_k = [ForwardDiff.derivative(
                      tt -> obj.lagrange_fn(p._Y[:,k], p._U[:,k], tt), t_nodes[k])
                  for k in 1:N]
        if objectid(var) == objectid(p.t0_var)
            g = sum(-qs[k]/dt * L_k[k] + qs[k] * dLdt_k[k] * (1.0 - τ[k]) for k in 1:N)
        else
            g = sum( qs[k]/dt * L_k[k] + qs[k] * dLdt_k[k] * τ[k]         for k in 1:N)
        end
        return [g]
    end
    return zeros(nlp_length(p, var))
end

# AD gradient of a MayerObjective w.r.t. one variable block.
# Returns a Vector{Float64} of length nlp_length(p, var).
# Only called when no analytic Jacobian has been registered.
function _mayer_gradient_ad_chunk(p::CollocationPhase, obj::MayerObjective,
                                   var::DirectSolverVariable)
    y0   = p._Y[:, 1]
    yf   = p._Y[:, end]
    t0   = p._t0
    tf   = p._tf
    prms = p._params
    ns   = p._n_states
    N    = size(p._Y, 2)
    if var.var isa AbstractStateArray
        # ∂J/∂y0 and ∂J/∂yf are each (ns,); embed into full (ns*N,) vector
        g_y0 = ForwardDiff.gradient(
            y0d -> obj.fn(BoundaryContext(_state_named(p, y0d), _state_named(p, yf),  t0, tf, prms)), y0)
        g_yf = ForwardDiff.gradient(
            yfd -> obj.fn(BoundaryContext(_state_named(p, y0),  _state_named(p, yfd), t0, tf, prms)), yf)
        g = zeros(ns * N)
        g[1:ns] .= g_y0
        g[(N-1)*ns+1 : N*ns] .= g_yf
        return g
    elseif var.var isa AbstractTime
        if objectid(var) == objectid(p.t0_var)
            return [ForwardDiff.derivative(
                t0d -> obj.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0d, tf,  prms)), t0)]
        else
            return [ForwardDiff.derivative(
                tfd -> obj.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0,  tfd, prms)), tf)]
        end
    elseif var.var isa AbstractParameter
        isempty(prms) && return zeros(length(var.lower_bounds))
        return ForwardDiff.gradient(
            pv -> obj.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0, tf, pv)), prms)
    end
    return zeros(nlp_length(p, var))
end

# Evaluate the scalar objective for a single phase.
# Sign convention: negated for :Max so the NLP always minimizes.
function get_objective(p::CollocationPhase)
    obj = p.objective
    if obj isa MayerObjective
        ctx = BoundaryContext(_state_named(p, p._Y[:, 1]), _state_named(p, p._Y[:, end]), p._t0, p._tf, p._params)
        raw = applicable(obj.fn, ctx) ? obj.fn(ctx) : obj.fn()
        val = Float64(raw)
        return obj.sense === :Max ? -val : val
    elseif obj isa BolzaObjective
        mayer_val    = Float64(obj.mayer_fn())
        mesh         = build_mesh(p.transcription)
        qs           = quadrature_weights(mesh, p._t0, p._tf)
        t_nodes      = node_times(mesh, p._t0, p._tf)
        lagrange_val = sum(qs[k] * obj.lagrange_fn(p._Y[:, k], p._U[:, k], t_nodes[k])
                          for k in 1:mesh.N)
        val = mayer_val + lagrange_val
        return obj.sense === :Max ? -val : val
    else
        throw(ArgumentError(
            "phase :$(p.name) must carry an objective before it can be solved; " *
            "declare one with Objective(quantity, phase; sense = Min())"))
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# jacobian_chunk
#
# Returns a dense Matrix{Float64} of shape (pf.n_nlp, nlp_length(phase, var))
# for one (PhaseFunction, variable) Jacobian block.  Only called for (i,j)
# pairs where sparsity_structure(phase, pf)[j] == true.
#
# Dispatch:
#   DefectBlock        → the transcription's own defect Jacobians (defect_jacobian_*)
#   BoundaryConstraint → user-registered analytic Jacobian via add_jacobian!
#                        (AutoDiff fallback — not yet implemented)
# ─────────────────────────────────────────────────────────────────────────────

function jacobian_chunk(p::CollocationPhase, pf::PhaseFunction,
                        var::DirectSolverVariable)
    if pf.source isa DefectBlock
        return _defect_jacobian_chunk(p, var)
    elseif pf.source isa BoundaryConstraint
        bf = pf.source.calc
        has_jacobian(bf, var) && return get_jacobian(bf, var)
        return _bc_jacobian_ad_chunk(p, bf, var)
    elseif pf.source isa PathConstraintBlock
        return _path_jacobian_chunk(p, pf.source.pc, var)
    end
    throw(ArgumentError(
        "a phase function's source must be a declared quantity or a registered partial; " *
        "got $(typeof(pf.source))"))
end

# AD Jacobian of a BoundaryConstraint w.r.t. one variable block.
# Returns a Matrix{Float64} of shape (n_bc, nlp_length(p, var)).
# Only called when no analytic Jacobian has been registered.
# Requires the BoundaryFunction.fn to use the new (ctx::BoundaryContext) API.
function _bc_jacobian_ad_chunk(p::CollocationPhase, bf::BoundaryFunction,
                               var::DirectSolverVariable)
    y0   = p._Y[:, 1]
    yf   = p._Y[:, end]
    t0   = p._t0
    tf   = p._tf
    prms = p._params
    ns   = p._n_states
    N    = size(p._Y, 2)
    # Evaluate once to get output size
    n_bc = length(bf.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0, tf, prms)))
    if var.var isa AbstractStateArray
        # Compact blocks: ∂bc/∂y0 (n_bc×ns) and ∂bc/∂yf (n_bc×ns)
        # Embed into full (n_bc × ns*N) layout: y0 at cols 1:ns, yf at cols end-ns+1:end
        J_y0 = ForwardDiff.jacobian(
            y0d -> bf.fn(BoundaryContext(_state_named(p, y0d), _state_named(p, yf),  t0, tf, prms)), y0)
        J_yf = ForwardDiff.jacobian(
            yfd -> bf.fn(BoundaryContext(_state_named(p, y0),  _state_named(p, yfd), t0, tf, prms)), yf)
        J = zeros(n_bc, ns * N)
        J[:, 1:ns] .= J_y0
        J[:, (N-1)*ns+1 : N*ns] .= J_yf
        return J
    elseif var.var isa AbstractTime
        if objectid(var) == objectid(p.t0_var)
            return ForwardDiff.jacobian(
                t0v -> bf.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0v[1], tf,     prms)), [t0])
        else
            return ForwardDiff.jacobian(
                tfv -> bf.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0,     tfv[1], prms)), [tf])
        end
    elseif var.var isa AbstractParameter
        isempty(prms) && return zeros(n_bc, length(var.lower_bounds))
        return ForwardDiff.jacobian(
            pv -> bf.fn(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0, tf, pv)), prms)
    end
    return zeros(n_bc, nlp_length(p, var))
end

"""Emit the nonzero blocks of a path constraint Jacobian instead of assembling it.

A path condition is evaluated at one node, so the chunk is block diagonal: `n_pc` by `n_cols` on the
diagonal, `N` times. Assembling it costs `(n_pc·N)` by `(n_cols·N)`, which is the mesh squared for
the same information, so this is the one chunk besides the defects that cannot exist at a large
problem size.

Follows `_path_jacobian_chunk` on which derivative to use: the registered analytic block where there
is one, automatic differentiation where the constraint has no analytic blocks at all, and nothing
where some are registered and this one is not — which is that function's rule that a partly analytic
constraint treats its unregistered blocks as zero. Emitting nothing is how a zero chunk is spelled
here, since the value vector starts at zero.

The block is a buffer reused between calls, so a sink reads it before returning.

# Returns
`nothing`. Everything is delivered through `sink`.
"""
function path_blocks(sink, p::CollocationPhase, pc::PathConstraint,
                     var::DirectSolverVariable)
    mesh   = build_mesh(p.transcription)
    N      = mesh.N
    n_pc   = pc.n_pc
    n_cols = var.var isa AbstractStateArray ? p._n_states : p._n_controls
    n_cols == 0 && return nothing

    analytic = has_path_jacobian(pc, var)
    analytic || isempty(registered_jacobians(pc)) || return nothing

    jac_fn! = analytic ? get_path_jacobian(pc, var) : nothing
    dg      = zeros(n_pc, n_cols)

    for k in 1:N
        y_k = p._Y[:, k]
        t_k = p._t0 + (p._tf - p._t0) * mesh.τ_global[k]
        ctx = EvalContext(p.model, p._U[:, k], p._params)
        if analytic
            fill!(dg, 0.0)
            jac_fn!(dg, y_k, ctx, t_k)
        elseif var.var isa AbstractStateArray
            dg .= ForwardDiff.jacobian(
                yd -> (g = similar(yd, n_pc); pc.fn(g, yd, ctx, t_k); g), y_k)
        else
            dg .= ForwardDiff.jacobian(
                ud -> (g = similar(ud, n_pc);
                       pc.fn(g, y_k, EvalContext(p.model, ud, p._params), t_k); g),
                p._U[:, k])
        end
        sink(dg, (k - 1) * n_pc, (k - 1) * n_cols)
    end
    return nothing
end

function _path_jacobian_chunk(p::CollocationPhase, pc::PathConstraint,
                               var::DirectSolverVariable)
    if !has_path_jacobian(pc, var)
        # If the user has registered ANY analytic Jacobians on this path constraint,
        # treat unregistered blocks as exactly zero (user is in analytic mode).
        # Only fall back to AD if no analytic Jacobians at all are registered.
        if !isempty(registered_jacobians(pc))
            mesh = build_mesh(p.transcription)
            N    = mesh.N
            if var.var isa AbstractStateArray
                return zeros(pc.n_pc * N, p._n_states * N)
            elseif var.var isa AbstractControlArray
                return zeros(pc.n_pc * N, p._n_controls * N)
            else
                return zeros(pc.n_pc * N, nlp_length(p, var))
            end
        end
        return _path_jacobian_ad_chunk(p, pc, var)
    end
    mesh = build_mesh(p.transcription)
    N    = mesh.N
    n_pc = pc.n_pc
    jac_fn! = get_path_jacobian(pc, var)
    if var.var isa AbstractStateArray
        n_cols = p._n_states
    elseif var.var isa AbstractControlArray
        n_cols = p._n_controls
    else
        throw(ArgumentError(
        "a path constraint's Jacobian is taken with respect to state or control only, " *
        "since a path condition is evaluated at a node; got $(typeof(var.var))"))
    end
    J  = zeros(n_pc * N, n_cols * N)
    dg = zeros(n_pc, n_cols)
    for k in 1:N
        y_k = p._Y[:, k]
        t_k = p._t0 + (p._tf - p._t0) * mesh.τ_global[k]
        ctx = EvalContext(p.model, p._U[:, k], p._params)
        fill!(dg, 0.0)
        jac_fn!(dg, y_k, ctx, t_k)
        rows = (k-1)*n_pc+1 : k*n_pc
        cols = (k-1)*n_cols+1 : k*n_cols
        J[rows, cols] .= dg
    end
    J
end

# AD fallback for path constraint Jacobians.  Called when no analytic
# set_path_jacobian! has been registered for this (pc, var) pair.
# Differentiates pc.fn node-by-node, assembling the block-diagonal result
# into a dense (n_pc*N) × (n_col*N) matrix.
function _path_jacobian_ad_chunk(p::CollocationPhase, pc::PathConstraint,
                                  var::DirectSolverVariable)
    mesh  = build_mesh(p.transcription)
    N     = mesh.N
    n_pc  = pc.n_pc
    if var.var isa AbstractStateArray
        n_cols = p._n_states
    elseif var.var isa AbstractControlArray
        n_cols = p._n_controls
    else
        return zeros(n_pc * N, nlp_length(p, var))
    end
    J = zeros(n_pc * N, n_cols * N)
    for k in 1:N
        y_k  = p._Y[:, k]
        u_k  = p._U[:, k]
        t_k  = p._t0 + (p._tf - p._t0) * mesh.τ_global[k]
        rows = (k-1)*n_pc  + 1 : k*n_pc
        cols = (k-1)*n_cols + 1 : k*n_cols
        if var.var isa AbstractStateArray
            ctx_k = EvalContext(p.model, u_k, p._params)
            J[rows, cols] = ForwardDiff.jacobian(
                yd -> (g = similar(yd, n_pc); pc.fn(g, yd, ctx_k, t_k); g),
                y_k)
        else   # AbstractControlArray
            J[rows, cols] = ForwardDiff.jacobian(
                ud -> (g = similar(ud, n_pc);
                       pc.fn(g, y_k, EvalContext(p.model, ud, p._params), t_k);
                       g),
                u_k)
        end
    end
    J
end

# Non-autonomous correction for defect/t0 and defect/tf Jacobians.
# The analytic defect_jacobian_t0/tf functions assume autonomous dynamics
# (∂F/∂t = 0).  For non-autonomous dynamics the missing term at each node k is:
#   t0: J[:,1] -= ∂F(y_k,u_k,t_k)/∂t · (1 - τ_k)
#   tf: J[:,1] -= ∂F(y_k,u_k,t_k)/∂t · τ_k
# ∂F/∂t is computed via ForwardDiff.derivative; automatically zero for
# autonomous dynamics so this is safe to always apply.
function _add_nonauto_time_correction!(J, p::CollocationPhase, mesh,
                                        which::Symbol)
    ns = p._n_states
    N  = mesh.N
    t0 = p._t0; tf = p._tf
    for k in 1:N
        y_k   = p._Y[:, k]
        ctx_k = EvalContext(p.model, p._U[:, k], p._params)
        τ_k   = mesh.τ_global[k]
        w_k   = which === :t0 ? (1.0 - τ_k) : τ_k
        t_k   = t0 + τ_k * (tf - t0)
        dy_dt = ForwardDiff.derivative(
            tt -> (dy = zeros(typeof(tt), ns); p.dynamics(dy, y_k, ctx_k, tt); dy), t_k)
        J[(k-1)*ns+1 : k*ns, 1] .-= w_k .* dy_dt
    end
end

"""
    has_jacobian_blocks(phase, pf, var) -> Bool

Whether this function and variable pair can deliver its Jacobian as sub-blocks.

A producer that can is used in preference to [`jacobian_chunk`](@ref), which returns the whole
`pf.n_nlp` by `nlp_length(phase, var)` matrix. For a defect block against the state that matrix is
the trajectory squared, so on a large problem it is the difference between a solve and an
out-of-memory error.

False is the safe answer and the default: the assembler falls back to the chunk, which is correct
and merely expensive. So a transcription gains this one producer at a time rather than all at once.

# Returns
`true` when [`jacobian_blocks`](@ref) is defined for the pair.

# Example
```julia
has_jacobian_blocks(phase, first(function_list(phase)), phase.state_var)
```
"""
has_jacobian_blocks(::CollocationPhase, ::PhaseFunction, ::Any) = false

function has_jacobian_blocks(p::CollocationPhase, pf::PhaseFunction, var::DirectSolverVariable)
    state_or_control = var.var isa AbstractStateArray || var.var isa AbstractControlArray
    state_or_control || return false
    # A path constraint is evaluated at one node, so its chunk is block diagonal whatever the
    # transcription. Only the defects need the stencil, and only Hermite-Simpson has one here.
    pf.source isa PathConstraintBlock && return true
    pf.source isa DefectBlock         && return has_defect_blocks(p.transcription)
    return false
end

"""
    jacobian_blocks(sink, phase, pf, var)

Deliver a Jacobian chunk as sub-blocks, without assembling the chunk.

`sink(block, row_offset, col_offset)` is called once per nonzero block, with offsets counted from
the top left of the chunk so the producer needs to know nothing about where the assembler places it.
The block is a buffer reused between calls, so a sink reads it before returning.

Defined for the Hermite-Simpson defects against state and control, which are the two chunks that
grow with the square of the mesh. The parameter and time chunks are one column or a few, so they
grow with the mesh rather than its square and keep the [`jacobian_chunk`](@ref) path.

# Returns
`nothing`. Everything is delivered through `sink`.
"""
function jacobian_blocks(sink, p::CollocationPhase, pf::PhaseFunction, var::DirectSolverVariable)
    pf.source isa PathConstraintBlock && return path_blocks(sink, p, pf.source.pc, var)
    mesh = build_mesh(p.transcription)
    ns   = p._n_states
    if var.var isa AbstractStateArray
        jac_fn! = has_dynamics_jacobian(p, var) ?
            get_dynamics_jacobian(p, var) :
            _make_dynamics_jac_y_ad(p.dynamics, ns)
        defect_blocks_Y(sink, jac_fn!, p._Y, p._U, p._params, p.model, mesh, p._t0, p._tf)
    elseif var.var isa AbstractControlArray
        jac_fn! = has_dynamics_jacobian(p, var) ?
            get_dynamics_jacobian(p, var) :
            _make_dynamics_jac_u_ad(p.dynamics, ns)
        defect_blocks_U(sink, jac_fn!, p._Y, p._U, p._params, p.model, mesh, p._t0, p._tf)
    else
        throw(ArgumentError(
            "jacobian_blocks is defined for the Hermite-Simpson defects against state or " *
            "control; got $(typeof(var.var)). has_jacobian_blocks says which pairs have it."))
    end
    return nothing
end

function _defect_jacobian_chunk(p::CollocationPhase, var::DirectSolverVariable)
    mesh = build_mesh(p.transcription)
    ns   = p._n_states
    if var.var isa AbstractStateArray
        jac_fn! = has_dynamics_jacobian(p, var) ?
            get_dynamics_jacobian(p, var) :
            _make_dynamics_jac_y_ad(p.dynamics, ns)
        return defect_jacobian_Y(jac_fn!, p._Y, p._U, p._params, p.model,
                                 mesh, p._t0, p._tf)
    elseif var.var isa AbstractControlArray
        jac_fn! = has_dynamics_jacobian(p, var) ?
            get_dynamics_jacobian(p, var) :
            _make_dynamics_jac_u_ad(p.dynamics, ns)
        return defect_jacobian_U(jac_fn!, p._Y, p._U, p._params, p.model,
                                 mesh, p._t0, p._tf)
    elseif var.var isa AbstractParameter
        jac_fn! = has_dynamics_jacobian(p, var) ?
            get_dynamics_jacobian(p, var) :
            _make_dynamics_jac_p_ad(p.dynamics, ns)
        return defect_jacobian_P(jac_fn!, p._Y, p._U, p._params,
                                 p.model, mesh, p._t0, p._tf)
    elseif var.var isa AbstractTime && objectid(var) == objectid(p.t0_var)
        J = defect_jacobian_t0(p._Y, p._U, p._params, p.model, p.dynamics,
                               mesh, p._t0, p._tf)
        _add_nonauto_time_correction!(J, p, mesh, :t0)
        return J
    elseif var.var isa AbstractTime && objectid(var) == objectid(p.tf_var)
        J = defect_jacobian_tf(p._Y, p._U, p._params, p.model, p.dynamics,
                               mesh, p._t0, p._tf)
        _add_nonauto_time_correction!(J, p, mesh, :tf)
        return J
    end
    throw(ArgumentError(
        "a defect Jacobian is taken with respect to state, control, parameter, initial " *
        "time or final time; got $(typeof(var.var))"))
end

# ─────────────────────────────────────────────────────────────────────────────
# Non-autonomous time correction for HS
#
# Adds ∂F/∂t contributions to the t0 or tf Jacobian columns already computed
# by defect_jacobian_t0 / defect_jacobian_tf.
#
# For each point gi with normalized time τ_gi:
#   ∂r_i/∂tf correction: += h · HS_B[i,p] · (∂f_gi/∂t) · τ_gi
#   ∂r_i/∂t0 correction: += h · HS_B[i,p] · (∂f_gi/∂t) · (1 - τ_gi)
# ─────────────────────────────────────────────────────────────────────────────

function _add_nonauto_time_correction!(J, p::CollocationPhase, mesh::HermiteSimpsonMesh,
                                        which::Symbol)
    ns     = p._n_states
    t0     = p._t0;  tf = p._tf
    h_phys = mesh.h * (tf - t0)

    for k in 1:mesh.n_steps
        L = 2k - 1;  M = 2k;  R = 2k + 1
        for (p_local, gi) in enumerate((L, M, R))
            τ_gi = mesh.τ_global[gi]
            t_gi = t0 + τ_gi * (tf - t0)
            ctx  = EvalContext(p.model, p._U[:, gi], p._params)
            # w is the τ-weight for the chosen endpoint
            w = which === :t0 ? (1.0 - τ_gi) : τ_gi
            dy_dt = ForwardDiff.derivative(
                tt -> (dy = zeros(typeof(tt), ns);
                       p.dynamics(dy, p._Y[:, gi], ctx, tt); dy),
                t_gi)
            for i in 1:2
                b = HS_B[i, p_local]
                b == 0.0 && continue
                d_rows = (2(k-1) + i - 1) * ns + 1 : (2(k-1) + i) * ns
                J[d_rows, 1] .-= h_phys * b * w .* dy_dt
            end
        end
    end
end
