# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The sequence-level user API, solve_trajectory!, and CollocationManager.
#

# ─────────────────────────────────────────────────────────────────────────────
# add_sequence! overload for CollocationPhase
# Stores phase in a module-level side-channel dict keyed by objectid(seq).
# Sequence itself is not modified.
# ─────────────────────────────────────────────────────────────────────────────

# The collocation code below asks a sequence what shooting phases it holds and,
# when it holds some, hands them to the shooting solver. Both live in
# shooting_sequence.jl, which this file precedes, so they are declared here and
# given methods there. A transcription in another package adds its methods the
# same way.
# `solve!` is the verb a user writes, and what it solves depends on what it is
# given: a sequence, a single phase, an estimation problem. Those live in
# different packages, so the generic has to be declared somewhere all of them
# can see, and adding a method is how each one joins.
function solve! end

function sequence_sf_phases end
function solve_trajectory_shooting! end
function initialize_shooting! end
function _invalidate_sm end

const _sequence_phases = Dict{UInt64, Vector{CollocationPhase}}()

function add_sequence!(seq::Sequence, phase::CollocationPhase)
    phases = get!(() -> CollocationPhase[], _sequence_phases, objectid(seq))
    push!(phases, phase)
    return nothing
end

sequence_phases(seq::Sequence) =
    get(_sequence_phases, objectid(seq), CollocationPhase[])

# ─────────────────────────────────────────────────────────────────────────────
# Sequence-level user API
#
# CollocationManager is an internal assembly object.  Users interact only
# with Sequence.  A CollocationManager is built lazily on the first call after
# initialize!(seq) and cached.  The cache is only rebuilt when initialize! is
# called again (e.g. after meshrefinement).  set_decision_vector! updates phase
# arrays in-place without invalidating the cached manager.
# ─────────────────────────────────────────────────────────────────────────────

const _cm_cache = Dict{UInt64, Any}()   # values are CollocationManager; Any avoids forward-ref

function _get_or_build_cm(seq::Sequence)
    get!(() -> CollocationManager(seq), _cm_cache, objectid(seq))
end

function _invalidate_cm(seq::Sequence)
    delete!(_cm_cache, objectid(seq))
end

# initialize! — unified dispatch for collocation and shooting sequences
"""
Attach any partials declared for a phase's dynamics.

Deferred to `initialize!` because `set_dynamics_jacobian!` needs the phase's
state and control variables, and those arrive with `Vary` — after the
constructor has run. Declaring a partial and never solving is therefore
harmless, and declaration order does not matter.
"""
function _attach_dynamics_partials!(p)
    f = get(_raw_dynamics, objectid(p), nothing)
    f === nothing && return nothing
    for (q, jac) in partials_of(f)
        tag = _path_tag(q)
        tag === nothing && continue
        # The framework hands over (dF, y, u, ctx, t) with no model; a partial
        # takes the same arguments as the function it differentiates.
        set_dynamics_jacobian!(p, tag,
            (dF, y, u, pp, t) -> jac(dF, y, u, pp, t, p.model))
    end
    return nothing
end

# Things that must happen once the Sequence exists but are written before it
# does — a Link's constraints, for one. Filled in by the layer above.
const _INIT_HOOKS = Vector{Any}()

function initialize!(seq::Sequence)
    for h in _INIT_HOOKS
        h(seq)
    end
    for p in sequence_phases(seq)
        p isa CollocationPhase && _attach_dynamics_partials!(p)
    end
    sf_phases = sequence_sf_phases(seq)
    if !isempty(sf_phases)
        initialize_shooting!(seq)
        return nothing
    end
    for p in sequence_phases(seq)
        initialize!(p)
    end
    _invalidate_cm(seq)        # force rebuild with fresh sizes
    _get_or_build_cm(seq)
    return nothing
end

"""Where each declared position sits in the value vector, as compressed columns.

`get_jacobian` cannot be the route to the value vector at any real problem size. It allocates an
`n_funs` by `n_vars` matrix, which is 80 GB for a hundred thousand variables and a hundred thousand
constraints, and it is allocated again on every call. This is what replaces it: the pattern is
indexed once, and the assembly writes each chunk entry straight to its slot.

Columns are compressed because the pattern is already sorted by column and then by row, so the
lookup is a binary search inside one column rather than a hash of a coordinate pair.
"""
struct JacobianIndex
    colptr::Vector{Int}      # length nx + 1; column c occupies colptr[c] : colptr[c+1] - 1
    rowval::Vector{Int}      # the declared rows, ascending within each column
    ng::Int
    nx::Int
end

function JacobianIndex(pattern::SNOW.SparsePattern, ng::Int, nx::Int)
    order  = sortperm(collect(zip(pattern.cols, pattern.rows)))
    rowval = [pattern.rows[k] for k in order]
    cols   = [pattern.cols[k] for k in order]
    colptr = Vector{Int}(undef, nx + 1)
    k = 1
    for c in 1:nx
        colptr[c] = k
        while k <= length(cols) && cols[k] == c
            k += 1
        end
    end
    colptr[nx + 1] = k
    return JacobianIndex(colptr, rowval, ng, nx)
end

"""The slot for `(r, c)`, or zero when the pattern does not declare it."""
@inline function _slot(ix::JacobianIndex, r::Int, c::Int)
    (1 <= c <= ix.nx) || return 0
    lo, hi = ix.colptr[c], ix.colptr[c + 1] - 1
    lo <= hi || return 0
    k = searchsortedfirst(view(ix.rowval, lo:hi), r) + lo - 1
    return (k <= hi && ix.rowval[k] == r) ? k : 0
end

get_decision_vector(seq::Sequence)          = get_decision_vector(_get_or_build_cm(seq))
get_functions(seq::Sequence)                = get_functions(_get_or_build_cm(seq))
get_constraint_bounds(seq::Sequence)        = get_constraint_bounds(_get_or_build_cm(seq))
get_variable_bounds(seq::Sequence)          = get_variable_bounds(_get_or_build_cm(seq))
get_jacobian(seq::Sequence)                 = get_jacobian(_get_or_build_cm(seq))
get_objective_gradient(seq::Sequence)       = get_objective_gradient(_get_or_build_cm(seq))
get_objective(seq::Sequence)                = get_objective(_get_or_build_cm(seq))
build_sparsity_pattern(seq::Sequence)       = build_sparsity_pattern(_get_or_build_cm(seq))

function get_jacobian_values!(vals::AbstractVector{Float64}, seq::Sequence,
                              ix::JacobianIndex)
    get_jacobian_values!(vals, _get_or_build_cm(seq), ix)
end

function set_decision_vector!(seq::Sequence, x::Vector{Float64})
    # Update phase arrays in-place.  DO NOT invalidate the cache — the
    # CollocationManager offsets are invariant to x; only _Y/_U/_t0/_tf change.
    set_decision_vector!(_get_or_build_cm(seq), x)
    return nothing
end

# evaluate!(seq, x) — single entry point for solver callbacks.
# Sets the decision vector once; all subsequent accessors read updated state.
function evaluate!(seq::Sequence, x::Vector{Float64})
    set_decision_vector!(seq, x)
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# solve_trajectory! — thin SNOW wrapper for collocation NLP
#
# Uses SNOW's UserDeriv mode: fun!(g, df, dg, x) fills constraints g,
# objective gradient df, and Jacobian matrix dg (ng×nx), and returns
# the scalar objective.  All derivatives are analytic.
# ─────────────────────────────────────────────────────────────────────────────

# Declaring the dense pattern again, for comparing an answer against the one the declared pattern
# gives. The arithmetic is identical either way; what changes is how much of it the solver is told
# about, so a disagreement between the two is a mistake in the pattern and not a tolerance.
const _DENSE_JACOBIAN = Ref(false)

"""
    dense_jacobian!(flag) -> Bool

Declare the constraint Jacobian dense on every subsequent solve, or stop doing so.

Off by default. Turning it on is how a solve is compared against itself: the same problem solved
with the declared pattern and with everything declared has to reach the same answer, and a
difference beyond the solver's own tolerance means the pattern omits a derivative. That check is
worth more than a containment test, which can only say the positions are present.

# Returns
The flag, as set.

# Example
```julia
dense_jacobian!(true)
```
"""
dense_jacobian!(flag::Bool) = (_DENSE_JACOBIAN[] = flag)

"""
    dense_jacobian() -> Bool

Whether solves are currently declaring the constraint Jacobian dense.

# Returns
`true` when [`dense_jacobian!`](@ref) has been turned on.

# Example
```julia
dense_jacobian()
```
"""
dense_jacobian() = _DENSE_JACOBIAN[]

"""Every position of an `ng` by `nx` Jacobian, in column-major order."""
function _dense_pattern(ng::Int, nx::Int)
    n = ng * nx
    rows = Vector{Int}(undef, n)
    cols = Vector{Int}(undef, n)
    k = 0
    for c in 1:nx, r in 1:ng
        k += 1
        rows[k] = r
        cols[k] = c
    end
    return SNOW.SparsePattern(rows, cols)
end

"""The constraint Jacobian sparsity pattern declared to the solver, as SNOW's `SparsePattern`.

Built from the transcription's stencil rather than from an evaluation. The rule an omitted entry
breaks is not a tolerance: the solver never asks for that derivative, nothing reports it, and the
run converges to the wrong answer. `test_correctness_jacobian_sparsity.jl` is the check, and it is
why `build_sparsity_pattern` below is unused — that function reads the dense Jacobian at one point
and keeps what is nonzero there, which on the brachistochrone drops 210 of 331 real entries.

Where refinement is cheap and the structure is certain, this declares the structure. Everywhere else
it over-declares, because an entry that is present and always zero costs one multiply and an absent
one costs the answer:

- A Hermite-Simpson defect for step k reads nodes 2k-1, 2k and 2k+1, so the defect rows are banded.
  That band is the whole saving. Two smaller zeros inside it are not claimed: the first defect's
  midpoint state block is the identity rather than a full block, and its midpoint control block is
  empty, and neither is worth the chance of indexing it wrong.
- A path constraint row reads one node.
- A boundary constraint reads the first and last nodes. Both are declared whichever end the
  constraint names, because the automatic-differentiation chunk writes both and a `Boundary`
  constraint can genuinely read each.
- Params and both time columns are declared for every row. The path chunks return zeros there today,
  but the residual does depend on t0 and tf through the node time, so declaring them costs two
  entries per row and survives that hook being wired up.
- Any transcription other than Hermite-Simpson keeps the dense declaration for its defect rows,
  since the stencil above is Hermite-Simpson's and nothing else here has one.
- Linkage rows are declared dense against both phases they join.
"""
function jacobian_pattern(seq::Sequence, ng::Int, nx::Int; dense::Bool = _DENSE_JACOBIAN[])
    dense && return _dense_pattern(ng, nx)
    cm      = _get_or_build_cm(seq)
    entries = Set{Tuple{Int, Int}}()
    function add!(rows, cols)
        for c in cols, r in rows
            push!(entries, (r, c))
        end
        return nothing
    end

    for (i, p) in enumerate(cm.phases)
        col_off = cm.phase_var_offsets[i]
        row_off = cm.phase_fun_offsets[i]
        ns, nc  = p._n_states, p._n_controls
        nnp     = p._n_params
        N       = n_control_nodes(p)
        ranges  = variable_ranges(p)

        # variable_list order is state, [control], [params], t0, tf, and control and params are
        # present only when declared, so the blocks are picked out rather than computed by offset.
        k_block   = 1
        state_r   = ranges[k_block]; k_block += 1
        control_r = p.control_var !== nothing ? ranges[k_block] : nothing
        p.control_var !== nothing && (k_block += 1)
        param_r   = p.param_var !== nothing ? ranges[k_block] : nothing
        p.param_var !== nothing && (k_block += 1)
        t_cols    = [col_off + only(ranges[k_block]), col_off + only(ranges[k_block + 1])]
        p_cols    = param_r === nothing ? Int[] : collect(col_off .+ param_r)

        y0 = col_off + first(state_r) - 1
        u0 = control_r === nothing ? 0 : col_off + first(control_r) - 1
        y_at(k) = (y0 + (k - 1) * ns + 1):(y0 + k * ns)
        u_at(k) = control_r === nothing ? (1:0) : ((u0 + (k - 1) * nc + 1):(u0 + k * nc))
        all_y   = collect(col_off .+ state_r)
        all_u   = control_r === nothing ? Int[] : collect(col_off .+ control_r)

        fun_row = 0
        for pf in function_list(p)
            base = row_off + fun_row
            if pf.source isa DefectBlock
                groups = ns > 0 ? defect_row_groups(p.transcription, ns) :
                                  Tuple{UnitRange{Int}, Vector{Int}}[]
                if !isempty(groups)
                    for (local_rows, nodes) in groups
                        rows = (base + first(local_rows)):(base + last(local_rows))
                        for g in nodes
                            add!(rows, y_at(g))
                            add!(rows, u_at(g))
                        end
                        add!(rows, p_cols)
                        add!(rows, t_cols)
                    end
                else
                    rows = (base + 1):(base + pf.n_nlp)
                    add!(rows, all_y); add!(rows, all_u)
                    add!(rows, p_cols); add!(rows, t_cols)
                end
            elseif pf.source isa BoundaryConstraint
                rows = (base + 1):(base + pf.n_nlp)
                add!(rows, y_at(1))
                add!(rows, y_at(N))
                add!(rows, p_cols); add!(rows, t_cols)
            elseif pf.source isa PathConstraintBlock
                n_pc = pf.source.pc.n_pc
                for k in 1:N
                    rows = (base + (k - 1) * n_pc + 1):(base + k * n_pc)
                    add!(rows, y_at(k)); add!(rows, u_at(k))
                    add!(rows, p_cols);  add!(rows, t_cols)
                end
            else
                rows = (base + 1):(base + pf.n_nlp)
                add!(rows, 1:nx)
            end
            fun_row += pf.n_nlp
        end
    end

    # A linkage joins two phases through a closure that can read any of either phase's variables.
    for (k, lc) in enumerate(cm.linkages)
        rows = (cm.linkage_fun_offsets[k] + 1):(cm.linkage_fun_offsets[k] + length(lc.lower_bounds))
        for ph in (lc.phase1, lc.phase2)
            j = findfirst(q -> q === ph, cm.phases)
            j === nothing && continue
            for r in variable_ranges(ph)
                add!(rows, cm.phase_var_offsets[j] .+ r)
            end
        end
    end

    # Column-major, so the value vector walks the dense Jacobian in memory order.
    sorted = sort!(collect(entries), by = e -> (e[2], e[1]))
    return SNOW.SparsePattern([e[1] for e in sorted], [e[2] for e in sorted])
end

# =============================================================================
# check_sparsity — an omitted derivative is worse than a slow one
# =============================================================================

"""
    check_sparsity(seq; n_samples = 24, seed = 20260922, verbose = true)

Compare the declared constraint Jacobian sparsity pattern against what the Jacobian actually holds.

The pattern may over-declare freely: an entry that is present and always zero costs one
multiplication. Omitting one costs the answer, and silently, because the solver never asks for that
derivative and nothing reports its absence. This is the check for that, and it is the reason
[`dense_jacobian!`](@ref) exists as its end-to-end counterpart.

The Jacobian is evaluated at several points inside the variable bounds rather than at one, because a
structurally nonzero entry can be numerically zero anywhere. The partial of a thrust cone with
respect to the thrust is exactly zero wherever the thrust is, which is where a minimum-fuel guess
starts, so a pattern read from a single evaluation drops it and never recovers it.

# Arguments
- `seq`: the sequence, already assembled. Nothing is solved and the decision vector is left as found.

# Keyword arguments
- `n_samples`: how many points beyond the current one to sample.
- `seed`: fixes the sampling, so a run is repeatable.
- `verbose`: prints the comparison.

# Returns
The number of omitted positions, which is zero when the pattern is sound. Warns when it is not.

# Example
<!-- doc-fragment -->
```julia
check_sparsity(Sequence(phase))
```
"""
function check_sparsity(seq::Sequence; n_samples::Int = 24, seed::Int = 20260922,
                        verbose::Bool = true)
    rng    = Random.MersenneTwister(seed)
    x0     = get_decision_vector(seq)
    lx, ux = get_variable_bounds(seq)
    found  = Set{Tuple{Int, Int}}()
    ng, nx = 0, 0

    for k in 0:n_samples
        # Sample inside the bounds. An out-of-bounds point is not wrong for a structural question,
        # but it can put a user function where it was never meant to be evaluated, and a NaN there
        # reads as a structural zero.
        x = k == 0 ? copy(x0) :
            clamp.(x0 .+ max.(abs.(x0), 1.0) .* (2 .* rand(rng, length(x0)) .- 1), lx, ux)
        evaluate!(seq, x)
        J = get_jacobian(seq)
        ng, nx = size(J)
        for c in axes(J, 2), r in axes(J, 1)
            J[r, c] != 0.0 && push!(found, (r, c))
        end
    end
    evaluate!(seq, x0)          # leave the sequence as it was found

    pattern  = jacobian_pattern(seq, ng, nx; dense = false)
    declared = Set(zip(pattern.rows, pattern.cols))
    omitted  = setdiff(found, declared)
    outside  = count(r -> !(1 <= r <= ng), pattern.rows) +
               count(c -> !(1 <= c <= nx), pattern.cols)

    if verbose
        @printf("  %-22s %s
", "constraint Jacobian", "$(ng) by $(nx)")
        @printf("  %-22s %d
", "dense", ng * nx)
        @printf("  %-22s %d  (%.2f%% of dense)
", "declared", length(declared),
                100 * length(declared) / (ng * nx))
        @printf("  %-22s %d  (over %d sampled points)
", "nonzero in practice",
                length(found), n_samples + 1)
        @printf("  %-22s %d
", "omitted", length(omitted))
        @printf("  %-22s %d
", "outside the problem", outside)
        println()
        println(isempty(omitted) && outside == 0 ?
                "  the declared pattern holds every entry the Jacobian produced" :
                "  the declared pattern is not sound")
    end

    isempty(omitted) || @warn "check_sparsity: the declared pattern omits " *
                              "$(length(omitted)) position(s) the Jacobian produced. The solver " *
                              "will not ask for those derivatives. Compare against " *
                              "dense_jacobian!(true)."
    outside == 0 || @warn "check_sparsity: the declared pattern names $outside position(s) " *
                          "outside the $(ng) by $(nx) Jacobian."
    return length(omitted)
end

function _solve_phases!(seq::Sequence, options::SNOW.Options)
    !isempty(sequence_sf_phases(seq)) && return solve_trajectory_shooting!(seq, options)
    x0       = get_decision_vector(seq)
    lx, ux   = get_variable_bounds(seq)
    lg, ug   = get_constraint_bounds(seq)
    ng       = length(lg)

    # SNOW defaults to DensePattern, which tells IPOPT every constraint depends on every variable.
    # Declaring the pattern changes what `fun!` is handed: `dg` arrives as a value vector in pattern
    # order rather than as an ng x nx matrix.
    #
    # Only the user-derivative path is declared. Under `derivatives = :fd` SNOW computes the
    # Jacobian itself through a different cache, and a pattern there is separate work.
    # The dense route is kept deliberately, and it is the one `dense_jacobian!` selects. It assembles
    # the whole matrix, so it cannot serve a large problem, but it is a second implementation of the
    # same answer — which is what makes the differential in
    # `test_correctness_jacobian_sparsity_differential.jl` worth running. If both routes shared the
    # fill, a mistake in the fill would move both answers the same way and the comparison would see
    # nothing.
    if options.derivatives isa SNOW.UserDeriv && !_DENSE_JACOBIAN[]
        pattern = jacobian_pattern(seq, ng, length(x0); dense = false)
        index   = JacobianIndex(pattern, ng, length(x0))
        options = SNOW.Options(sparsity    = pattern,
                               derivatives = options.derivatives,
                               solver      = options.solver)

        function fun_sparse!(g, df, dg, x)
            evaluate!(seq, x)
            g  .= get_functions(seq)
            df .= get_objective_gradient(seq)
            get_jacobian_values!(dg, seq, index)
            return get_objective(seq)
        end

        xopt, fopt, info = SNOW.minimize(fun_sparse!, x0, ng, lx, ux, lg, ug, options)
        evaluate!(seq, xopt)   # leave phase state at optimum
        return (variables = xopt, objective = fopt, info = info)
    end

    function fun!(g, df, dg, x)
        evaluate!(seq, x)
        g  .= get_functions(seq)
        df .= get_objective_gradient(seq)
        dg .= get_jacobian(seq)      # ng×nx dense matrix
        return get_objective(seq)
    end

    xopt, fopt, info = SNOW.minimize(fun!, x0, ng, lx, ux, lg, ug, options)
    evaluate!(seq, xopt)   # leave phase state at optimum
    return (variables = xopt, objective = fopt, info = info)
end

function _solve_phases!(seq::Sequence)
    !isempty(sequence_sf_phases(seq)) && return solve_trajectory_shooting!(seq)
    ip_options = Dict(
        "max_iter"         => 2000,
        "tol"              => 1e-6,
        "print_level"      => 5,
        "output_file"      => tempname() * "_ipopt.out",
        "file_print_level" => 5,
    )
    options = SNOW.Options(derivatives = SNOW.UserDeriv(),
                           solver      = SNOW.IPOPT(ip_options))
    _solve_phases!(seq, options)
end

# ─────────────────────────────────────────────────────────────────────────────
# CollocationManager
#
# Assembles the global NLP from one or more CollocationPhases.  Each phase
# owns its local variable layout; CollocationManager assigns global offsets
# and concatenates decision vectors, constraint residuals, Jacobians, and
# the objective gradient.
# ─────────────────────────────────────────────────────────────────────────────

struct CollocationManager
    phases             ::Vector{CollocationPhase}
    phase_var_offsets  ::Vector{Int}    # global column start (0-based) per phase
    phase_fun_offsets  ::Vector{Int}    # global row start  (0-based) per phase
    n_vars             ::Int            # total NLP variable count
    n_funs             ::Int            # total NLP constraint count
    linkages           ::Vector{LinkageConstraint}
    linkage_fun_offsets::Vector{Int}    # global row start (0-based) per linkage
end

function CollocationManager(seq::Sequence)
    phases = sequence_phases(seq)
    isempty(phases) && throw(ArgumentError("a Sequence must hold at least one CollocationPhase to transcribe. " *
                             "Did you call add_sequence!(seq, phase)?"))
    for p in phases
        nlp_length(p) > 0 || throw(ArgumentError("phase :$(p.name) must declare at least one variable. " *
                                   "Call initialize!(phase) before building CollocationManager."))
    end
    n = length(phases)
    var_offsets = Vector{Int}(undef, n)
    fun_offsets = Vector{Int}(undef, n)
    col = 0
    row = 0
    for (i, p) in enumerate(phases)
        var_offsets[i] = col
        fun_offsets[i] = row
        col += nlp_length(p)
        row += sum(pf.n_nlp for pf in function_list(p))
    end
    # Linkage constraints
    lcs = sequence_linkages(seq)
    nl  = length(lcs)
    lc_fun_offsets = Vector{Int}(undef, nl)
    for (k, lc) in enumerate(lcs)
        lc_fun_offsets[k] = row
        row += length(lc.lower_bounds)
    end
    CollocationManager(phases, var_offsets, fun_offsets, col, row, lcs, lc_fun_offsets)
end

# ── Global decision vector ────────────────────────────────────────────────────

function get_decision_vector(cm::CollocationManager)
    x = Vector{Float64}(undef, cm.n_vars)
    for (i, p) in enumerate(cm.phases)
        xp = get_decision_vector(p)
        x[cm.phase_var_offsets[i]+1 : cm.phase_var_offsets[i]+length(xp)] .= xp
    end
    x
end

function set_decision_vector!(cm::CollocationManager, x::Vector{Float64})
    length(x) == cm.n_vars || throw(ArgumentError(
        "Decision vector length $(length(x)) does not match n_vars=$(cm.n_vars)"))
    for (i, p) in enumerate(cm.phases)
        np = nlp_length(p)
        set_decision_vector!(p, x[cm.phase_var_offsets[i]+1 : cm.phase_var_offsets[i]+np])
    end
    nothing
end

# ── Global constraint residuals ───────────────────────────────────────────────

function get_functions(cm::CollocationManager)
    F = Vector{Float64}(undef, cm.n_funs)
    for (i, p) in enumerate(cm.phases)
        fp = get_functions(p)
        F[cm.phase_fun_offsets[i]+1 : cm.phase_fun_offsets[i]+length(fp)] .= fp
    end
    # Linkage constraints
    for (k, lc) in enumerate(cm.linkages)
        p1   = lc.phase1
        p2   = lc.phase2
        ctx1 = BoundaryContext(_state_named(p1, p1._Y[:,1]),   _state_named(p1, p1._Y[:,end]),
                               p1._t0, p1._tf, p1._params)
        ctx2 = BoundaryContext(_state_named(p2, p2._Y[:,1]),   _state_named(p2, p2._Y[:,end]),
                               p2._t0, p2._tf, p2._params)
        res  = lc.fn(ctx1, ctx2)
        n_lc = length(res)
        F[cm.linkage_fun_offsets[k]+1 : cm.linkage_fun_offsets[k]+n_lc] .= res
    end
    F
end

function get_constraint_bounds(cm::CollocationManager)
    lb = Float64[]
    ub = Float64[]
    for p in cm.phases
        lbp, ubp = get_constraint_bounds(p)
        append!(lb, lbp)
        append!(ub, ubp)
    end
    for lc in cm.linkages
        append!(lb, lc.lower_bounds)
        append!(ub, lc.upper_bounds)
    end
    lb, ub
end

function get_variable_bounds(cm::CollocationManager)
    lx = Float64[]
    ux = Float64[]
    for (i, p) in enumerate(cm.phases)
        for v in variable_list(p)
            lbv, ubv = nlp_bounds(p, v)
            append!(lx, lbv)
            append!(ux, ubv)
        end
    end
    lx, ux
end

# ── Global Jacobian ───────────────────────────────────────────────────────────

# AD Jacobian of a LinkageConstraint w.r.t. one variable in one phase.
# `which_phase` is :phase1 or :phase2; the other phase's context is held fixed
# at its current nominal values.  Returns Matrix{Float64} of shape
# (n_lc, nlp_length(p, var)).
function _linkage_jacobian_ad_chunk(lc::LinkageConstraint, which_phase::Symbol,
                                    var::DirectSolverVariable)
    if which_phase == :phase1
        p    = lc.phase1
        pfix = lc.phase2
        ctx_fix = BoundaryContext(
            _state_named(pfix, pfix._Y[:,1]), _state_named(pfix, pfix._Y[:,end]),
            pfix._t0, pfix._tf, pfix._params)
        wrapper = ctx -> lc.fn(ctx, ctx_fix)
    else   # :phase2
        p    = lc.phase2
        pfix = lc.phase1
        ctx_fix = BoundaryContext(
            _state_named(pfix, pfix._Y[:,1]), _state_named(pfix, pfix._Y[:,end]),
            pfix._t0, pfix._tf, pfix._params)
        wrapper = ctx -> lc.fn(ctx_fix, ctx)
    end
    y0   = p._Y[:, 1]
    yf   = p._Y[:, end]
    t0   = p._t0
    tf   = p._tf
    prms = p._params
    ns   = p._n_states
    N    = size(p._Y, 2)
    n_lc = length(lc.lower_bounds)
    if var.var isa AbstractStateArray
        J_y0 = ForwardDiff.jacobian(
            y0d -> wrapper(BoundaryContext(_state_named(p, y0d), _state_named(p, yf),  t0, tf, prms)), y0)
        J_yf = ForwardDiff.jacobian(
            yfd -> wrapper(BoundaryContext(_state_named(p, y0),  _state_named(p, yfd), t0, tf, prms)), yf)
        J = zeros(n_lc, ns * N)
        J[:, 1:ns] .= J_y0
        J[:, (N-1)*ns+1 : N*ns] .= J_yf
        return J
    elseif var.var isa AbstractTime
        if objectid(var) == objectid(p.t0_var)
            return ForwardDiff.jacobian(
                t0v -> wrapper(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0v[1], tf,     prms)), [t0])
        else
            return ForwardDiff.jacobian(
                tfv -> wrapper(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0,     tfv[1], prms)), [tf])
        end
    elseif var.var isa AbstractParameter
        isempty(prms) && return zeros(n_lc, length(var.lower_bounds))
        return ForwardDiff.jacobian(
            pv -> wrapper(BoundaryContext(_state_named(p, y0), _state_named(p, yf), t0, tf, pv)), prms)
    end
    return zeros(n_lc, nlp_length(p, var))
end

function get_jacobian(cm::CollocationManager)
    J = zeros(cm.n_funs, cm.n_vars)
    for (i, p) in enumerate(cm.phases)
        col_off = cm.phase_var_offsets[i]
        row_off = cm.phase_fun_offsets[i]
        vars   = variable_list(p)
        funs   = function_list(p)
        ranges = variable_ranges(p)
        fun_row = 0
        for pf in funs
            sparsity = sparsity_structure(p, pf)
            for (j, var) in enumerate(vars)
                sparsity[j] || continue
                chunk = jacobian_chunk(p, pf, var)
                grow = row_off + fun_row + 1 : row_off + fun_row + size(chunk, 1)
                gcol = col_off .+ ranges[j]
                J[grow, gcol] .= chunk
            end
            fun_row += pf.n_nlp
        end
    end
    # Linkage Jacobians
    for (k, lc) in enumerate(cm.linkages)
        row_off = cm.linkage_fun_offsets[k]
        n_lc    = length(lc.lower_bounds)
        grows   = row_off+1 : row_off+n_lc
        # Phase 1 contribution
        i1       = findfirst(p -> p === lc.phase1, cm.phases)
        col_off1 = cm.phase_var_offsets[i1]
        ranges1  = variable_ranges(lc.phase1)
        for (j, var) in enumerate(variable_list(lc.phase1))
            chunk = _linkage_jacobian_ad_chunk(lc, :phase1, var)
            J[grows, col_off1 .+ ranges1[j]] .+= chunk
        end
        # Phase 2 contribution
        i2       = findfirst(p -> p === lc.phase2, cm.phases)
        col_off2 = cm.phase_var_offsets[i2]
        ranges2  = variable_ranges(lc.phase2)
        for (j, var) in enumerate(variable_list(lc.phase2))
            chunk = _linkage_jacobian_ad_chunk(lc, :phase2, var)
            J[grows, col_off2 .+ ranges2[j]] .+= chunk
        end
    end
    J
end

# ── Global objective value ───────────────────────────────────────────────────

function get_objective(cm::CollocationManager)
    # Sum across phases (supports future multi-phase Bolza objectives).
    total = 0.0
    for p in cm.phases
        isnothing(p.objective) || (total += get_objective(p))
    end
    total
end

# ── Global objective gradient ─────────────────────────────────────────────────

function get_objective_gradient(cm::CollocationManager)
    g = zeros(cm.n_vars)
    for (i, p) in enumerate(cm.phases)
        col_off = cm.phase_var_offsets[i]
        ranges  = variable_ranges(p)
        for (j, var) in enumerate(variable_list(p))
            chunk = objective_gradient_chunk(p, var)
            g[col_off .+ ranges[j]] .+= chunk
        end
    end
    g
end

# ── Sparse Jacobian support ───────────────────────────────────────────────────
#
# build_sparsity_pattern  — call once at problem setup (uses current x).
#   Returns (rows, cols, vals) where rows/cols are the fixed Int vectors
#   expected by sparse NLP solvers, and vals is the initial value vector.
#
# get_jacobian_values!    — call every iteration; fills vals in-place
#   without re-allocating rows/cols.

function build_sparsity_pattern(cm::CollocationManager)
    J = get_jacobian(cm)        # dense; evaluated at current point
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for col in 1:size(J,2), row in 1:size(J,1)
        if J[row,col] != 0.0
            push!(rows, row)
            push!(cols, col)
            push!(vals, J[row,col])
        end
    end
    rows, cols, vals
end

"""Fill the declared positions of the constraint Jacobian, without assembling it.

Mirrors `get_jacobian` block for block, writing each chunk entry to its slot instead of into a
matrix. Positions the pattern declares and the assembly does not produce keep the zero they were
given, which is what makes over-declaring free.

A nonzero the pattern does not declare raises. That case is the one this whole mechanism exists to
prevent — the solver would otherwise never ask for that derivative and would converge somewhere else
without reporting anything — so it is worth a check on the path that would hide it.
"""
function get_jacobian_values!(vals::AbstractVector{Float64}, cm::CollocationManager,
                              ix::JacobianIndex)
    fill!(vals, 0.0)

    # Walk each chunk column against the rows declared in that column, rather than looking up every
    # entry. A chunk covers a contiguous range of global rows and columns, and the declared rows
    # inside a column are ascending, so the two are merged: one search per column to find where the
    # chunk starts, then a walk over the declared rows only.
    #
    # Looking each entry up instead costs a binary search per entry, which on a dense chunk is more
    # than the matrix it was meant to replace: it made every problem in the benchmark slower than
    # assembling the matrix had been.
    #
    # This walks declared positions, so a nonzero the pattern does not declare is skipped rather than
    # raised. `check_sparsity` is what catches that, offline and over sampled points, because a
    # per-entry check here is the cost this function exists to avoid.
    @inline function place_chunk!(chunk, r0, c0)
        nrows = size(chunk, 1)
        @inbounds for cj in axes(chunk, 2)
            c = c0 + cj
            (1 <= c <= ix.nx) || continue
            lo, hi = ix.colptr[c], ix.colptr[c + 1] - 1
            lo > hi && continue
            k = searchsortedfirst(view(ix.rowval, lo:hi), r0 + 1) + lo - 1
            while k <= hi && ix.rowval[k] <= r0 + nrows
                v = chunk[ix.rowval[k] - r0, cj]
                v == 0.0 || (vals[k] += v)
                k += 1
            end
        end
        return nothing
    end

    for (i, p) in enumerate(cm.phases)
        col_off = cm.phase_var_offsets[i]
        row_off = cm.phase_fun_offsets[i]
        vars    = variable_list(p)
        ranges  = variable_ranges(p)
        fun_row = 0
        for pf in function_list(p)
            sparsity = sparsity_structure(p, pf)
            for (j, var) in enumerate(vars)
                sparsity[j] || continue
                r0 = row_off + fun_row
                c0 = col_off + first(ranges[j]) - 1
                if has_jacobian_blocks(p, pf, var)
                    jacobian_blocks(p, pf, var) do blk, lr, lc
                        place_chunk!(blk, r0 + lr, c0 + lc)
                    end
                else
                    place_chunk!(jacobian_chunk(p, pf, var), r0, c0)
                end
            end
            fun_row += pf.n_nlp
        end
    end

    for (k, lc) in enumerate(cm.linkages)
        r0 = cm.linkage_fun_offsets[k]
        for (which, ph) in ((:phase1, lc.phase1), (:phase2, lc.phase2))
            i = findfirst(q -> q === ph, cm.phases)
            i === nothing && continue
            c_off  = cm.phase_var_offsets[i]
            ranges = variable_ranges(ph)
            for (j, var) in enumerate(variable_list(ph))
                chunk = _linkage_jacobian_ad_chunk(lc, which, var)
                place_chunk!(chunk, r0, c_off + first(ranges[j]) - 1)
            end
        end
    end
    return nothing
end

Base.show(io::IO, cm::CollocationManager) = begin
    lk_str = isempty(cm.linkages) ? "" : ", $(length(cm.linkages)) linkage(s)"
    print(io, "CollocationManager($(length(cm.phases)) phase(s), " *
              "n_vars=$(cm.n_vars), n_funs=$(cm.n_funs)$lk_str)")
end
# ─────────────────────────────────────────────────────────────────────────────
# Compact show methods — prevent recursive struct explosion in the REPL
# ─────────────────────────────────────────────────────────────────────────────

Base.show(io::IO, lc::LinkageConstraint) = begin
    name_str = isempty(lc.name) ? "" : ", \"$(lc.name)\""
    print(io, "LinkageConstraint(:$(lc.phase1.name) → :$(lc.phase2.name)" *
              ", n=$(length(lc.lower_bounds))$name_str)")
end

Base.show(io::IO, p::CollocationPhase) =
    print(io, "CollocationPhase(:$(p.name))")

Base.show(io::IO, v::DirectSolverVariable) =
    print(io, "SolverVariable($(typeof(v.var).name.name), lb=$(v.lower_bounds), ub=$(v.upper_bounds))")

Base.show(io::IO, bf::BoundaryFunction) = begin
    phases_str = join([":$(p.name)" for p in bf.phases], ", ")
    name_str   = isempty(bf.name) ? "" : ", \"$(bf.name)\""
    jac_str    = isempty(registered_jacobians(bf)) ? "" : ", $(length(registered_jacobians(bf))) jac(s)"
    print(io, "BoundaryFunction($phases_str$name_str$jac_str)")
end

Base.show(io::IO, bc::BoundaryConstraint) =
    print(io, "BoundaryConstraint($(bc.calc), lb=$(bc.lower_bounds), ub=$(bc.upper_bounds))")

Base.show(io::IO, obj::MayerObjective) = begin
    phases_str = join([":$(p.name)" for p in obj.phases], ", ")
    jac_str    = isempty(registered_jacobians(obj)) ? "" : ", $(length(registered_jacobians(obj))) jac(s)"
    print(io, "MayerObjective($phases_str, sense=:$(obj.sense)$jac_str)")
end

Base.show(io::IO, ::PathConstraintBlock) =
    print(io, "PathConstraintBlock()")

Base.show(io::IO, pc::PathConstraint) = begin
    jac_str = isempty(registered_jacobians(pc)) ? "" : ", $(length(registered_jacobians(pc))) jac(s)"
    print(io, "PathConstraint(\"$(pc.name)\", n=$(pc.n_pc)$jac_str)")
end

Base.show(io::IO, ::DefectBlock) =
    print(io, "DefectBlock()")

# A transcription's show method lives with the transcription.

Base.show(io::IO, pf::PhaseFunction) =
    print(io, "PhaseFunction(\"$(pf.name)\", n_nlp=$(pf.n_nlp))")

# ─────────────────────────────────────────────────────────────────────────────
# solve_trajectory! dispatches a Sequence to event, collocation, or shooting machinery.
#
# An event graph and a transcribed phase are different problems with different machinery, and a
# Sequence can hold either. They were two functions of the same name in two places for a while,
# which meant whichever loaded second silently won. There is one function now, and it asks the
# sequence what it holds.
#
# Both phase registries have to be asked before an empty sequence is treated as an event graph. A
# sequence of shooting phases looked empty when only the collocation registry was checked, went
# down the event graph path, and died in Ipopt on untyped bounds.
# ─────────────────────────────────────────────────────────────────────────────

"""
    solve_trajectory!(seq::Sequence; record_iterations = false)
    solve_trajectory!(seq::Sequence, options::SNOW.Options; record_iterations = false)

Solve the trajectory optimization problem the sequence describes.

A sequence holds either events or phases, and which one decides the machinery
used. A sequence of events runs the event graph, propagating and applying each
event in order. A sequence of phases assembles the transcribed nonlinear
program and hands it to the solver with analytic derivatives.

# Arguments
- `seq`: the sequence to solve.
- `options`: SNOW solver options. Defaults differ by problem: an event graph
  uses finite differences, a transcribed phase uses its analytic derivatives.

# Keyword Arguments
- `record_iterations`: keep each solver iteration in spacecraft history, for
  plotting how the trajectory converged. Event sequences only.

# Returns
A named tuple carrying `variables`, `objective` and `info`, where `info` is the
solver's convergence status.

# Examples
```julia
seq = Sequence()
add_sequence!(seq, toi_event, prop_event)
result = solve_trajectory!(seq)
```
"""
function solve_trajectory!(seq::Sequence, options::SNOW.Options;
                           record_iterations::Bool = false)
    if isempty(sequence_phases(seq)) && isempty(sequence_sf_phases(seq))
        return _solve_event_graph!(seq, options; record_iterations = record_iterations)
    end
    return _solve_phases!(seq, options)
end

function solve_trajectory!(seq::Sequence; record_iterations::Bool = false)
    if isempty(sequence_phases(seq)) && isempty(sequence_sf_phases(seq))
        return _solve_event_graph!(seq; record_iterations = record_iterations)
    end
    return _solve_phases!(seq)
end
