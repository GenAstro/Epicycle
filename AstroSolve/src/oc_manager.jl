# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# OCManager — unified optimal-control manager for mixed-transcription sequences
#
# Supports any mix of CollocationPhase + AbstractShootingPhase phases in a single
# NLP, with cross-type continuity via `add_continuity!(fn, seq, p1, p2; ...)`.
#
# Phase interface contract assumed (already implemented by CollocationPhase,
# every phase type that implements the element interface):
#
#   variable_list(p)             → Vector{DirectSolverVariable}
#   variable_ranges(p)           → Vector{UnitRange{Int}}     (local cols per var)
#   nlp_length(p[, v])           → Int
#   nlp_bounds(p, v)             → (lb::Vector, ub::Vector)
#   function_list(p)             → Vector{PhaseFunction}
#   sparsity_structure(p, pf)    → Vector{Bool}
#   jacobian_chunk(p, pf, v)     → Matrix{Float64}
#   get_decision_vector(p)       → Vector{Float64}
#   set_decision_vector!(p, x)   → nothing
#   get_functions(p)             → Vector{Float64}
#   get_constraint_bounds(p)     → (lb, ub)
#   get_objective(p)             → Float64
#   objective_gradient_chunk(p,v)→ Vector{Float64}
#
# Cross-phase residuals receive `UnifiedBoundaryCtx` objects with the
# transcription-agnostic fields `(y0::Vector, yf::Vector, t0, tf)` so user
# closures can mix any two phase types without referencing transcription-
# specific ctx fields.
#
# Cross-phase Jacobian rows are computed via global-NLP central differences,
# which is correct for any phase whose `set_decision_vector!` updates the
# cached state read by `_unified_boundary_context`.
# =============================================================================

# ── Unified boundary context (passed to add_continuity! closures) ────────────
struct UnifiedBoundaryCtx{T<:Real}
    y0 :: Vector{T}
    yf :: Vector{T}
    t0 :: T
    tf :: T
end

# Dispatch — each phase type maps its native cache fields to the unified ctx.
function _unified_boundary_context(p::CollocationPhase)
    UnifiedBoundaryCtx(Vector{Float64}(p._Y[:, 1]),
                       Vector{Float64}(p._Y[:, end]),
                       Float64(p._t0), Float64(p._tf))
end

function _unified_boundary_context(p::AbstractShootingPhase)
    c = _shoot_boundary_context(p)
    UnifiedBoundaryCtx(collect(Float64.(c.y0)), collect(Float64.(c.yf)),
                       Float64(c.t0), Float64(c.tf))
end

# Sims-Flanagan and MGAnDSMs name their endpoint states x0/xf, and both keep
# mass outside the state. A collocation phase keeps mass in it. So the unified
# context reports state-and-mass as one vector, which is what makes a link
# between a custom transcription and a general one a plain equality over every
# component rather than a special case.
function _unified_boundary_context(p::SimsFlanaganPhase)
    c = _sf_boundary_context(p)
    UnifiedBoundaryCtx(vcat(collect(Float64.(c.x0)), Float64(c.m0)),
                       vcat(collect(Float64.(c.xf)), Float64(c.mf)),
                       Float64(c.t0), Float64(c.tf))
end

function _unified_boundary_context(p::MGAnDSMsPhase)
    c = _mga_boundary_context(p)
    UnifiedBoundaryCtx(vcat(collect(Float64.(c.x0)), Float64(c.m0)),
                       vcat(collect(Float64.(c.xf)), Float64(c.mf)),
                       Float64(c.t0), Float64(c.tf))
end

# ── ContinuityLink — unified cross-phase constraint ──────────────────────────
mutable struct ContinuityLink
    phase1       :: Any
    phase2       :: Any
    fn           :: Function                  # (ctx1, ctx2) -> Vector
    lower_bounds :: Vector{Float64}
    upper_bounds :: Vector{Float64}
    name         :: String
end

Base.show(io::IO, lk::ContinuityLink) =
    print(io, "ContinuityLink(:$(_name_of(lk.phase1)) → :$(_name_of(lk.phase2))" *
              ", n=$(length(lk.lower_bounds))" *
              (isempty(lk.name) ? "" : ", \"$(lk.name)\"") * ")")

_name_of(p) = hasproperty(p, :name) ? string(p.name) : "phase"

# ── OCManager struct ─────────────────────────────────────────────────────────
struct OCManager
    phases            :: Vector{Any}
    ordered_vars      :: Vector{DirectSolverVariable}
    global_var_ranges :: Dict{UInt64, UnitRange{Int}}
    phase_fun_offsets :: Vector{Int}              # 1-based row start per phase
    links             :: Vector{ContinuityLink}
    link_fun_offsets  :: Vector{Int}              # 1-based row start per link
    n_vars            :: Int
    n_funs            :: Int
end

function OCManager(phases::Vector,
                   links::Vector{ContinuityLink} = ContinuityLink[])
    isempty(phases) && throw(ArgumentError(
        "OCManager: phases must hold at least one phase; got an empty collection"))

    # 1. Deduplicated variable list (objectid identity → shared junction vars)
    seen     = Set{UInt64}()
    ord_vars = DirectSolverVariable[]
    for p in phases, v in variable_list(p)
        oid = objectid(v)
        if oid ∉ seen
            push!(ord_vars, v)
            push!(seen, oid)
        end
    end

    # 2. Global column ranges
    gvr = Dict{UInt64, UnitRange{Int}}()
    col = 1
    for v in ord_vars
        n = 0
        for p in phases
            for vv in variable_list(p)
                if objectid(vv) == objectid(v)
                    n = nlp_length(p, v)
                    break
                end
            end
            n > 0 && break
        end
        gvr[objectid(v)] = col : col + n - 1
        col += n
    end
    n_vars = col - 1

    # 3. Per-phase function offsets
    fun_offs = Vector{Int}(undef, length(phases))
    row = 1
    for (i, p) in enumerate(phases)
        fun_offs[i] = row
        row += sum(pf.n_nlp for pf in function_list(p); init = 0)
    end

    # 4. Link function offsets
    lk_offs = Vector{Int}(undef, length(links))
    for (k, lk) in enumerate(links)
        lk_offs[k] = row
        row += length(lk.lower_bounds)
    end
    n_funs = row - 1

    OCManager(Vector{Any}(phases), ord_vars, gvr, fun_offs,
              links, lk_offs, n_vars, n_funs)
end

Base.show(io::IO, om::OCManager) = print(io,
    "OCManager($(length(om.phases)) phase(s), n_vars=$(om.n_vars), " *
    "n_funs=$(om.n_funs), $(length(om.links)) link(s))")

# ── Decision vector ──────────────────────────────────────────────────────────
function get_decision_vector(om::OCManager)
    x = Vector{Float64}(undef, om.n_vars)
    for p in om.phases
        xp     = get_decision_vector(p)
        ranges = variable_ranges(p)
        for (j, v) in enumerate(variable_list(p))
            gcol = om.global_var_ranges[objectid(v)]
            x[gcol] .= xp[ranges[j]]
        end
    end
    x
end

function set_decision_vector!(om::OCManager, x::Vector{Float64})
    length(x) == om.n_vars || throw(ArgumentError(
        "OCManager.set_decision_vector!: x must have one entry per NLP variable, so " *
        "length $(om.n_vars); got $(length(x))"))
    for p in om.phases
        xp     = Vector{Float64}(undef, nlp_length(p))
        ranges = variable_ranges(p)
        for (j, v) in enumerate(variable_list(p))
            gcol = om.global_var_ranges[objectid(v)]
            xp[ranges[j]] .= x[gcol]
        end
        set_decision_vector!(p, xp)
    end
    nothing
end

# ── Constraint residuals ─────────────────────────────────────────────────────
function get_functions(om::OCManager)
    F = Vector{Float64}(undef, om.n_funs)
    for (i, p) in enumerate(om.phases)
        fp   = get_functions(p)
        roff = om.phase_fun_offsets[i]
        F[roff : roff + length(fp) - 1] .= fp
    end
    for (k, lk) in enumerate(om.links)
        roff = om.link_fun_offsets[k]
        n_lk = length(lk.lower_bounds)
        c1   = _unified_boundary_context(lk.phase1)
        c2   = _unified_boundary_context(lk.phase2)
        F[roff : roff + n_lk - 1] .= lk.fn(c1, c2)
    end
    F
end

# ── Bounds ────────────────────────────────────────────────────────────────────
function get_constraint_bounds(om::OCManager)
    lb = Float64[];  ub = Float64[]
    for p in om.phases
        lbp, ubp = get_constraint_bounds(p)
        append!(lb, lbp);  append!(ub, ubp)
    end
    for lk in om.links
        append!(lb, lk.lower_bounds);  append!(ub, lk.upper_bounds)
    end
    lb, ub
end

function get_variable_bounds(om::OCManager)
    lx = Vector{Float64}(undef, om.n_vars)
    ux = Vector{Float64}(undef, om.n_vars)
    for v in om.ordered_vars
        gcol = om.global_var_ranges[objectid(v)]
        for p in om.phases
            plist = variable_list(p)
            idx   = findfirst(vv -> objectid(vv) == objectid(v), plist)
            if !isnothing(idx)
                lb, ub  = nlp_bounds(p, v)
                lx[gcol] .= lb
                ux[gcol] .= ub
                break
            end
        end
    end
    lx, ux
end

# ── Objective ────────────────────────────────────────────────────────────────
function get_objective(om::OCManager)
    total = 0.0
    for p in om.phases
        isnothing(p.objective) || (total += get_objective(p))
    end
    total
end

function get_objective_gradient(om::OCManager)
    g = zeros(om.n_vars)
    for p in om.phases
        isnothing(p.objective) && continue
        for v in variable_list(p)
            gcol  = om.global_var_ranges[objectid(v)]
            chunk = objective_gradient_chunk(p, v)
            g[gcol] .+= chunk
        end
    end
    g
end

# ── Jacobian ─────────────────────────────────────────────────────────────────
function get_jacobian(om::OCManager)
    J = zeros(om.n_funs, om.n_vars)

    # Per-phase blocks
    for (i, p) in enumerate(om.phases)
        roff   = om.phase_fun_offsets[i]
        flist  = function_list(p)
        vlist  = variable_list(p)
        ranges = variable_ranges(p)
        fun_row = 0
        for pf in flist
            sp = sparsity_structure(p, pf)
            for (j, var) in enumerate(vlist)
                sp[j] || continue
                chunk = jacobian_chunk(p, pf, var)
                grows = roff + fun_row : roff + fun_row + size(chunk, 1) - 1
                gcol  = om.global_var_ranges[objectid(var)]
                J[grows, gcol] .+= chunk
            end
            fun_row += pf.n_nlp
        end
    end

    # Cross-phase continuity blocks — global-NLP central differences.
    for (k, lk) in enumerate(om.links)
        roff = om.link_fun_offsets[k]
        n_lk = length(lk.lower_bounds)
        rows = roff : roff + n_lk - 1
        for phs in (lk.phase1, lk.phase2)
            for var in variable_list(phs)
                haskey(om.global_var_ranges, objectid(var)) || continue
                gcol  = om.global_var_ranges[objectid(var)]
                chunk = _oc_link_jac_fd_chunk(lk, om, var)
                J[rows, gcol] .+= chunk
            end
        end
    end
    J
end

# Avoid double-counting when phase1 and phase2 share variables (currently rare
# but supported for safety): collect each unique (link, var) only once.
function _oc_link_jac_fd_chunk(lk::ContinuityLink, om::OCManager,
                                var::DirectSolverVariable;
                                h::Float64 = 1e-7)
    n_lk  = length(lk.lower_bounds)
    gcol  = om.global_var_ranges[objectid(var)]
    n_var = length(gcol)
    J     = zeros(n_lk, n_var)

    x_saved = get_decision_vector(om)
    for (k, gj) in enumerate(gcol)
        xp = copy(x_saved);  xp[gj] += h
        set_decision_vector!(om, xp)
        rp = lk.fn(_unified_boundary_context(lk.phase1),
                   _unified_boundary_context(lk.phase2))
        xm = copy(x_saved);  xm[gj] -= h
        set_decision_vector!(om, xm)
        rm = lk.fn(_unified_boundary_context(lk.phase1),
                   _unified_boundary_context(lk.phase2))
        @inbounds for i in 1:n_lk
            J[i, k] = (rp[i] - rm[i]) / (2h)
        end
    end
    set_decision_vector!(om, x_saved)
    return J
end

# Sparse-Jacobian wrappers (mirror the other managers)
function build_sparsity_pattern(om::OCManager)
    J    = get_jacobian(om)
    rows = Int[];  cols = Int[];  vals = Float64[]
    for col in 1:size(J, 2), row in 1:size(J, 1)
        if J[row, col] != 0.0
            push!(rows, row);  push!(cols, col);  push!(vals, J[row, col])
        end
    end
    rows, cols, vals
end

function get_jacobian_values!(vals::Vector{Float64}, om::OCManager,
                              rows::Vector{Int}, cols::Vector{Int})
    J = get_jacobian(om)
    @inbounds for k in eachindex(vals)
        vals[k] = J[rows[k], cols[k]]
    end
    nothing
end

# =============================================================================
# Sequence-level API
# =============================================================================

const _sequence_oc_links = Dict{UInt64, Vector{ContinuityLink}}()
const _oc_cache          = Dict{UInt64, OCManager}()

sequence_oc_links(seq::Sequence) =
    get(_sequence_oc_links, objectid(seq), ContinuityLink[])

function _invalidate_oc(seq::Sequence)
    delete!(_oc_cache, objectid(seq))
end

function _get_or_build_oc(seq::Sequence)
    get!(_oc_cache, objectid(seq)) do
        # Collocation phases come first, then shooting phases.  Phase order
        # affects only block layout — it is independent of continuity wiring,
        # which references phases by identity.
        phases = Vector{Any}()
        append!(phases, sequence_phases(seq))
        append!(phases, sequence_sf_phases(seq))
        links  = sequence_oc_links(seq)
        OCManager(phases, links)
    end
end

"""
    add_continuity!(fn, seq, phase1, phase2; lower_bounds, upper_bounds, name="")

Register a cross-phase continuity (or coupling) constraint between any two
phases (collocation ↔ shooting, shoot ↔ shoot, coll ↔ coll).

The user closure `fn(ctx1, ctx2)` receives two [`UnifiedBoundaryCtx`](@ref)
objects with fields `y0::Vector, yf::Vector, t0, tf` — transcription-agnostic.

For the default state+time match across mixed types, use the no-closure form
[`add_continuity!(seq, phase1, phase2)`](@ref).
"""
function add_continuity!(fn::Function, seq::Sequence, phase1, phase2;
                          lower_bounds::AbstractVector{<:Real},
                          upper_bounds::AbstractVector{<:Real},
                          name::String = "")
    length(lower_bounds) == length(upper_bounds) ||
        throw(ArgumentError(
            "add_continuity!: lower_bounds and upper_bounds must have the same length; " *
            "got $(length(lower_bounds)) and $(length(upper_bounds))"))
    lk = ContinuityLink(phase1, phase2, fn,
                        Float64.(lower_bounds), Float64.(upper_bounds), name)
    push!(get!(() -> ContinuityLink[], _sequence_oc_links, objectid(seq)), lk)
    _invalidate_oc(seq)
    return lk
end

# Default state+time continuity — works for any pair as long as both phases
# expose y0/yf vectors of equal length.  Dispatched when at least one side is
# a shooting phase; pure collocation↔collocation continues to use the existing
# `add_continuity!(::Collocation, ::Collocation)` (which goes through
# `add_linkage!`/`CollocationManager`).
function add_continuity!(seq::Sequence, phase1::AbstractShootingPhase, phase2)
    _generic_continuity!(seq, phase1, phase2)
end
function add_continuity!(seq::Sequence, phase1::CollocationPhase,
                          phase2::AbstractShootingPhase)
    _generic_continuity!(seq, phase1, phase2)
end

# Collocation↔collocation: route to OCManager only when the sequence already
# contains a shooting phase (mixed sequence).  Pure-collocation sequences
# continue to use the legacy `add_linkage!`/CollocationManager pipeline.
function add_continuity!(seq::Sequence, phase1::CollocationPhase,
                          phase2::CollocationPhase)
    if isempty(sequence_sf_phases(seq))
        # Legacy pure-collocation path (replicates the body in OptControlStubs.jl).
        ns = phase1._n_states
        ns == phase2._n_states ||
            throw(ArgumentError(
                "add_continuity!: both phases must carry the same number of states to be " *
                "linked; got $ns and $(phase2._n_states)"))
        lb = zeros(Float64, ns + 1)
        ub = zeros(Float64, ns + 1)
        add_linkage!(seq, phase1, phase2;
                     lower_bounds = lb, upper_bounds = ub,
                     name = "continuity_$(phase1.name)_to_$(phase2.name)",
        ) do ctx1, ctx2
            yf1 = [getfield(ctx1.yf, i) for i in 1:fieldcount(typeof(ctx1.yf))]
            y02 = [getfield(ctx2.y0, i) for i in 1:fieldcount(typeof(ctx2.y0))]
            [yf1 .- y02; ctx1.tf - ctx2.t0]
        end
    else
        _generic_continuity!(seq, phase1, phase2)
    end
end

_state_dim(p::CollocationPhase)      = p._n_states

# Read the dimension off the unified context rather than off a field. A shooting phase keeps mass
# beside its state and the unified context appends it, so the number a link has to match is the
# context's length and nothing else — deriving it here is what stops the two from disagreeing.
#
# This replaced `p.n_states`, which no shooting phase has: SimsFlanaganPhase and MGAnDSMsPhase are
# the only two subtypes and neither declares that field, so every add_continuity! involving a
# shooting phase raised a FieldError naming internal fields, which §9.9 exists to prevent.
_state_dim(p::AbstractShootingPhase) = length(_unified_boundary_context(p).y0)

function _generic_continuity!(seq::Sequence, phase1, phase2)
    n1 = _state_dim(phase1)
    n2 = _state_dim(phase2)
    n1 == n2 || throw(ArgumentError(
        "add_continuity!: both phases must carry the same number of states to be linked; " *
        "got $n1 and $n2"))
    n = n1
    add_continuity!(seq, phase1, phase2;
        lower_bounds = zeros(n + 1), upper_bounds = zeros(n + 1),
        name = "continuity_$(_name_of(phase1))_to_$(_name_of(phase2))",
    ) do ctx1, ctx2
        [ctx1.yf .- ctx2.y0; ctx1.tf - ctx2.t0]
    end
end

"""
    initialize_oc!(seq::Sequence)

Initialize every phase in `seq` (collocation phases via `initialize!(phase)`,
shooting phases are already self-initialized by their `set_*!` calls) and
build a fresh [`OCManager`](@ref).  Mirror of `initialize!`/`initialize_shooting!`.
"""
function initialize_oc!(seq::Sequence)
    for p in sequence_phases(seq)
        initialize!(p)
    end
    _invalidate_oc(seq)
    _get_or_build_oc(seq)
    return nothing
end

# =============================================================================
# solve_trajectory_oc! — SNOW driver for the unified manager
# =============================================================================

"""
    solve_trajectory_oc!(seq[, options])

SNOW + IPOPT wrapper that drives the unified [`OCManager`](@ref).  Supports
any mix of collocation and shooting phases linked via `add_continuity!`.
"""
function solve_trajectory_oc!(seq::Sequence, options::SNOW.Options)
    om = _get_or_build_oc(seq)

    x0       = get_decision_vector(om)
    lx, ux   = get_variable_bounds(om)
    lg, ug   = get_constraint_bounds(om)
    ng       = length(lg)

    function fun!(g, df, dg, x)
        set_decision_vector!(om, x)
        g  .= get_functions(om)
        df .= get_objective_gradient(om)
        dg .= get_jacobian(om)
        return get_objective(om)
    end

    xopt, fopt, info = SNOW.minimize(fun!, x0, ng, lx, ux, lg, ug, options)
    set_decision_vector!(om, xopt)
    return (variables = xopt, objective = fopt, info = info, manager = om)
end

function solve_trajectory_oc!(seq::Sequence)
    ip_options = Dict(
        "max_iter"         => 2000,
        "tol"              => 1e-6,
        "print_level"      => 5,
        "output_file"      => tempname() * "_ipopt.out",
        "file_print_level" => 5,
    )
    options = SNOW.Options(derivatives = SNOW.UserDeriv(),
                           solver      = SNOW.IPOPT(ip_options))
    solve_trajectory_oc!(seq, options)
end
