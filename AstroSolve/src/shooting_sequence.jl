# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# AbstractShootingPhase and ShootingManager: the shooting NLP.
#

# ═══════════════════════════════════════════════════════════════════════════════
# AbstractShootingPhase
#
# Protocol supertype for all shooting phases (Sims-Flanagan, etc).
# Defined here so OptControlStubs can dispatch on it in add_sequence! and
# ShootingManager before SimsFlanagan.jl is included.
# ═══════════════════════════════════════════════════════════════════════════════

# AbstractShootingPhase is declared in element_interface.jl: a
# sequence has to name the type whether or not a shooting package is loaded.

# ═══════════════════════════════════════════════════════════════════════════════
# ShootingManager
#
# Assembles the global NLP from one or more SimsFlanaganPhases, with full
# support for shared variables across phases (e.g. a junction time registered
# as tf of phase k and t0 of phase k+1).
#
# Variable deduplication:
#   Variables are identified by objectid(var).  Each unique variable occupies
#   exactly one contiguous block in the global decision vector, at the column
#   range given by global_var_ranges[objectid(var)].  Shared junction variables
#   contribute one column that receives Jacobian contributions from all phases
#   that reference it.
#
# Interface mirrors CollocationManager exactly:
#   get_decision_vector / set_decision_vector!
#   get_functions / get_constraint_bounds / get_variable_bounds
#   get_jacobian / build_sparsity_pattern / get_jacobian_values!
#   get_objective / get_objective_gradient
# ═════════════════════════════════════════════════════════════════════════════

struct ShootingManager
    phases            ::Vector                  # Vector{SimsFlanaganPhase}
    ordered_vars      ::Vector{DirectSolverVariable}   # unique, ordered by first encounter
    global_var_ranges ::Dict{UInt64, UnitRange{Int}}   # objectid(var) → global col range
    phase_fun_offsets ::Vector{Int}             # global row start (1-based) per phase
    seq_constraints   ::Vector{SequenceConstraint}     # sequence-level cross-phase constraints
    seq_con_offsets   ::Vector{Int}             # global row start (1-based) per seq constraint
    n_vars            ::Int
    n_funs            ::Int
end

function ShootingManager(phases::Vector,
                          seq_constraints::Vector{SequenceConstraint} = SequenceConstraint[])
    isempty(phases) && throw(ArgumentError(
        "ShootingManager: phases must hold at least one phase; got an empty collection"))

    # ── 1. Deduplicated global variable list ─────────────────────────────────
    seen      = Set{UInt64}()
    ord_vars  = DirectSolverVariable[]
    for p in phases
        for v in variable_list(p)
            oid = objectid(v)
            if oid ∉ seen
                push!(ord_vars, v)
                push!(seen, oid)
            end
        end
    end

    # ── 2. Global column ranges ───────────────────────────────────────────────
    gvr   = Dict{UInt64, UnitRange{Int}}()
    col   = 1
    for v in ord_vars
        n           = nlp_length(phases[1], v)   # all phases share same nlp_length dispatch
        # find any phase that contains this var for accurate nlp_length
        for p in phases
            if any(objectid(vv) == objectid(v) for vv in variable_list(p))
                n = nlp_length(p, v)
                break
            end
        end
        gvr[objectid(v)] = col : col+n-1
        col += n
    end
    n_vars = col - 1

    # ── 3. Function offsets ───────────────────────────────────────────────────
    fun_offs = Vector{Int}(undef, length(phases))
    row      = 1
    for (i, p) in enumerate(phases)
        fun_offs[i] = row
        row += n_constraints(p)
    end

    # ── 4. Sequence constraint offsets ────────────────────────────────────────
    sc_offs = Vector{Int}(undef, length(seq_constraints))
    for (k, sc) in enumerate(seq_constraints)
        sc_offs[k] = row
        row += length(sc.lower_bounds)
    end
    n_funs = row - 1

    ShootingManager(phases, ord_vars, gvr, fun_offs, seq_constraints, sc_offs, n_vars, n_funs)
end

# ── Decision vector ───────────────────────────────────────────────────────────

function get_decision_vector(sm::ShootingManager)
    x = Vector{Float64}(undef, sm.n_vars)
    for p in sm.phases
        xp     = get_decision_vector(p)
        ranges = variable_ranges(p)
        for (j, v) in enumerate(variable_list(p))
            gcol = sm.global_var_ranges[objectid(v)]
            x[gcol] .= xp[ranges[j]]
        end
    end
    x
end

function set_decision_vector!(sm::ShootingManager, x::Vector{Float64})
    length(x) == sm.n_vars || throw(ArgumentError(
        "ShootingManager.set_decision_vector!: x must have one entry per NLP variable, so " *
        "length $(sm.n_vars); got $(length(x))"))
    for p in sm.phases
        # Build a local decision vector from the global one
        xp     = Vector{Float64}(undef, nlp_length(p))
        ranges = variable_ranges(p)
        for (j, v) in enumerate(variable_list(p))
            gcol = sm.global_var_ranges[objectid(v)]
            xp[ranges[j]] .= x[gcol]
        end
        set_decision_vector!(p, xp)
    end
    nothing
end

# ── Constraint residuals ──────────────────────────────────────────────────────

# Generic boundary-context dispatch for sequence constraints.
# Each AbstractShootingPhase subtype defines a method returning a value-typed
# context (any struct).  The user closure `sc.fn(ctx1, ctx2)` receives whichever
# contexts the two phases produce — homogeneous (MGA↔MGA, ZOH↔ZOH) or mixed.
#
# Subtypes provide:  _shoot_boundary_context(p::MySubtype) -> MyCtx
function _shoot_boundary_context end

function get_functions(sm::ShootingManager)
    F = Vector{Float64}(undef, sm.n_funs)
    for (i, p) in enumerate(sm.phases)
        fp   = get_functions(p)
        roff = sm.phase_fun_offsets[i]
        F[roff : roff+length(fp)-1] .= fp
    end
    # Sequence constraints
    for (k, sc) in enumerate(sm.seq_constraints)
        roff  = sm.seq_con_offsets[k]
        n_con = length(sc.lower_bounds)
        ctx1  = _shoot_boundary_context(sc.phase1)
        ctx2  = _shoot_boundary_context(sc.phase2)
        F[roff : roff + n_con - 1] .= sc.fn(ctx1, ctx2)
    end
    F
end

function get_constraint_bounds(sm::ShootingManager)
    lb = Float64[];  ub = Float64[]
    for p in sm.phases
        lbp, ubp = get_constraint_bounds(p)
        append!(lb, lbp);  append!(ub, ubp)
    end
    for sc in sm.seq_constraints
        append!(lb, sc.lower_bounds)
        append!(ub, sc.upper_bounds)
    end
    lb, ub
end

function get_variable_bounds(sm::ShootingManager)
    lx = Vector{Float64}(undef, sm.n_vars)
    ux = Vector{Float64}(undef, sm.n_vars)
    for v in sm.ordered_vars
        gcol  = sm.global_var_ranges[objectid(v)]
        # find the phase that owns this var for its bounds
        for p in sm.phases
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

# ── Jacobian ──────────────────────────────────────────────────────────────────

function get_jacobian(sm::ShootingManager)
    J = zeros(sm.n_funs, sm.n_vars)
    for (i, p) in enumerate(sm.phases)
        roff  = sm.phase_fun_offsets[i]
        flist = function_list(p)
        vlist = variable_list(p)
        fun_row = 0
        for pf in flist
            sp = sparsity_structure(p, pf)
            for (j, var) in enumerate(vlist)
                sp[j] || continue
                chunk = jacobian_chunk(p, pf, var)
                grows = roff + fun_row : roff + fun_row + size(chunk,1) - 1
                gcol  = sm.global_var_ranges[objectid(var)]
                J[grows, gcol] .+= chunk   # .+= handles shared junction vars
            end
            fun_row += pf.n_nlp
        end
    end
    # Sequence constraint Jacobians
    for (k, sc) in enumerate(sm.seq_constraints)
        roff  = sm.seq_con_offsets[k]
        n_con = length(sc.lower_bounds)
        rows  = roff : roff + n_con - 1
        for phs in (sc.phase1, sc.phase2)
            for var in variable_list(phs)
                haskey(sm.global_var_ranges, objectid(var)) || continue
                gcol  = sm.global_var_ranges[objectid(var)]
                chunk = _seq_jacobian_chunk(sc, sm, phs, var)
                J[rows, gcol] .+= chunk
            end
        end
    end
    J
end

# ─────────────────────────────────────────────────────────────────────────────
# Sequence constraint Jacobian helpers
# ─────────────────────────────────────────────────────────────────────────────

function _seq_jacobian_chunk(sc::SequenceConstraint, sm::ShootingManager,
                              p_owner, var::DirectSolverVariable)
    # Analytic path
    if has_registered_jacobian(sc, var)
        chunk_phys = registered_jacobian(sc, var)()
        var_sc = length(var.scale) == 1 ? fill(var.scale[1], nlp_length_sm(sm, var)) :
                 repeat(var.scale, div(nlp_length_sm(sm, var), length(var.scale)))
        return chunk_phys .* transpose(var_sc)
    end
    # Phase-specific AD path — fast when available.
    if isdefined(@__MODULE__, :MGAnDSMsPhase) && p_owner isa MGAnDSMsPhase
        return _seq_jac_ad_chunk(sc, sm, p_owner, var)
    end
    # Generic finite-difference fallback — works for any AbstractShootingPhase
    # that implements `_shoot_boundary_context`.  Slow (O(n_var) sm-level
    # set/propagate cycles per chunk) but correct and transcription-agnostic.
    return _seq_jac_fd_chunk(sc, sm, p_owner, var)
end

function _seq_jac_fd_chunk(sc::SequenceConstraint, sm::ShootingManager,
                            p_owner, var::DirectSolverVariable;
                            h::Float64 = 1e-7)
    n_con = length(sc.lower_bounds)
    gcol  = sm.global_var_ranges[objectid(var)]
    n_var = length(gcol)
    J     = zeros(n_con, n_var)

    x_saved = get_decision_vector(sm)
    for (k, gj) in enumerate(gcol)
        xp = copy(x_saved);  xp[gj] += h
        set_decision_vector!(sm, xp)
        rp = sc.fn(_shoot_boundary_context(sc.phase1),
                   _shoot_boundary_context(sc.phase2))
        xm = copy(x_saved);  xm[gj] -= h
        set_decision_vector!(sm, xm)
        rm = sc.fn(_shoot_boundary_context(sc.phase1),
                   _shoot_boundary_context(sc.phase2))
        @inbounds for i in 1:n_con
            J[i, k] = (rp[i] - rm[i]) / (2h)
        end
    end
    set_decision_vector!(sm, x_saved)
    return J
end

function nlp_length_sm(sm::ShootingManager, var::DirectSolverVariable)
    length(sm.global_var_ranges[objectid(var)])
end

function _seq_jac_ad_chunk(sc::SequenceConstraint, sm::ShootingManager,
                            p_owner, var::DirectSolverVariable)
    ctx1_f64 = _mga_boundary_context(sc.phase1)
    ctx2_f64 = _mga_boundary_context(sc.phase2)
    is_phase1 = (p_owner === sc.phase1)
    sc_v = var.scale
    sh_v = var.shift
    phys_cur = copy(var.value)
    x_nlp    = (phys_cur .- sh_v) ./ sc_v
    n_nlp    = length(x_nlp)
    J = ForwardDiff.jacobian(xv -> begin
        T  = eltype(xv)
        phys_T = xv .* T.(sc_v) .+ T.(sh_v)
        if is_phase1
            ctx1 = _mga_context_replace(ctx1_f64, p_owner, var, phys_T)
            ctx2 = MGABoundaryContext{T}(ctx2_f64)
        else
            ctx1 = MGABoundaryContext{T}(ctx1_f64)
            ctx2 = _mga_context_replace(ctx2_f64, p_owner, var, phys_T)
        end
        sc.fn(ctx1, ctx2)
    end, x_nlp)
    J
end

function build_sparsity_pattern(sm::ShootingManager)
    J    = get_jacobian(sm)
    rows = Int[];  cols = Int[];  vals = Float64[]
    for col in 1:size(J,2), row in 1:size(J,1)
        if J[row,col] != 0.0
            push!(rows, row);  push!(cols, col);  push!(vals, J[row,col])
        end
    end
    rows, cols, vals
end

function get_jacobian_values!(vals::Vector{Float64}, sm::ShootingManager,
                               rows::Vector{Int}, cols::Vector{Int})
    J = get_jacobian(sm)
    @inbounds for k in eachindex(vals)
        vals[k] = J[rows[k], cols[k]]
    end
    nothing
end

# ── Objective ─────────────────────────────────────────────────────────────────

function get_objective(sm::ShootingManager)
    total = 0.0
    for p in sm.phases
        isnothing(p.objective) || (total += get_objective(p))
    end
    total
end

function get_objective_gradient(sm::ShootingManager)
    g = zeros(sm.n_vars)
    for p in sm.phases
        isnothing(p.objective) && continue
        for v in variable_list(p)
            gcol  = sm.global_var_ranges[objectid(v)]
            chunk = objective_gradient_chunk(p, v)
            g[gcol] .+= chunk
        end
    end
    g
end

Base.show(io::IO, sm::ShootingManager) =
    print(io, "ShootingManager($(length(sm.phases)) phase(s), " *
              "n_vars=$(sm.n_vars), n_funs=$(sm.n_funs))")

# ─────────────────────────────────────────────────────────────────────────────
# add_sequence! overload for SimsFlanaganPhase
# ─────────────────────────────────────────────────────────────────────────────

const _sequence_sf_phases = Dict{UInt64, Vector{Any}}()

function add_sequence!(seq::Sequence, phase::AbstractShootingPhase)
    phases = get!(() -> Any[], _sequence_sf_phases, objectid(seq))
    push!(phases, phase)
    return nothing
end

sequence_sf_phases(seq::Sequence) =
    get(_sequence_sf_phases, objectid(seq), Any[])

# ─────────────────────────────────────────────────────────────────────────────
# Sequence-level API for shooting phases
# ─────────────────────────────────────────────────────────────────────────────

const _sm_cache = Dict{UInt64, ShootingManager}()

function _get_or_build_sm(seq::Sequence)
    get!(() -> ShootingManager(sequence_sf_phases(seq),
                               sequence_mga_constraints(seq)),
         _sm_cache, objectid(seq))
end

function _invalidate_sm(seq::Sequence)
    delete!(_sm_cache, objectid(seq))
end

# initialize! for shooting sequences (no LGL mesh, but flushes cache)
function initialize_shooting!(seq::Sequence)
    _invalidate_sm(seq)
    _get_or_build_sm(seq)
    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# solve_trajectory! for shooting sequences
# ─────────────────────────────────────────────────────────────────────────────

function solve_trajectory_shooting!(seq::Sequence, options::SNOW.Options)
    sm = _get_or_build_sm(seq)

    x0       = get_decision_vector(sm)
    lx, ux   = get_variable_bounds(sm)
    lg, ug   = get_constraint_bounds(sm)
    ng       = length(lg)

    function fun!(g, df, dg, x)
        set_decision_vector!(sm, x)
        g  .= get_functions(sm)
        df .= get_objective_gradient(sm)
        dg .= get_jacobian(sm)      # ng×nx dense matrix
        return get_objective(sm)
    end

    xopt, fopt, info = SNOW.minimize(fun!, x0, ng, lx, ux, lg, ug, options)
    set_decision_vector!(sm, xopt)
    return (variables = xopt, objective = fopt, info = info,
            manager = sm)
end

function solve_trajectory_shooting!(seq::Sequence)
    ip_options = Dict(
        "max_iter"         => 2000,
        "tol"              => 1e-6,
        "print_level"      => 5,
        "output_file"      => tempname() * "_ipopt.out",
        "file_print_level" => 5,
    )
    options = SNOW.Options(derivatives = SNOW.UserDeriv(),
                           solver      = SNOW.IPOPT(ip_options))
    solve_trajectory_shooting!(seq, options)
end

# =============================================================================
# NEW PROPOSED API LAYER
# =============================================================================
#
# Tag structs, typed-closure registration, and solve! shim.
# Sits on top of the existing internals — old code unchanged.
#
# Design rules:
#   • Tag structs are the universal identity token — no string/symbol keys
#   • Closures receive typed structs directly — no wrapper context objects
#   • Jacobian tag selects which arg arrives in the closure
#   • Omitting add_jacobian! for a tag = zero Jacobian (framework assumes sparse zero)
#   • solve! is an alias for solve_trajectory!
