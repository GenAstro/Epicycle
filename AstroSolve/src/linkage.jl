# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# LinkageConstraint and SequenceConstraint: the two ways two phases join.
#

# ─────────────────────────────────────────────────────────────────────────────
# Sequence-level cross-phase constraints
#
# BoundaryFunctions that reference more than one phase cannot be registered
# on any single phase — they belong on the Sequence.  SequenceManager then
# has full knowledge of all phase column offsets and can place these rows
# correctly in the global Jacobian.
#
# Storage: module-level dict keyed by objectid(seq) so Sequence itself needs
# no new fields (Sequence lives in AstroSolve which we don't modify here).
# ─────────────────────────────────────────────────────────────────────────────

const _sequence_constraints = Dict{UInt64, Vector{Any}}()

function add_constraint!(seq::Sequence, con)
    key = objectid(seq)
    if !haskey(_sequence_constraints, key)
        _sequence_constraints[key] = Any[]
    end
    push!(_sequence_constraints[key], con)
    return nothing
end

sequence_constraints(seq::Sequence) =
    get(_sequence_constraints, objectid(seq), Any[])

# ─────────────────────────────────────────────────────────────────────────────
# LinkageConstraint — connects two CollocationPhases at their boundaries.
#
# The user-supplied function `fn(ctx1, ctx2) -> AbstractVector` receives a
# BoundaryContext for each phase and returns a residual vector that must lie
# within [lower_bounds, upper_bounds] at the solution.
# AD Jacobians are computed automatically via ForwardDiff.
# ─────────────────────────────────────────────────────────────────────────────

struct LinkageConstraint
    phase1       ::CollocationPhase
    phase2       ::CollocationPhase
    fn           ::Function          # (ctx1::BoundaryContext, ctx2::BoundaryContext) -> Vector
    lower_bounds ::Vector{Float64}
    upper_bounds ::Vector{Float64}
    name         ::String
end

const _sequence_linkages = Dict{UInt64, Vector{LinkageConstraint}}()

"""
    add_linkage!(fn, seq, phase1, phase2; lower_bounds, upper_bounds, name="")

Register a linkage constraint between two phases.  Do-block form:

    add_linkage!(seq, phase1, phase2;
                 lower_bounds = [...], upper_bounds = [...], name = "...") do ctx1, ctx2
        [ ctx1.yf.r - ctx2.y0.r,
          ctx1.tf   - ctx2.t0   ]
    end

`fn(ctx1, ctx2)` receives a `BoundaryContext` for each phase and returns a
residual vector whose length must equal that of `lower_bounds`/`upper_bounds`.
Equality constraints use identical lower and upper bounds (typically zeros).
"""
function add_linkage!(fn::Function, seq::Sequence,
                      phase1::CollocationPhase, phase2::CollocationPhase;
                      lower_bounds::AbstractVector{<:Real},
                      upper_bounds::AbstractVector{<:Real},
                      name::String = "")
    length(lower_bounds) == length(upper_bounds) ||
        throw(ArgumentError(
            "add_linkage!: lower_bounds and upper_bounds must have the same length; " *
            "got $(length(lower_bounds)) and $(length(upper_bounds))"))
    lc = LinkageConstraint(phase1, phase2, fn,
                           Float64.(lower_bounds), Float64.(upper_bounds), name)
    key = objectid(seq)
    if !haskey(_sequence_linkages, key)
        _sequence_linkages[key] = LinkageConstraint[]
    end
    push!(_sequence_linkages[key], lc)
    return lc
end

sequence_linkages(seq::Sequence) =
    get(_sequence_linkages, objectid(seq), LinkageConstraint[])

# ─────────────────────────────────────────────────────────────────────────────
# SequenceConstraint — sequence-level constraint spanning two shooting phases.
#
# Similar to LinkageConstraint for collocation, but for AbstractShootingPhase.
# fn(ctx1::MGABoundaryContext, ctx2::MGABoundaryContext) -> Vector
#
# Jacobians keyed by objectid(var), matching the global_var_ranges layout in
# ShootingManager.  Analytic Jacobians registered via add_jacobian!(sc, var)
# or auto-differentiated via _seq_jac_ad_chunk if none registered.
# ─────────────────────────────────────────────────────────────────────────────

mutable struct SequenceConstraint
    phase1       ::Any    # AbstractShootingPhase
    phase2       ::Any    # AbstractShootingPhase
    fn           ::Function          # (ctx1, ctx2) -> Vector
    lower_bounds ::Vector{Float64}
    upper_bounds ::Vector{Float64}
    name         ::String
    jac_fns      ::Dict{UInt64, Function}   # keyed by objectid(var); zero-arg closures
end
