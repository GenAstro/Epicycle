# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The do-block layer for shooting phases: the same verbs the collocation half
# declares, dispatching on AbstractShootingPhase.
#

# =============================================================================
# NEW API LAYER — SHOOTING PHASES
#
# Mirrors the collocation do-block API so user-defined variables, constraints,
# and objectives are written the same way regardless of transcription type.
# Transcription-intrinsic variables (set_departure_vinf!, set_alpha!, etc.)
# and the boundary context fields (ctx.vinf_dep, ctx.mf, …) remain
# transcription-specific by design — they reflect what the method manages.
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# Handle types
# ─────────────────────────────────────────────────────────────────────────────

struct ShootingBoundaryConstraintHandle
    bc    ::BoundaryConstraint
    phase ::AbstractShootingPhase
end

mutable struct ShootingObjectiveHandle
    phase           ::AbstractShootingPhase
    sense           ::Union{Min,Max}
    _mayer_ref      ::Ref{Any}   # MayerObjective once add_mayer! is called
end
ShootingObjectiveHandle(phase, sense) =
    ShootingObjectiveHandle(phase, sense, Ref{Any}(nothing))

struct ShootingMayerHandle
    phase    ::AbstractShootingPhase
    internal ::MayerObjective
end

# ─────────────────────────────────────────────────────────────────────────────
# add_boundary_constraint! — shooting variant
#
#   bc = add_boundary_constraint!(phase; name, lower_bounds, upper_bounds,
#                                        scale, shift) do ctx
#       [...]
#   end
#
# fn(ctx) -> Vector  where ctx is the transcription's boundary context
# (SFBoundaryContext for Sims-Flanagan, MGABoundaryContext for MGAnDSMs, etc.)
# Returns a ShootingBoundaryConstraintHandle for add_jacobian!.
# ─────────────────────────────────────────────────────────────────────────────

function add_boundary_constraint!(fn::Function, phase::AbstractShootingPhase;
                                  name::Symbol,
                                  equality     = nothing,
                                  lower_bounds = nothing,
                                  upper_bounds = nothing,
                                  scale        = nothing,
                                  shift        = nothing)
    if equality !== nothing
        lb = Float64.(equality isa Number ? [equality] : equality)
        ub = copy(lb)
    else
        lb = Float64.(lower_bounds isa Number ? [lower_bounds] : lower_bounds)
        ub = Float64.(upper_bounds isa Number ? [upper_bounds] : upper_bounds)
    end
    n  = length(lb)
    sc = scale === nothing ? ones(n)  : Float64.(scale)
    sh = shift === nothing ? zeros(n) : Float64.(shift)

    bf = BoundaryFunction(phase; name = string(name)) do ctx
        fn(ctx)
    end
    bc = BoundaryConstraint(bf, lb, ub, sc, sh)
    add_constraint!(phase, bc)
    return ShootingBoundaryConstraintHandle(bc, phase)
end

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! — for ShootingBoundaryConstraintHandle
#
#   add_jacobian!(bc_handle, var) do val
#       [n_bc × length(val)] matrix     # val is copy of var.value
#   end
#
# The framework wraps fn into a zero-arg closure that reads var.value at
# evaluation time, matching the zero-arg convention used in jacobian_chunk.
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::ShootingBoundaryConstraintHandle,
                       var::DirectSolverVariable, fn::Function)
    zero_arg = () -> fn(copy(var.value))
    add_jacobian!(h.bc.calc, var, zero_arg)
    return nothing
end

add_jacobian!(fn::Function, h::ShootingBoundaryConstraintHandle,
              var::DirectSolverVariable) =
    add_jacobian!(h, var, fn)

# ─────────────────────────────────────────────────────────────────────────────
# set_objective! — shooting variant
#   Returns a ShootingObjectiveHandle; sense is Min() or Max().
# ─────────────────────────────────────────────────────────────────────────────

function set_objective!(phase::AbstractShootingPhase, sense::Union{Min,Max})
    return ShootingObjectiveHandle(phase, sense)
end

# ─────────────────────────────────────────────────────────────────────────────
# add_mayer! — shooting variant
#
#   myr = add_mayer!(obj) do ctx
#       scalar_value
#   end
#
# fn(ctx) -> scalar  where ctx is the transcription's boundary context.
# Installs a MayerObjective on the phase and returns a ShootingMayerHandle
# for add_jacobian!.
# ─────────────────────────────────────────────────────────────────────────────

function add_mayer!(fn::Function, obj::ShootingObjectiveHandle)
    ph        = obj.phase
    sense_sym = obj.sense isa Min ? :Min : :Max

    internal_obj = MayerObjective(ph; sense = sense_sym) do ctx
        fn(ctx)
    end
    set_objective!(ph, internal_obj)
    obj._mayer_ref[] = internal_obj
    return ShootingMayerHandle(ph, internal_obj)
end

add_mayer!(obj::ShootingObjectiveHandle, fn::Function) = add_mayer!(fn, obj)

# ─────────────────────────────────────────────────────────────────────────────
# add_jacobian! — for ShootingMayerHandle
#
#   add_jacobian!(myr, var) do val
#       Vector of length nlp_length(var)   # val is copy of var.value
#   end
# ─────────────────────────────────────────────────────────────────────────────

function add_jacobian!(h::ShootingMayerHandle,
                       var::DirectSolverVariable, fn::Function)
    zero_arg = () -> fn(copy(var.value))
    add_objective_jacobian!(h.internal, var, zero_arg)
    return nothing
end

add_jacobian!(fn::Function, h::ShootingMayerHandle,
              var::DirectSolverVariable) =
    add_jacobian!(h, var, fn)

# ─────────────────────────────────────────────────────────────────────────────
# Sequence(phase) constructor for AbstractShootingPhase
# ─────────────────────────────────────────────────────────────────────────────

function Sequence(phase::AbstractShootingPhase)
    seq = Sequence()
    add_sequence!(seq, phase)
    return seq
end
