# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# The element extension interface, the abstract array types and the tag structs.
#

# ─────────────────────────────────────────────────────────────────────────────
# The element extension interface
#
# Sequence assembles the NLP from its elements, and an element answers only for
# what it can see. That division is the point: an element knows its own
# contents, and Sequence knows what spans elements. Neither can do the other's
# half, so neither can be written in terms of the other.
#
# Four methods. A phase type that supplies them is assembled like any other,
# and nothing in Sequence names it.
#
# These are declared rather than left implicit because five phase types
# already implement them and none of them was written against a stated
# contract. Without the declarations a third party had nothing to implement
# against.
# ─────────────────────────────────────────────────────────────────────────────

"""
    variable_list(element) -> Vector

The solver variables this element contributes, in order.

Order fixes the column layout, so it must not change between calls. A variable
that appears in two elements' lists is one column that both reach, and neither
element knows the other exists; Sequence resolves that when it composes the
global index.
"""
function variable_list end

"""
    nlp_length(element, var) -> Int

How many NLP columns one of this element's variables occupies.

A state block over a mesh occupies one column per state per node; a scalar
parameter occupies one. Asked once per variable when the layout is built, so it
may be computed rather than stored.
"""
function nlp_length end

"""
    function_list(element) -> Vector{PhaseFunction}

The NLP functions this element contributes, each carrying its row count and
bounds.

Order fixes the row layout, the same way `variable_list` fixes the columns.
"""
function function_list end

"""
    sparsity_structure(element, f) -> Vector{Bool}
    sparsity_structure(element)    -> Matrix{Bool}

Which of this element's own variables the function `f` can touch, in
`variable_list` order. The one-argument form assembles the whole block by
calling the two-argument form for every function.

Derive this from the evaluation context the element hands `f` rather than from
the variable's type. The element builds that context, so it knows what `f` can
read, and a function depends on what it is given and on nothing else. That
makes the context's fields an upper bound which is exact at the block level and
needs no declaration from the user.

**Over-declare rather than drop.** A caller skips computing any block this
returns `false` for, so a dropped dependence never reaches the solver and
nothing downstream can tell a derivative that is zero from one that was never
asked for. An over-declared block is computed and found to be zero, which costs
time and nothing else.

Answers for this element alone. Coupling that spans elements is Sequence's, and
an element that tried to report it would have to know which other elements
exist.
"""
function sparsity_structure end


# ─────────────────────────────────────────────────────────────────────────────
# Abstract array types
# ─────────────────────────────────────────────────────────────────────────────

abstract type AbstractStateArray{V<:AbstractState} end
abstract type AbstractControlArray{V<:AbstractControl} end
abstract type AbstractParameter end

# AbstractModel — marker supertype only.  T lives on the concrete struct, not here.
abstract type AbstractModel end

# The extension interface for a shooting transcription. A package providing one
# subtypes this, and the sequence machinery dispatches on it without knowing
# what the transcription is. Declared here rather than with the shooting code so
# that a sequence can name the type with no such package loaded.
abstract type AbstractShootingPhase end

# eltype_state lives next to AbstractStateArray, which is the type it reads from.
eltype_state(::AbstractStateArray{S}) where {S} = S

# ─────────────────────────────────────────────────────────────────────────────
# Tag structs, defined before the add_jacobian! methods that dispatch on them
# ─────────────────────────────────────────────────────────────────────────────

# Path-level (per-node)
struct State     end
struct Control   end
struct Parameter end
struct TimeTag   end   # not `Time`, which is AstroEpochs' epoch type and
                       # would collide the moment a script uses both

# Sentinel key for time-partial registrations.
# Stored inside the same Dict{UInt64,Function} as variable-keyed entries.
# objectid(TimeTag) is a DataType objectid — distinct from any DirectSolverVariable's
# objectid, which is always a heap-allocated struct instance.
const _TIME_JAC_KEY = objectid(TimeTag)

# Boundary-level (endpoints)
struct InitialState end
struct FinalState   end
struct InitialTime  end
struct FinalTime    end

# Objective sense
"""
    Min()

The `sense` that asks an `Objective` to be minimized, and the default:
`Objective(f, phase; sense = Min())`. See also [`Max`](@ref).

# Example
<!-- doc-fragment -->
```julia
Objective(tracking_error, phase; sense = Min(), at = Path())
```
"""
struct Min end

"""
    Max()

The `sense` that asks an `Objective` to be maximized, as in `Objective(final_mass, phase; sense
= Max())`. The solver minimizes the negated value. See also [`Min`](@ref).

# Example
<!-- doc-fragment -->
```julia
Objective(final_mass, phase; sense = Max())
```
"""
struct Max end

# ─────────────────────────────────────────────────────────────────────────────
# Transcription
# ─────────────────────────────────────────────────────────────────────────────

# The transcription extension interface is in src/transcription.jl. It is
# twelve names: EvalContext, four generic functions on the transcription type,
# and seven on the mesh it builds. Nothing else about a transcription is known
# here, which is what lets an implementation live in another package on either
# side of the open/Enterprise boundary.
#
# Hermite-Simpson is the open tier's reference transcription and lives in this
# package, beside the interface it implements. Enterprise transcriptions reach
# the same generic functions from their own package, which is what makes this an
# interface rather than a description of what happens to be in one file.

# ─────────────────────────────────────────────────────────────────────────────
# BoundaryContext
#
# Carries y0, yf, t0, tf, params for boundary function and Mayer objective
# evaluation.  Each field has its own type parameter so ForwardDiff can seed
# any one of them with Dual numbers while the rest stay Float64.
#
# Convention (single-phase):
#   y0     = p._Y[:, 1]     initial state at first collocation node
#   yf     = p._Y[:, end]   final state at last collocation node
#   t0, tf = p._t0, p._tf   current phase times
#   params = p._params       free parameter vector
#
# For multi-phase linkage the framework will build one BoundaryContext per
# participating phase and pass them as separate arguments.
# ─────────────────────────────────────────────────────────────────────────────

struct BoundaryContext{Y0, Yf, T0<:Number, Tf<:Number, P<:AbstractVector}
    y0::Y0
    yf::Yf
    t0::T0
    tf::Tf
    params::P
end

BoundaryContext(y0, yf, t0, tf, params) =
    BoundaryContext{typeof(y0), typeof(yf), typeof(t0), typeof(tf), typeof(params)}(
        y0, yf, t0, tf, params)

BoundaryContext(y0, yf, t0, tf) = BoundaryContext(y0, yf, t0, tf, Float64[])

# ─────────────────────────────────────────────────────────────────────────────
# State named-tuple helpers
#
# _state_named constructs the user's parametric state struct from a plain vector,
# preserving the element type so ForwardDiff Dual numbers flow through cleanly.
# All user state types satisfy the parametric contract (S{T} <: AbstractState),
# so S is always a UnionAll and S(v...) always works.
#
# ─────────────────────────────────────────────────────────────────────────────

function _state_named(p, v::AbstractVector)
    S = eltype_state(p.state_var.var)
    S(v...)
end
