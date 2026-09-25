# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# DirectSolverVariable, and the keyword split that keeps AstroSolve's SolverVariable.
#

# ─────────────────────────────────────────────────────────────────────────────
# SolverVariable: one keyword constructor for two kinds of variable.
# The event-graph form takes calc= and singular bounds; the transcription form
# takes var= and plural bounds.
# ─────────────────────────────────────────────────────────────────────────────

mutable struct DirectSolverVariable
    var          ::Any
    lower_bounds ::Vector{Float64}
    upper_bounds ::Vector{Float64}
    scale        ::Vector{Float64}   # x_nlp = (x_phys - shift) / scale  (default: 1)
    shift        ::Vector{Float64}   # x_nlp = (x_phys - shift) / scale  (default: 0)
    name         ::String
    value        ::Vector{Float64}   # current physical value (always dimensional)
end

# Two kinds of variable want this name. An event graph varies a quantity on a
# subject, given as `calc` with singular bounds; a transcribed phase varies a
# block of the decision vector, given as `var` with plural bounds. Julia
# dispatches on positional arguments and both forms take none, so they cannot be
# two methods: this is the one method, and it decides from the keywords.
#
# Neither should be a name a user types. `Vary` is the verb for the first and
# `set_state!`, `set_control!` and `set_parameter!` for the second, but call
# sites still say SolverVariable.
function SolverVariable(; var=nothing, calc=nothing, name::String="",
                          lower_bounds=nothing, upper_bounds=nothing,
                          lower_bound=nothing, upper_bound=nothing,
                          scale=nothing, shift=nothing,
                          value=nothing,
                          role::AbstractRole = SolveFor(),
                          covariance=nothing, process_noise=nothing)
    if var === nothing && calc !== nothing &&
       lower_bounds === nothing && upper_bounds === nothing
        # The parameterized inner constructor, not the keyword outer one: this
        # method extends the same generic, so calling the outer form recurses.
        return SolverVariable{typeof(calc)}(; calc = calc, name = name,
                                                       lower_bound = lower_bound,
                                                       upper_bound = upper_bound,
                                                       scale = scale, shift = shift,
                                                       role = role,
                                                       covariance = covariance,
                                                       process_noise = process_noise)
    end
    # A transcribed block is bounded, never estimated, so the estimation
    # keywords have no meaning on this branch.
    (covariance === nothing && process_noise === nothing && role === SolveFor()) ||
        throw(ArgumentError(
            "SolverVariable: role, covariance and process_noise belong to an " *
            "estimated quantity, which is the calc form. This is the var form, " *
            "which a transcription bounds and scales."))
    n  = length(something(lower_bounds, upper_bounds))
    sc = scale === nothing ? ones(n)  : Float64.(scale)
    sh = shift === nothing ? zeros(n) : Float64.(shift)
    # default initial value: the shift (= physical zero of the NLP parameterisation)
    v0 = value === nothing ? copy(sh) : Float64.(value)
    return DirectSolverVariable(
        something(var, calc),
        Float64.(lower_bounds),
        Float64.(upper_bounds),
        sc, sh,
        name,
        v0,
    )
end
