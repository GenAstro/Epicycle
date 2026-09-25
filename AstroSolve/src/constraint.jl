# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Constraint — require a quantity to reach a value
#
# `Constraint` and `Vary` are the two verbs a user writes side by side, and they
# take the same first arguments: the quantity, its subject, then whatever the
# quantity depends on. They belong in the same package for that reason.
#
# AstroCallbacks evaluates quantities; AstroSolve applies solver bounds to them.
#
# This lived in AstroCallbacks. Deciding what a solver does with a quantity is this package's job,
# but the calc machinery it builds on is still AstroCallbacks' and is reached through that package's
# public names rather than its internals.
# =============================================================================

"""
    Constraint(quantity, subject, deps...; equals)
    Constraint(quantity, subject, deps...; lower_bound, upper_bound, scale)

Require a quantity to reach a value.

Built by the constructors below rather than field by field. The fields are listed because a solved
sequence hands these back and a report reads them.

# Fields
- `calc::C`: Quantity, subject, and dependencies evaluated at each iteration.
- `lower_bound::Vector{T}`: Lower bound for each quantity component.
- `upper_bound::Vector{T}`: Upper bound for each quantity component.
- `scale::Vector{T}`: Solver scale for each quantity component.
- `numvars::Int`: Number of scalar components.
- `name::String`: Name used in reports.

# Notes
[`Vary`](@ref) is the other half of a solver spec.
"""
struct Constraint{C, T<:Real}
    calc::C
    lower_bound::Vector{T}
    upper_bound::Vector{T}
    scale::Vector{T}
    numvars::Int
    name::String        # what to call it when it is reported

    # Inner constructor: single point of truth for validation/invariants
    function Constraint{C,T}(
        calc::C,
        lb::Vector{T},
        ub::Vector{T},
        sc::Vector{T},
        numvars::Integer,
        name::AbstractString = "",
    ) where {C, T<:Real}
        n = Int(numvars)
        if length(lb) != n || length(ub) != n || length(sc) != n
            throw(ArgumentError(
                "Constraint: lower/upper/scale lengths must equal numvars=$(n); " *
                "got lower=$(length(lb)), upper=$(length(ub)), scale=$(length(sc))."
            ))
        end
        return new{C,T}(calc, lb, ub, sc, n, String(name))
    end
end

# TODO: Define a compact display for Constraint.

"""
    Constraint(quantity, subject, deps...; equals)
    Constraint(quantity, subject, deps...; lower, upper, scale = 1.0)

Require a quantity to reach a value.

Its first arguments are the same as every other spec — the quantity, its
subject, then whatever the quantity depends on — and the goal is a keyword, so
it reads like the stopping condition beside it:

# Arguments
- `quantity::Function`: Quantity to constrain.
- `subject`: Object on which the quantity is evaluated.
- `deps...`: Additional dependencies required by the quantity.
- `equals`: Value required for an equality constraint.
- `lower_bound`: Scalar or component-wise lower bound.
- `upper_bound`: Scalar or component-wise upper bound.
- `scale`: Scalar or component-wise solver scale.
- `name::AbstractString`: Name used in reports.
- `at`: Unsupported for event constraints; use it only with a phase constraint.

# Notes
`equals` sets both bounds and cannot be combined with `lower_bound` or
`upper_bound`. Scalars are applied to every component. The constructor throws
an `ArgumentError` when no bound is supplied, when bound sizes do not match the
quantity, or when `at` is given.

# Returns
A registered `Constraint`.

# Example
<!-- doc-fragment -->
```julia
StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
Constraint(position_magnitude, sat; equals = 85000.0)
Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0))
```
"""
function Constraint(quantity::Function, subject, deps...;
                    equals = nothing, lower_bound = nothing, upper_bound = nothing,
                    scale = nothing, name::AbstractString = "", at = nothing)
    # `at` says where a constraint on a phase is evaluated — a path, an
    # endpoint. Here the subject already says: this constrains that maneuver,
    # that spacecraft, at the point the sequence puts it. Accepted so the
    # failure is this sentence rather than a MethodError.
    at === nothing || throw(ArgumentError(
        "Constraint: `at` has no meaning here — the subject you named is the " *
        "place. `at` belongs on a constraint attached to a phase, where a " *
        "path and an endpoint are different points."))
    calc = Calc(quantity, subject, deps...)
    n    = calc_numvars(calc)

    if equals !== nothing
        (lower_bound === nothing && upper_bound === nothing) || throw(ArgumentError(
            "Constraint: give `equals`, or `lower_bound`/`upper_bound`, not both."))
        lower_bound = upper_bound = equals
    end
    (lower_bound === nothing && upper_bound === nothing) && throw(ArgumentError(
        "Constraint on $(EpicycleBase.label(quantity)) needs `equals`, or a " *
        "`lower_bound` or `upper_bound`."))

    _vec(x, default) = x === nothing ? fill(default, n) :
                       (x isa AbstractVector ? Float64.(collect(x)) : fill(Float64(x), n))

    # Straight to the inner constructor, which validates. `_vec` gives all three
    # the same Float64 eltype, so there is nothing to promote. Bounds are
    # constants and never carry a dual number, whatever the state does.
    con = Constraint{typeof(calc), Float64}(
        calc, _vec(lower_bound, -Inf), _vec(upper_bound, Inf),
        _vec(scale, 1.0), n, name)
    return _register_constraint!(con)   # checks the step before it inside a target! block
end

"""
    func_eval(constraint::Constraint)

Evaluate a constraint and preserve the result's element type.

# Arguments
- `constraint::Constraint`: Constraint to evaluate.

# Notes
The element type is carried through rather than converted to `Float64`, so a dual number survives
evaluation and the constraint can be differentiated. Narrowing the return type breaks every
derivative taken through a constraint.

# Returns
A vector containing the scalar or vector value of the constrained quantity.

Throws an `ArgumentError` when the quantity returns neither a number nor an
abstract vector.

# Example
<!-- doc-fragment -->
```julia
func_eval(Constraint(position_magnitude, sat; equals = 7000.0))
```
"""
function func_eval(constraint::Constraint)
    val = get_calc(constraint.calc)
    if isa(val, Number)
        return [val]
    elseif isa(val, AbstractVector)
        return collect(val)
    else
        throw(ArgumentError(
            "func_eval: a constrained quantity must return a Number or an AbstractVector, " *
            "since a constraint is a scalar or a vector of scalars; got $(typeof(val))"))
    end
end
