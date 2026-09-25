# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# `Calc` — a quantity, unapplied.
#
#     raan(sat, cs)          the value now
#     Calc(raan, sat, cs)    the same tokens, deferred
#
# Two things make it worth having rather than reaching for a closure.
#
# It is a **zero-arg callable**, so anything that accepts "a callable taking no
# arguments" accepts it with no new protocol: a `StopAt` here, or a solver
# variable or a `Constraint` in AstroSolve.
#
# And it **keeps its parts**, so it can be re-applied to a different subject.
# That is what a history walk needs: every recorded sample is a snapshot of the
# subject, and the walk evaluates the same quantity against each in turn. A
# closure cannot do this — `() -> raan(sat, cs)` has `sat` sealed inside and
# will return the same number for every sample.
#
# So: **structured where it must be re-applied, opaque where it need not be.**
# `history` needs structure. A constraint does not — it is evaluated against
# whatever the sequence has done to the object, never re-pointed at another one
# — which is why a closure is perfectly legal there and an expression like
# `sma(a) - sma(b)` can only ever be a closure.
# =============================================================================

"""
    Calc(quantity, subject, deps...)

A quantity together with everything it needs, held unapplied.

Call it with no arguments to evaluate. `Calc` carries the subject, so nothing
is left to context. Reporting one quantity for two spacecraft uses two `Calc`s,
each naming its own.

# Fields
- `f`: Quantity function.
- `args`: Subject followed by any dependencies, such as a coordinate system.

# Examples
```julia
q = Calc(raan, sat, EarthMJ2000Ec)
q()                                    # evaluate now

history(q)                             # the series over sat's recorded samples
AstroSolve.Constraint(q; equals = 0.0) # a constraint on it, from AstroSolve
```

Several columns in different frames come from one walk of the history, because
each `Calc` carries its own dependency:

```julia
t, r_eq, r_ec = history(Calc(epoch,           sat),
                        Calc(position_vector, sat, EarthMJ2000Eq),
                        Calc(position_vector, sat, EarthMJ2000Ec))
```

See also [`reapply`](@ref).
"""
struct Calc{F,A<:Tuple}
    f::F
    args::A
end

Calc(f, args...) = Calc(f, args)

(c::Calc)() = c.f(c.args...)

"""
    reapply(c::Calc, subject) -> value

Evaluate `c` against a different subject, keeping its dependencies.

This is what a history walk does with each recorded sample. It works because a
`Calc` keeps its parts; there is no equivalent for a closure.

# Returns
The value of the stored quantity evaluated on `subject`.

# Example
```julia
q = Calc(position_magnitude, sat)
reapply(q, another_sat)
```
"""
reapply(c::Calc, subject) = c.f(subject, Base.tail(c.args)...)

"""
Evaluating a `Calc` is calling it, so it works anywhere `get_calc` is used —
which is what lets a stopping condition take one with no other change.
"""
get_calc(c::Calc) = c()

"""
Write a value back through the quantity's own setter.

This is the half of the interface a solver needs and a report does not. There
is no second vocabulary for it: reading dispatches on the subject and so does
writing, through `set_quantity!(subject, quantity, deps...; to = value)`, so
`Vary(delta_v, toi)` names the quantity once and gets both directions.
"""
function set_calc!(c::Calc, value)
    EpicycleBase.is_settable(first(c.args), c.f, Base.tail(c.args)...) ||
        throw(ArgumentError(
            "$(EpicycleBase.label(c.f)) cannot be set on a " *
            "$(nameof(typeof(first(c.args)))); a solver can vary only a quantity " *
            "that has a `set_quantity!` method for its subject."))
    return EpicycleBase.set_quantity!(first(c.args), c.f, Base.tail(c.args)...; to = value)
end

"""
    calc_numvars(calc) -> Int

Number of scalar values produced by a calculation.

# Returns
The vector length, or `1` for a scalar quantity.

# Example
```julia
calc_numvars(Calc(position_vector, sat))
```
"""
calc_numvars(c::Calc) = (v = c(); v isa AbstractVector ? length(v) : 1)

"""
    calc_is_settable(calc) -> Bool

Whether the calculation can write a value back to its subject.

# Returns
`true` when a matching setter exists, otherwise `false`.

# Example
```julia
calc_is_settable(Calc(delta_v, toi))
```
"""
calc_is_settable(c::Calc) =
    EpicycleBase.is_settable(first(c.args), c.f, Base.tail(c.args)...)

# The same two questions asked of the older stack, so a caller need not know
# which kind of calc it holds.
calc_numvars(c::AbstractCalc) =
    Base.hasproperty(c, :var) ? calc_numvars(getproperty(c, :var)) : 1
calc_is_settable(c::AbstractCalc) =
    Base.hasproperty(c, :var) ? calc_is_settable(getproperty(c, :var)) : false

"""The quantity's own label, so a `Calc` names itself in a report."""
EpicycleBase.label(c::Calc) = EpicycleBase.label(c.f)

"""Traits pass through to the quantity, so a `Calc` on an angle remains cyclic."""
EpicycleBase.is_cyclic(c::Calc) = EpicycleBase.is_cyclic(c.f)
EpicycleBase.cycle(c::Calc)     = EpicycleBase.cycle(c.f)
EpicycleBase.tag(c::Calc)       = EpicycleBase.tag(c.f)

function Base.show(io::IO, c::Calc)
    print(io, "Calc(", EpicycleBase.label(c.f))
    length(c.args) > 1 && print(io, ", ", join(map(_show_arg, Base.tail(c.args)), ", "))
    print(io, ")")
end

_show_arg(x) = string(x)
_show_arg(cs::AbstractCoordinateSystem) = string(nameof(typeof(cs.axes)), "@", _origin_name(cs))
_origin_name(cs) = hasproperty(cs.origin, :name) ? String(cs.origin.name) : "origin"
