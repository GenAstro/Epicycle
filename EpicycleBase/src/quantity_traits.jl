# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Traits on quantity functions.
#
# A quantity is a plain function, such as `semi_major_axis(sat)`, and everything it
# needs to carry beyond its value is declared here rather than stored in a
# wrapper type. A Julia function is a singleton with a type, so a trait on
# `typeof(f)` dispatches exactly as one on a struct would, and resolves at
# compile time.
#
#     label(::typeof(semi_major_axis))                    = "Semi-major axis"
#     set_quantity!(sat, ::typeof(semi_major_axis); to)   = …
#     is_cyclic(::typeof(raan))                           = true
#     cycle(::typeof(raan))                               = 2π
#
# Declaring these is the whole opt-in. A user's own quantity becomes
# indistinguishable from a shipped one by adding these methods. There is no
# registry and nothing to subtype.
#
# Every trait has an `::Any` default, so a function that declares none still
# works everywhere it is only read.
# =============================================================================

"""
    label(quantity) -> String

What to call this quantity in a report or an error message.

Domain language, capitalised as a title: `"Semi-major axis"`, not `"sma"` and
not `"semi-major axis (km)"`. Units belong to the value, not the name.

Defaults to `"quantity"`, so an undeclared function still prints something.

# Returns
The label as a `String`.

# Examples
```julia
using EpicycleBase
beta_angle(sat) = 0.0                                    # a quantity of your own
EpicycleBase.label(::typeof(beta_angle)) = "Beta angle"  # declaring its label
label(beta_angle)                                        # "Beta angle"
```
"""
label(::Any) = "quantity"

"""
    set_quantity!(subject, quantity, deps...; to)

Write `to` as the value of `quantity` on `subject`. The generic writer a solver calls, because a
solver holds the quantity as a value and cannot name its setter.

A user writes the quantity's own setter, its name with a bang: `semi_major_axis!(sat; to = v)`.
A settable quantity defines that setter and one method of this function that calls it, so the
two forms reach the same code.

Writing dispatches on the subject and the quantity together. A quantity is not settable in the
abstract; it is settable on a particular kind of subject, and only a method on the pair can say
so. There is no table mapping a quantity to a writer.

# Arguments
- `subject`: what the quantity is read from and written to.
- `quantity`: the quantity's reader function.
- `deps`: the extra arguments the reader takes, such as a coordinate system, in the same order.
  Declare them with defaults, one arity at a time, rather than as `deps...`, so that
  [`is_settable`](@ref) answers only for argument shapes the method accepts.
- `to`: the value to write, in the units the reader returns.

# Returns
Whatever the method returns; callers use it for its effect on `subject`, which it mutates.

# Example
```jldoctest
mutable struct Burn
    dv::Float64
end
burn_size(b::Burn) = b.dv
burn_size!(b::Burn; to::Real) = (b.dv = to; b)                        # what a user writes
EpicycleBase.set_quantity!(b::Burn, ::typeof(burn_size); to) = burn_size!(b; to = to)

b = Burn(0.1)
set_quantity!(b, burn_size; to = 0.25)                                # what a solver writes
burn_size(b)

# output
0.25
```
"""
function set_quantity! end

"""
    is_settable(subject, quantity, deps...) -> Bool

Whether this quantity can be written on this subject. Derived from whether a
[`set_quantity!`](@ref) method takes these positional arguments, so there is no second place to
declare it and no way for the two to disagree.

Resolves at compile time without evaluating anything, so a solver can reject an unsettable
variable when it is declared rather than at the fortieth iteration, and pays nothing to ask again
on every write. It does not check for the keyword `to`, because a keyword probe cannot be resolved
at compile time; a method that omits `to` fails on the first write instead.

# Returns
`true` when a `set_quantity!` method matches the argument types, `false` otherwise.

# Example
```jldoctest
mutable struct SettableValue
    value::Float64
end
example_value(::SettableValue) = 0.0
EpicycleBase.set_quantity!(subject::SettableValue, ::typeof(example_value); to) =
    (subject.value = to; nothing)
is_settable(SettableValue(0.0), example_value)

# output
true
```
"""
is_settable(subject, quantity, deps...) =
    hasmethod(set_quantity!, Tuple{typeof(subject), typeof(quantity), map(typeof, deps)...})

"""
    output_partial(subject, quantity, deps...) -> AbstractMatrix

The derivative of this quantity with respect to the subject's state, one row per
component of the quantity and one column per component of the state.

This is the derivative a solver cannot work out for itself and the quantity
can. Semi-major axis is a function of position and velocity and of nothing else,
so the six numbers here depend on nothing outside the quantity. What the solver
actually needs is the derivative with respect to whatever it is varying, and it
gets there by the chain rule, multiplying this by the state transition matrix the
propagation already produced. A quantity that tried to supply that product
instead would have to know which solver variable was upstream of it, which is
exactly the coupling this split keeps out of a user's code.

Declared as a method, the same shape as [`set_quantity!`](@ref), so a reader
learns one rule for both. Deps are the extra arguments the reader took and are
matched by type, so a quantity read in a named frame declares a different method
from the same quantity read in the subject's own.

Not every quantity has one. A gravitational parameter is a scalar field on a
model rather than a function of the state; it declares [`tag`](@ref) instead and
its derivative comes from `param_jac!`. And a quantity that declares neither is
differentiated, which is correct and slower and requires its reader to accept
dual numbers.

# Returns
A matrix with one row per component of the quantity and one column per state component.

# Example
```jldoctest
struct PointMass
    r::Vector{Float64}
end
radius(p::PointMass) = sqrt(sum(abs2, p.r))
EpicycleBase.output_partial(p::PointMass, ::typeof(radius)) = reshape(p.r ./ radius(p), 1, :)
J = output_partial(PointMass([3.0, 4.0, 0.0]), radius)
size(J), J[1, 2]

# output
((1, 3), 0.8)
```
"""
function output_partial end

"""
    has_output_partial(subject, quantity, deps...) -> Bool

Whether this quantity can hand a solver its own derivative on this subject.
Derived from whether an [`output_partial`](@ref) method exists, so there is no
second place to declare it and no way for the two to disagree.

Resolves without evaluating anything, so the sparsity pattern can be built
before the first iteration rather than discovered during it.

# Returns
`true` when an `output_partial` method matches the argument types, `false` otherwise.

# Example
```jldoctest
struct PartialSubject end
partial_quantity(::PartialSubject) = 0.0
EpicycleBase.output_partial(::PartialSubject, ::typeof(partial_quantity)) = zeros(1, 6)
has_output_partial(PartialSubject(), partial_quantity)

# output
true
```
"""
has_output_partial(subject, quantity, deps...) =
    hasmethod(output_partial, Tuple{typeof(subject), typeof(quantity),
                                    map(typeof, deps)...})

"""
    is_cyclic(quantity) -> Bool

Whether this quantity wraps, as an angle does and a distance or an energy does not.

A stopping condition on `raan = 0` brackets a sign change across a step, and at
the 2π wrap boundary the residual jumps by a full turn, so it either misses the
crossing or converges onto the discontinuity. This trait is how a stopping
condition would know to unwrap first. AstroProp does not read it yet, so
stopping on an angle is not supported.

See also [`cycle`](@ref), which gives the period.

# Returns
`true` for a quantity that has declared itself cyclic; `false` by default.

# Examples
```julia
using EpicycleBase
is_cyclic(sin)    # false: nothing has declared it
```
"""
is_cyclic(::Any) = false

"""
    cycle(quantity) -> Real or nothing

The period over which a cyclic quantity repeats, such as `2π` for an angle in
radians.

# Returns
The period, in the quantity's own units, or `nothing` when [`is_cyclic`](@ref) is `false`.

# Example
```jldoctest
cycle(sin) === nothing

# output
true
```
"""
cycle(::Any) = nothing

"""
    tag(quantity) -> AbstractVarTag or nothing

The variable tag this quantity corresponds to, when it corresponds to one.

This is the join between a quantity and the analytic partials already
registered in `AstroProp`, which are selected by a [`ModelVariable`](@ref)'s tag:
`param_jac!(out, force, mv.tag, t, y, sc)`. A quantity that is a genuine
scalar field on a model, such as a gravitational parameter, has a tag. One
derived from a state representation, such as semi-major axis, does not: its
partial is `∂g/∂x`, an output partial, which is a different derivative from
`∂f/∂p` and does not live here.

# Returns
The quantity's `AbstractVarTag` singleton, or `nothing` by default.

# Example
```jldoctest
tag(sin) === nothing

# output
true
```
"""
tag(::Any) = nothing
