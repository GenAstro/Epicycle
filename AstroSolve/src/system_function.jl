# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# A function bundled with the analytic Jacobians someone declared for it.
#
# An estimator is handed dynamics and a measurement model and has to ask each of
# them what it knows how to differentiate, falling back to automatic
# differentiation where the answer is nothing. That question is the same one
# add_jacobian! and has_jacobian already answer for a constraint or an
# objective, so these are methods on those verbs rather than a second set.

abstract type AbstractSystemFunction end

"""
    MeasurementFunction(f; name = :meas)
    MeasurementFunction(name = :meas) do y, u, p, t, model … end

A measurement function and its registered analytic Jacobians.

# Fields
- `f::F`: Measurement closure.
- `jacs::Dict{Any,Any}`: Analytic Jacobians keyed by differentiation tag.
- `name::Symbol`: Name used in diagnostics.

# Notes
The closure takes the universal signature `f(y, u, p, t, model)` and returns the predicted
measurement as a vector. Every estimator in this package calls it with all five arguments.

# Example
```julia
using AstroSolve
measurement = AstroSolve.MeasurementFunction((y, u, p, t, model) -> [y[1]]; name = :position)
measurement([2.0], (;), (;), 0.0, nothing)
```
"""
struct MeasurementFunction{F} <: AbstractSystemFunction
    f    :: F
    jacs :: Dict{Any,Any}
    name :: Symbol
end

MeasurementFunction(f; name::Symbol = :meas) =
    MeasurementFunction(f, Dict{Any,Any}(), name)

"""
    DynamicsFunction(f; name = :dyn)
    DynamicsFunction(name = :dyn) do y, u, p, t, model … end

A dynamics function and its registered analytic Jacobians.

# Fields
- `f::F`: Dynamics closure.
- `jacs::Dict{Any,Any}`: Analytic Jacobians keyed by differentiation tag.
- `name::Symbol`: Name used in diagnostics.

# Notes
The closure takes the universal signature `f(y, u, p, t, model)` and returns the state rate as a
vector. A Jacobian registered with [`add_jacobian!`](@ref) takes the same five arguments.

# Example
```julia
using AstroSolve
dynamics = AstroSolve.DynamicsFunction((y, u, p, t, model) -> [y[2], -y[1]]; name = :oscillator)
dynamics([1.0, 0.0], (;), 0.0)
```
"""
struct DynamicsFunction{F} <: AbstractSystemFunction
    f    :: F
    jacs :: Dict{Any,Any}
    name :: Symbol
end

DynamicsFunction(f; name::Symbol = :dyn) =
    DynamicsFunction(f, Dict{Any,Any}(), name)

# Make any system function callable: sf(state, p, t) → out
(sf::AbstractSystemFunction)(args...; kwargs...) = sf.f(args...; kwargs...)

"""
    add_jacobian!(jac, sf, tag)
    add_jacobian!(sf, tag) do y, u, p, t, model … end

Register an analytic Jacobian for one argument of a system function.

# Arguments
- `jac`: Jacobian closure with the same arguments as the wrapped function.
- `sf::AbstractSystemFunction`: Measurement or dynamics function to extend.
- `tag`: Differentiation tag. The three are `State()`, `Control()` and `Parameter()`.

# Returns
The modified system function, so registrations can be chained.

# Example
```julia
using AstroSolve
dynamics = AstroSolve.DynamicsFunction((y, u, p, t, model) -> [y[2], -y[1]])
AstroSolve.add_jacobian!((y, u, p, t, model) -> [0.0 1.0; -1.0 0.0],
                         dynamics, AstroSolve.State())
```
"""
function add_jacobian!(jac, sf::AbstractSystemFunction, tag)
    sf.jacs[tag] = jac
    return sf
end

"Look up the analytic Jacobian closure registered on `sf` for `tag`."
get_jacobian(sf::AbstractSystemFunction, tag) = sf.jacs[tag]

"Has an analytic Jacobian been registered on `sf` for `tag`?"
has_jacobian(sf::AbstractSystemFunction, tag) = haskey(sf.jacs, tag)
