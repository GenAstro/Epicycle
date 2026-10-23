# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Lagrange and Hermite interpolation of tabulated values, over a window of points around the
# time asked for. Nothing is cached or mutated, so one table can be evaluated from many threads.

"""
    AbstractInterpolationMethod

Supertype of the interpolation methods: [`Lagrange`](@ref) and [`Hermite`](@ref).

# Examples
```jldoctest
Hermite(7) isa AbstractInterpolationMethod

# output

true
```
"""
abstract type AbstractInterpolationMethod end

"""
    Lagrange(degree::Integer = 7)

Lagrange interpolation of the given degree, through `degree + 1` points around the time asked for.

# Fields
- `degree`: Polynomial degree, at least 1.

# Notes
- Values and derivatives are interpolated independently: interpolating positions does not use
  velocities. To interpolate velocities, pass them as values.
- Throws `ArgumentError` if `degree < 1`.

# Examples
```jldoctest
Lagrange(7)

# output

Lagrange(7)
```
"""
struct Lagrange <: AbstractInterpolationMethod
    degree::Int
    function Lagrange(degree::Integer = 7)
        degree >= 1 || throw(ArgumentError("Lagrange: degree must be at least 1; got $degree"))
        return new(degree)
    end
end

"""
    Hermite(degree::Integer = 3)

Hermite interpolation of the given degree, through values and first derivatives at
`(degree + 1) ÷ 2` points around the time asked for.

# Fields
- `degree`: Polynomial degree, odd and at least 1. Each point contributes two conditions, a
  value and a derivative, so degree 7 needs four points.

# Notes
- The interpolated derivative is consistent with the interpolated value: it is the derivative
  of the same polynomial.
- Throws `ArgumentError` if `degree` is even or less than 1.

# Examples
```jldoctest
Hermite(7)

# output

Hermite(7)
```
"""
struct Hermite <: AbstractInterpolationMethod
    degree::Int
    function Hermite(degree::Integer = 3)
        (degree >= 1 && isodd(degree)) || throw(ArgumentError(
            "Hermite: degree must be odd and at least 1, since each point gives a value and a " *
            "derivative; got $degree"))
        return new(degree)
    end
end

"""
    required_points(method::AbstractInterpolationMethod)

The number of tabulated points `method` uses for one evaluation.

# Arguments
- `method`: A [`Lagrange`](@ref) or [`Hermite`](@ref) method.

# Returns
`degree + 1` for Lagrange; `(degree + 1) ÷ 2` for Hermite.

# Examples
```jldoctest
(required_points(Lagrange(7)), required_points(Hermite(7)))

# output

(8, 4)
```
"""
required_points(m::Lagrange) = m.degree + 1
required_points(m::Hermite) = (m.degree + 1) ÷ 2

"""
    interpolate(method::Lagrange, times, values, t)
    interpolate(method::Hermite, times, values, derivatives, t)

The value at `t` of the polynomial through the tabulated points nearest `t`.

# Arguments
- `method`: [`Lagrange`](@ref) or [`Hermite`](@ref).
- `times`: Strictly increasing sample times, in any consistent unit.
- `values`: The values at `times`: numbers or vectors, such as positions.
- `derivatives`: Hermite only. The derivatives of `values` with respect to time at `times`, in
  value units per time unit, such as velocities.
- `t`: The time to evaluate at, in the unit of `times`.

# Notes
- The window is [`required_points`](@ref)`(method)` consecutive points, centred on `t` where the
  table allows and shifted inward at its ends. With an odd count, the extra point is on the side
  nearer `t`.
- `t` must lie within `times[1]` to `times[end]`; there is no extrapolation. For data in
  segments, such as before and after a maneuver, pass one segment at a time so the window never
  spans a boundary.
- At a tabulated time the tabulated value is returned exactly.
- Differentiable with respect to `t`, `values` and `derivatives` (ForwardDiff). Nothing is
  cached or mutated, so the same table can be evaluated from many threads at once.
- Throws `ArgumentError` if `t` is outside the table, if the table has fewer points than the
  method needs, if the lengths differ, or if the times in the window do not increase.

# Returns
The interpolated value, of the same kind as `values` (a number, or a vector).

# Examples
```jldoctest
times  = [0.0, 1.0, 2.0, 3.0]
values = times .^ 3
interpolate(Lagrange(3), times, values, 1.5)

# output

3.375
```
"""
function interpolate(m::Lagrange, times::AbstractVector, values::AbstractVector, t::Real)
    first, n = _window(m, times, values, nothing, t)
    acc = _lagrange_weight(times, first, n, 1, t) * values[first]
    for j in 2:n
        acc += _lagrange_weight(times, first, n, j, t) * values[first + j - 1]
    end
    return acc
end

function interpolate(m::Hermite, times::AbstractVector, values::AbstractVector,
                     derivatives::AbstractVector, t::Real)
    return Base.first(_hermite(m, times, values, derivatives, t))
end

"""
    interpolate_with_derivative(method::Lagrange, times, values, t)
    interpolate_with_derivative(method::Hermite, times, values, derivatives, t)

The value at `t` and its derivative with respect to time, from the same polynomial as
[`interpolate`](@ref).

# Arguments
As for [`interpolate`](@ref).

# Notes
- The derivative is in value units per unit of `times`.
- For Lagrange, the derivative comes from differentiating the polynomial through the values.
  When the table also holds derivatives (velocities for positions), interpolating those directly
  is usually more accurate.

# Returns
A tuple `(value, derivative)`.

# Examples
```jldoctest
times  = [0.0, 1.0, 2.0]
values = times .^ 3
rates  = 3 .* times .^ 2
interpolate_with_derivative(Hermite(3), times, values, rates, 1.5)

# output

(3.375, 6.75)
```
"""
function interpolate_with_derivative(m::Lagrange, times::AbstractVector, values::AbstractVector, t::Real)
    first, n = _window(m, times, values, nothing, t)
    value = _lagrange_weight(times, first, n, 1, t) * values[first]
    rate = _lagrange_weight_derivative(times, first, n, 1, t) * values[first]
    for j in 2:n
        y = values[first + j - 1]
        value += _lagrange_weight(times, first, n, j, t) * y
        rate += _lagrange_weight_derivative(times, first, n, j, t) * y
    end
    return value, rate
end

interpolate_with_derivative(m::Hermite, times::AbstractVector, values::AbstractVector,
                            derivatives::AbstractVector, t::Real) =
    _hermite(m, times, values, derivatives, t)

# --- window ---------------------------------------------------------------------------------

# The value of t without any derivative parts. Which points are used must not depend on the
# partials of a dual number; the ForwardDiff extension adds the method for its duals.
_primal(x::Real) = x

function _window(m, times, values, derivatives, t)
    N = length(times)
    length(values) == N || throw(ArgumentError(
        "interpolate: $N times but $(length(values)) values"))
    derivatives === nothing || length(derivatives) == N || throw(ArgumentError(
        "interpolate: $N times but $(length(derivatives)) derivatives"))
    n = required_points(m)
    N >= n || throw(ArgumentError(
        "interpolate: $(nameof(typeof(m))) degree $(m.degree) needs $n points; the table has $N"))
    tv = _primal(t)                                        # decisions use the value alone
    times[begin] <= tv <= times[end] || throw(ArgumentError(
        "interpolate: t = $t is outside the table, which covers $(times[begin]) to $(times[end]); " *
        "there is no extrapolation"))

    o = firstindex(times) - 1
    k = clamp(searchsortedlast(times, tv) - o, 1, N - 1)      # times[k] <= t <= times[k+1]
    if iseven(n)
        first = k - n ÷ 2 + 1
    else
        nearer = (tv - times[k + o]) <= (times[k + 1 + o] - tv) ? k : k + 1
        first = nearer - (n - 1) ÷ 2
    end
    first = clamp(first, 1, N - n + 1) + o
    for i in first + 1:first + n - 1
        times[i] > times[i - 1] || throw(ArgumentError(
            "interpolate: times must increase strictly; times[$(i - 1)] = $(times[i - 1]) and " *
            "times[$i] = $(times[i])"))
    end
    return first, n
end

# --- Lagrange -------------------------------------------------------------------------------

# ℓⱼ(t) = Πₖ≠ⱼ (t − tₖ)/(tⱼ − tₖ) over the window. At a node every weight is exactly 0 or 1.
function _lagrange_weight(times, first, n, j, t)
    tj = times[first + j - 1]
    w = one(t - tj)
    for k in 1:n
        k == j && continue
        tk = times[first + k - 1]
        w *= (t - tk) / (tj - tk)
    end
    return w
end

# ℓⱼ'(t) = Σₘ≠ⱼ [1/(tⱼ − tₘ)] Πₖ≠ⱼ,ₘ (t − tₖ)/(tⱼ − tₖ), which stays finite at the nodes.
function _lagrange_weight_derivative(times, first, n, j, t)
    tj = times[first + j - 1]
    total = zero(t - tj)
    for m in 1:n
        m == j && continue
        tm = times[first + m - 1]
        term = one(t - tj) / (tj - tm)
        for k in 1:n
            (k == j || k == m) && continue
            tk = times[first + k - 1]
            term *= (t - tk) / (tj - tk)
        end
        total += term
    end
    return total
end

# --- Hermite --------------------------------------------------------------------------------

# Newton divided differences on doubled nodes z = (t₁, t₁, t₂, t₂, …): a repeated node takes the
# tabulated derivative in place of a first difference. Evaluated with its derivative by Horner.
function _hermite(m::Hermite, times, values, derivatives, t)
    first, n = _window(m, times, values, derivatives, t)
    # At a node, return the tabulated value and derivative exactly. Only for a plain float: a
    # ForwardDiff dual must go through the polynomial, or its derivative with respect to t is lost.
    if t isa AbstractFloat
        for i in first:first + n - 1
            times[i] == t && return values[i] * one(t), derivatives[i] * one(t)
        end
    end
    p = 2n
    z = [times[first + (i - 1) ÷ 2] for i in 1:p]
    V = promote_type(eltype(values), eltype(derivatives))
    c = Vector{V}(undef, p)
    for i in 1:p
        c[i] = values[first + (i - 1) ÷ 2]
    end
    for i in p:-1:2                                         # first differences
        c[i] = isodd(i) ? (c[i] - c[i - 1]) / (z[i] - z[i - 1]) : derivatives[first + (i - 1) ÷ 2]
    end
    for order in 2:p - 1                                    # higher differences
        for i in p:-1:order + 1
            c[i] = (c[i] - c[i - 1]) / (z[i] - z[i - order])
        end
    end
    value = c[p] * one(t - z[1])
    rate = zero(value)
    for i in p - 1:-1:1
        rate = rate * (t - z[i]) + value
        value = value * (t - z[i]) + c[i]
    end
    return value, rate
end
