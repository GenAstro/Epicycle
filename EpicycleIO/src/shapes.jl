# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Turning what the user passed into a list of series, and spreading attributes across them.
#
# This is one of the two things EpicycleIO adds to Plotly. Plotly draws one series
# per trace; `history` hands back a position column as a vector of 3-element vectors, and
# turning that into three traces is a loop the user should not write.

"""One series: paired x and y, both plain vectors, plus the component index it came from."""
struct Series
    # Element type is deliberately left open. An axis is not always numbers — a bar chart's
    # categories are strings, and a Plotly time axis takes date strings — so pinning these to
    # `Vector{Float64}` would refuse plots that work. Nothing here is on a hot path: this is
    # called once per `xyplot`, not per sample, so the missed specialization costs nothing worth
    # having (§14.4).
    x::Vector
    y::Vector
    component::Int      # 0 when the y argument was not multi-component
    ncomponents::Int
end

_isvecvec(v) = v isa AbstractVector && !isempty(v) && all(e -> e isa AbstractVector, v)

"""
    expand_series(args...) -> Vector{Series}

Turn the positional arguments of a plotting call into a list of series.

    xyplot(y)                 index on the horizontal axis
    xyplot(x, y)              one series
    xyplot(x, Y)              one series per column of a matrix, or per component of a
                            vector of vectors
    xyplot(x1, y1, x2, y2)    independent pairs, each with its own length

Lengths are checked pairwise: `x1` against `y1`, `x2` against `y2`. The pairs need nothing to do
with each other, which is why two spacecraft on one plot is not a special case.
"""
function expand_series(args...)
    isempty(args) && throw(ArgumentError(
        "plot needs data. Give it y, or x and y — for example xyplot(t, altitude). A leading " *
        "string names the panel: xyplot(\"Altitude\", t, altitude)."))

    if length(args) == 1
        y = args[1]
        return _pair(collect(1:_ylength(y)), y)
    end

    isodd(length(args)) && throw(ArgumentError(
        "plot got $(length(args)) positional arguments. Give one array, an x and a y, or " *
        "x/y pairs — an odd count above one leaves an array without a partner."))

    series = Series[]
    for i in 1:2:length(args)
        append!(series, _pair(args[i], args[i+1]))
    end
    return series
end

# Length of the y argument measured in samples, whatever its shape.
_ylength(y::AbstractMatrix) = size(y, 1)
_ylength(y) = length(y)

# One x with one y argument, which may itself carry several components.
function _pair(x, y)
    xs = collect(x)
    n = length(xs)

    if y isa AbstractMatrix
        size(y, 1) == n || throw(ArgumentError(
            "x has $n points but the matrix has $(size(y,1)) rows. Columns are series, so " *
            "rows must match x."))
        k = size(y, 2)
        return [Series(xs, collect(view(y, :, j)), j, k) for j in 1:k]

    elseif _isvecvec(y)
        # Sample-major: one entry per sample, each holding the components of that sample.
        # This is what `history` returns for a position column.
        if length(y) == n
            k = length(first(y))
            bad = findfirst(e -> length(e) != k, y)
            bad === nothing || throw(ArgumentError(
                "the components are ragged: sample 1 has $k values but sample $bad has " *
                "$(length(y[bad])). Every sample must carry the same number of components."))
            return [Series(xs, [e[j] for e in y], j, k) for j in 1:k]

        # Series-major: one entry per series, each already a full-length column.
        elseif all(e -> length(e) == n, y)
            k = length(y)
            return [Series(xs, collect(y[j]), j, k) for j in 1:k]

        else
            throw(ArgumentError(
                "x has $n points; the y argument holds $(length(y)) vectors of lengths " *
                "$(unique(length.(y))). Either give one entry per sample, or one full-length " *
                "vector per series."))
        end
    end

    ys = collect(y)
    length(ys) == n || throw(ArgumentError(
        "x has $n points but y has $(length(ys)). They must agree."))
    return [Series(xs, ys, 0, 1)]
end

"""
    plotdata(v)

Make a vector safe to serialise, turning anything not finite into `nothing`.

JSON has no NaN or Infinity — `JSON.json` refuses them outright, so a single gap in a series
used to take the whole call down. Plotly reads `null` as a break in the line, which is what a
gap means anyway: no measurement here, do not draw through it.

Gaps are ordinary in this domain. A quantity between ground station passes, argument of
periapsis on a circular orbit, a solver step that did not converge — none of those should cost
the user their plot.

Only walked when something is actually wrong, so a clean series is passed through untouched.
"""
plotdata(v::AbstractVector{<:Real}) =
    all(isfinite, v) ? v : Any[isfinite(x) ? x : nothing for x in v]
plotdata(v) = v

# ─── Attributes across several traces ─────────────────────────────────────────────────────────

"""
    distribute(key, value, i, ntraces, trace_type) -> value for trace i

Decide what one attribute means when a call expanded into several traces.

A scalar applies to every trace. A vector means one value per trace unless Plotly already reads
an array on that attribute as one value per data point, in which case its meaning is left alone
and per-trace is spelled by nesting.
"""
function distribute(key::Symbol, value, i::Int, ntraces::Int, trace_type::Symbol)
    ntraces == 1 && !(_isvecvec(value)) && return value

    if _isvecvec(value)
        # Nesting always means per-trace, whatever the schema says.
        length(value) == ntraces || throw(ArgumentError(
            "`$key` has $(length(value)) entries but the call draws $ntraces traces."))
        return value[i]
    end

    if value isa AbstractVector && !is_array_ok(trace_type, key)
        length(value) == ntraces || throw(ArgumentError(
            "`$key` has $(length(value)) entries but the call draws $ntraces traces. " *
            "Plotly does not read an array on `$key` as one value per point, so it can only " *
            "mean one per trace."))
        return value[i]
    end

    # Either a scalar, or an array Plotly reads per data point. Hand it over unchanged.
    return value
end

"""
    series_name(value, s::Series) -> String or nothing

`name` is the one attribute that rewrites a value rather than placing it. A scalar name across
several traces is suffixed with the component index, because three legend entries all reading
"position" is never what was meant.
"""
function series_name(value, s::Series)
    value === nothing && return nothing
    value isa AbstractVector && return value          # already distributed, one per trace
    s.ncomponents <= 1 && return value
    return string(value, "[", s.component, "]")
end
