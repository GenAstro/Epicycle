# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Moved here from AstroCallbacks, which is about quantities rather than output. It is gone
# from there — this is the only `report`.
#
# One thing changed on the way: `report` no longer knows what an epoch is. It used to
# special-case `Time` and write a bare Julian date, which says nothing about TT versus UTC and
# put a domain decision inside a text writer. The caller renders it now, and the error below
# says so. That removed this package's last dependency on any Astro package.

"""
    report(path; columns...)

Write named columns to a delimited text file.

The keywords become the column headings, so they are your words rather than whatever the
quantity happens to be called — someone who did not write the run has to read this. A
vector-valued column expands into components, so `position` becomes `position_1`, `position_2`,
`position_3`.

```julia
report("flight.txt"; time = [x.jd for x in t], position = r, velocity = v)
```

# Arguments
- `path`: the file to write. An existing file is overwritten.
- `delim`: what separates the columns. Two spaces by default.
- `columns...`: the data, one keyword per column. Every column must be the same length.

# Notes
Columns carry no units of their own — they are written as given, so name the keyword for what
the number is (`altitude_km`) if the reader will need to know.

Values must render on one line. An epoch does not, so convert it first — `t.jd`, `t.mjd` or
`t.isot` — and the file then records which scale and format you meant.

# Returns
The path written, so it can be passed straight on to something that reads it.

Throws an `ArgumentError` if no columns are given, if the columns differ in length, or if a
value's printed form spans lines and would break a row.
"""
function report(path::AbstractString; delim = "  ", columns...)
    names, cols = expand_columns(columns)
    n = length(first(cols))

    open(path, "w") do io
        println(io, join((rpad(nm, 24) for nm in names), delim))
        for i in 1:n
            println(io, join((rpad(_fmt(c[i]), 24) for c in cols), delim))
        end
    end
    return path
end

"""
    expand_columns(columns) -> (names, cols)

Turn named columns into flat, equal-length series ready to publish.

A vector-valued column becomes one series per component, named `base_1`, `base_2`, … Everything
else passes through under its own name.
"""
function expand_columns(columns)
    isempty(columns) && throw(ArgumentError(
        "report needs at least one column, named with a keyword: " *
        "report(\"flight.txt\"; time = t, position = r). The keywords become the column " *
        "headings."))

    names, cols = String[], Any[]
    for (name, col) in pairs(columns)
        if !isempty(col) && first(col) isa AbstractVector
            for j in eachindex(first(col))
                push!(names, "$(name)_$(j)")
                push!(cols, [row[j] for row in col])
            end
        else
            push!(names, String(name))
            push!(cols, col)
        end
    end

    n = length(first(cols))
    all(c -> length(c) == n, cols) || throw(ArgumentError(
        "columns have different lengths: " *
        join(("$(nm) $(length(c))" for (nm, c) in zip(names, cols)), ", ") *
        ". Columns from different subjects are different series and cannot share a report."))

    return names, cols
end

_fmt(x::Real) = string(round(float(x), sigdigits = 12))

# Guard the class rather than the case: any value whose printed form spans lines would turn one
# row into several, silently. An epoch is the one that turns up in practice, so the message
# names the way out.
function _fmt(x)
    s = string(x)
    if occursin('\n', s)
        throw(ArgumentError(
            "a $(typeof(x)) column cannot be written to a report: its printed form spans " *
            "several lines, which would break the row. If this is an epoch, render it first " *
            "— `t.jd`, `t.mjd` or `t.isot` — so the file records which scale and format you " *
            "meant. Otherwise convert it to a number or a string."))
    end
    return s
end
