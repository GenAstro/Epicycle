# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Problems found while reading an ephemeris file. A reader records every problem rather than
# stopping at the first, so one read reports everything that needs fixing.

"""
    Diagnostic

One problem found while reading an ephemeris file.

# Fields
- `severity`: `:fixed` for a departure from the format that the reader corrected, `:error` for
  one it could not, and `:note` for an assumption made reading a valid file, such as a
  caller's override or `ITRF2000` read as Epicycle's ITRF.
- `line`: The line number in the file, or 0 when the problem is not on one line.
- `field`: The keyword or part of the file involved, such as `"REF_FRAME"` or `"data"`.
- `message`: What was found, what was expected, and what was done or needs doing.

# Examples
```jldoctest
d = EpicycleIO.Diagnostic(:fixed, 3, "CREATION_DATE", "keyword in lower case; read as CREATION_DATE")
d.severity

# output

:fixed
```
"""
struct Diagnostic
    severity::Symbol
    line::Int
    field::String
    message::String
end

function Base.show(io::IO, d::Diagnostic)
    where_ = d.line > 0 ? rpad("line $(d.line)", 10) : rpad("", 10)
    print(io, where_, "  ", rpad(d.field, 14), "  ", d.message)
end

"""
    EphemerisReadError(path, problems)

The error thrown when an ephemeris file has problems the reader cannot correct. It lists every
problem found, not only the first.

# Fields
- `path`: The file read.
- `problems`: Every [`Diagnostic`](@ref) of severity `:error`, in file order.
- `fixed`: The problems that were corrected, shown after the errors for context.
- `hint`: A suggestion for how to proceed, or an empty string.

# Examples
```jldoctest
e = EpicycleIO.EphemerisReadError("a.oem",
        [EpicycleIO.Diagnostic(:error, 14, "REF_FRAME", "EME200 is not a known frame. Did you mean EME2000?")],
        EpicycleIO.Diagnostic[], "")
length(e.problems)

# output

1
```
"""
struct EphemerisReadError <: Exception
    path::String
    problems::Vector{Diagnostic}
    fixed::Vector{Diagnostic}
    hint::String
end

function Base.showerror(io::IO, e::EphemerisReadError)
    n = length(e.problems)
    print(io, "EphemerisReadError: ", basename(e.path), " has ", n,
          n == 1 ? " problem that needs fixing" : " problems that need fixing")
    for d in e.problems
        print(io, "\n  ", d)
    end
    if !isempty(e.fixed)
        print(io, "\n  and ", length(e.fixed), " corrected along the way:")
        for d in first(e.fixed, 5)
            print(io, "\n  ", d)
        end
        length(e.fixed) > 5 && print(io, "\n  …")
    end
    isempty(e.hint) || print(io, "\n  ", e.hint)
end

# Collects diagnostics during one read. `strict` turns every correction into an error.
mutable struct _Diagnostics
    path::String
    strict::Bool
    items::Vector{Diagnostic}
end
_Diagnostics(path::AbstractString; strict::Bool = false) = _Diagnostics(String(path), strict, Diagnostic[])

_fixed!(d::_Diagnostics, line, field, msg) =
    push!(d.items, Diagnostic(d.strict ? :error : :fixed, line, String(field), msg))
_error!(d::_Diagnostics, line, field, msg) = push!(d.items, Diagnostic(:error, line, String(field), msg))
_note!(d::_Diagnostics, line, field, msg) = push!(d.items, Diagnostic(:note, line, String(field), msg))
_has_errors(d::_Diagnostics) = any(x -> x.severity == :error, d.items)

# Throws if anything could not be corrected; otherwise warns once, listing the corrections, and
# returns the corrections and notes. Notes do not warn: they describe a valid file.
function _finish(d::_Diagnostics; hint::AbstractString = "")
    errors = filter(x -> x.severity == :error, d.items)
    fixed = filter(x -> x.severity == :fixed, d.items)
    isempty(errors) || throw(EphemerisReadError(d.path, errors, fixed, String(hint)))
    if !isempty(fixed)
        shown = join(("  " * sprint(show, x) for x in first(fixed, 10)), "\n")
        more = length(fixed) > 10 ? "\n  … and $(length(fixed) - 10) more" : ""
        @warn "$(basename(d.path)) read with $(length(fixed)) $(length(fixed) == 1 ? "fix" : "fixes"). " *
              "See diagnostics(eph).\n" * shown * more
    end
    return filter(x -> x.severity != :error, d.items)
end

# The closest of `candidates` to `value`, if it is close enough to be a likely typo.
function _suggest(value::AbstractString, candidates)
    v = uppercase(strip(value))
    best, dist = "", typemax(Int)
    for c in candidates
        k = _edit_distance(v, uppercase(c))
        k < dist && ((best, dist) = (c, k))
    end
    return dist <= max(2, length(v) ÷ 4) ? best : nothing
end

function _did_you_mean(value, candidates)
    s = _suggest(value, candidates)
    return s === nothing ? "" : " Did you mean $s?"
end

function _edit_distance(a::AbstractString, b::AbstractString)
    a, b = collect(a), collect(b)
    prev = collect(0:length(b))
    for i in 1:length(a)
        cur = similar(prev)
        cur[1] = i
        for j in 1:length(b)
            cur[j + 1] = min(prev[j + 1] + 1, cur[j] + 1, prev[j] + (a[i] == b[j] ? 0 : 1))
        end
        prev = cur
    end
    return prev[end]
end
