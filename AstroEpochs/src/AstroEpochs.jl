# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

__precompile__()
"""
Module containing time system implementations for astronomical times.
Supports high-precision time representations using dual-float Julian Date
storage and conversions between scales such as TT, TAI, UTC, TDB, TCB, TCG.

The API is inspired by Astropy's `Time`, and so are the algorithms: the scale transforms, TDB − TT,
the calendar and the ISOT formatting are ERFA's, the library Astropy computes time with, and a
two-part date is normalised as Astropy's is. AstroEpochs agrees with Astropy to the last bit of
the date across 1600-2500 (test/test_correctness_astropy_benchmark.jl), except UTC before 1972,
below. See THIRD_PARTY_NOTICES.md for what is adapted from ERFA, Astropy and Tempo.jl.

Organization. Every algorithm here is one of two kinds, and each has its own folder:
- `scales/` changes which clock an instant is measured on: `offsets.jl` (TT, TAI, TCG, TDB, TCB),
  `tdb.jl` (TDB − TT), `leap_seconds.jl` (UTC) and `graph.jl` (which pairs are joined, and how a
  conversion is routed).
- `formats/` changes how an instant is written: `calendar.jl` (Gregorian ↔ Julian date, time of
  day, leap seconds) and `formats.jl` (JD, MJD, ISOT).
This file holds the types, the `Time` struct, its construction, accessors and arithmetic.

Notes
- Leap seconds come from the IERS leap-second list, downloaded on first use. A stored copy past
  its expiry is replaced the next time a session loads the table; a running session keeps the table
  it loaded. See `scales/leap_seconds.jl`. UTC on a day that ends with a leap second is ERFA's
  quasi-Julian date, whose fraction of the day runs over 86401 s, so 23:59:60.5 is representable.
- UTC before 1972 uses TAI − UTC = 0, with a warning once per session, and steps to 10 s at
  1972-01-01T00:00:00. Astropy (ERFA) applies the drifting pre-1972 offsets, which differ by up
  to about 10 s; that is a known gap, not modelled here yet. UTC after the leap-second list's
  expiry assumes no further leap seconds; more than five years after it, as ERFA does, that
  warns once per session.
- TDB − TT is the Fairhead & Bretagnon series (ERFA `dtdb`) at the geocentre, as Astropy uses for a
  Time with no location. ERFA documents it as agreeing with ephemeris-based TDB to a few nanoseconds
  over 1600-2200; AstroEpochs has not checked that against an ephemeris. `_tdb_minus_tt_parts` in
  `scales/tdb.jl` is the one function scale conversions call for it. An observer's terms come from
  `tdb_minus_tt(t, position)`, and are not part of a scale conversion.
- Scales and formats are tags, such as `TT()` and `JD()`. Symbols are only how a `Time` stores them.
"""
module AstroEpochs

using Printf

export Time
export TAI, TT, TDB, UTC, TCB, TCG
export JD, MJD, ISOT
export refresh_leap_seconds!
export isot, tai_minus_utc, tdb_minus_tt

import Base: +,-

include("constants.jl")

# Format transforms that need no Time: the calendar, which the leap-second table also uses.
include("formats/calendar.jl")

# Scale transforms: the offset functions, then the graph that joins them.
include("scales/tdb.jl")
include("scales/offsets.jl")
include("scales/leap_seconds.jl")
include("scales/graph.jl")

abstract type AbstractTimeScale end
abstract type AbstractTimeFormat end

"""
    TT()

Terrestrial Time scale.

# Example
```jldoctest
Time(2451545.0, TT(), JD()).scale

# output

TT()
```
"""
struct TT   <: AbstractTimeScale end

"""
    TDB()

Barycentric Dynamical Time scale.

# Example
```jldoctest
Time(2451545.0, TDB(), JD()).scale

# output

TDB()
```
"""
struct TDB  <: AbstractTimeScale end

"""
    UTC()

Coordinated Universal Time scale.

# Example
```jldoctest
Time(2451545.0, UTC(), JD()).scale

# output

UTC()
```
"""
struct UTC  <: AbstractTimeScale end

"""
    TCB()

Barycentric Coordinate Time scale.

# Example
```jldoctest
Time(2451545.0, TCB(), JD()).scale

# output

TCB()
```
"""
struct TCB  <: AbstractTimeScale end

"""
    TCG()

Geocentric Coordinate Time scale.

# Example
```jldoctest
Time(2451545.0, TCG(), JD()).scale

# output

TCG()
```
"""
struct TCG  <: AbstractTimeScale end

"""
    TAI()

International Atomic Time scale.

# Example
```jldoctest
Time(2451545.0, TAI(), JD()).scale

# output

TAI()
```
"""
struct TAI  <: AbstractTimeScale end

"""
    JD()

Julian Date format.

# Example
```jldoctest
Time(2451545.0, TT(), JD()).format

# output

JD()
```
"""
struct JD     <: AbstractTimeFormat end

"""
    MJD()

Modified Julian Date format.

# Example
```jldoctest
Time(51544.5, TT(), MJD()).format

# output

MJD()
```
"""
struct MJD    <: AbstractTimeFormat end   

"""
    ISOT()

ISO 8601 Time format.

# Example
```jldoctest
Time("2000-01-01T12:00:00.000", TT(), ISOT()).format

# output

ISOT()
```
"""
struct ISOT   <: AbstractTimeFormat end

# Tag -> the symbol a Time stores
@inline _scale_symbol(::TT)  = :tt       # COV_EXCL_LINE
@inline _scale_symbol(::TDB) = :tdb      # COV_EXCL_LINE
@inline _scale_symbol(::UTC) = :utc      # COV_EXCL_LINE
@inline _scale_symbol(::TCB) = :tcb      # COV_EXCL_LINE
@inline _scale_symbol(::TCG) = :tcg      # COV_EXCL_LINE
@inline _scale_symbol(::TAI) = :tai      # COV_EXCL_LINE

@inline _format_symbol(::JD)     = :jd    # COV_EXCL_LINE
@inline _format_symbol(::MJD)    = :mjd   # COV_EXCL_LINE
@inline _format_symbol(::ISOT)   = :isot  # COV_EXCL_LINE

# The symbol a Time stores -> its tag
function _scale_tag(s::Symbol)
    s === :tt  && return TT()
    s === :tdb && return TDB()
    s === :utc && return UTC()
    s === :tcb && return TCB()
    s === :tcg && return TCG()
    s === :tai && return TAI()
    throw(ArgumentError("Time: unknown time scale :$s"))
end

function _format_tag(f::Symbol)
    f === :jd   && return JD()
    f === :mjd  && return MJD()
    f === :isot && return ISOT()
    throw(ArgumentError("Time: unknown time format :$f"))
end

""" 
    _typename_paren(x) = string(nameof(typeof(x)), "()")

Print a type name with parentheses for display.
"""
_typename_paren(x) = string(nameof(typeof(x)), "()")


# Marks the constructor that takes Julian-date parts whatever the format; see `_time_jd`.
struct _JDParts end

"""
    Time(value, scale, format)
    Time(part1, part2, scale, format)
    Time(string, scale, ISOT())
    Time(year, month, day, hour, minute, second, scale, format)

An astronomical epoch: a Julian date in two parts, on one time scale, with a format that sets how
it is written.

The scale says how to read what you pass in. It does not convert: the same string read as `UTC()`
and as `TT()` gives two different instants. To convert a time you hold, use `t.tt` or `Time(t, TT())`.

# Properties
- `t.jd1`, `t.jd2`: The two parts of the Julian date, a whole day and a fraction between -0.5 and 0.5.
- `t.scale`, `t.format`: The tags, such as `TT()` and `JD()`.
- `t.jd`, `t.mjd`: The date as one number. One number holds less precision than the two parts.
- `t.isot`: The date as a string, to the millisecond. See [`isot`](@ref) for other precisions.
- `t.tt`, `t.tai`, `t.utc`, `t.tdb`, `t.tcb`, `t.tcg`: The same instant in that scale, as a `Time`.

# Notes
- Numbers are Julian dates for `JD()` and Modified Julian Dates for `MJD()`. A whole number is
  promoted to floating point. Two parts are summed, which holds more precision than one number.
- A string has the form `YYYY-MM-DDTHH:MM:SS`, with optional decimal seconds and a four-digit year.
  The time of day is required, and there is no time zone suffix, since the scale says what clock it is.
- A second of 60 or more is refused, except in UTC in the last minute of a day that ends with a leap
  second.
- Conversion promotes to at least `Float64`, so a `Time{Float32}` comes back as `Time{Float64}`.
  `BigFloat` and derivative types carry through.
- A date that is not finite passes through arithmetic and conversion, and has no calendar date to write.

# Examples
```jldoctest
t = Time(2451545.0, 0.25, TT(), JD())

println(t.isot)
println(t.mjd)
println(t.scale)

# output

2000-01-01T18:00:00.000
51544.75
TT()
```
"""
struct Time{T<:Real}
    _jd1::T
    _jd2::T
    _scale::Symbol
    _format::Symbol

    # Julian-date parts, whatever the format; the format only sets how the time displays.
    # Reached through `_time_jd`, never directly by a caller.
    function Time{T}(::_JDParts, jd1::T, jd2::T, scale::Symbol, format::Symbol) where {T<:Real}
        jd1, jd2 = _rebalance(jd1, jd2)
        _validate_scale(scale)
        _validate_format(format)
        new{T}(jd1, jd2, scale, format)
    end
end


"""
    _time_jd(jd1, jd2, scale::Symbol, format::Symbol) -> Time

A `Time` from the two parts of a Julian date, displayed in `format`. The internal way to copy or
shift a time while keeping its format; the public two-part constructor reads its parts in the
format it names.
"""
function _time_jd(jd1::Real, jd2::Real, scale::Symbol, format::Symbol)
    T = promote_type(typeof(jd1), typeof(jd2))
    return Time{T}(_JDParts(), T(jd1), T(jd2), scale, format)
end

# The stored scale and format, for this package's own code. `t.scale` and `t.format` return tags.
@inline _scale_of(t::Time)  = getfield(t, :_scale)
@inline _format_of(t::Time) = getfield(t, :_format)

# A date cannot be held in integers, since its fraction of a day is what it carries.
_real_parts(val1::Real, val2::Real) = promote(val1, val2)
_real_parts(val1::Integer, val2::Integer) = (float(val1), float(val2))

# Two parts of a date in `format`, summed: Julian-date parts for `:jd`, Modified Julian Date parts
# for `:mjd`, as in AstroPy.
function _time_parts(val1::Real, val2::Real, scale::Symbol, format::Symbol)
    val1, val2 = _real_parts(val1, val2)
    format === :isot && throw(ArgumentError(
        "Time: an ISOT time is a string, such as Time(\"2024-01-01T00:00:00\", UTC(), ISOT()); " *
        "two numbers are the parts of a Julian date or a Modified Julian Date."))
    if format === :mjd
        # Normalise first, then add the zero point to the whole day, as Astropy's TimeMJD does.
        # Adding it to a fractional first part rounds the fraction at the magnitude of a Julian
        # date, which cost up to 15 µs.
        day, frac = _rebalance(val1, val2)
        return _time_jd(day + MJD_EPOCH, frac, scale, :mjd)
    end
    return _time_jd(val1, val2, scale, format)
end

# Anything the methods below do not take: name what was wrong, and show the form that works.
function _invalid_time_input(pieces...)
    msg = IOBuffer()
    for (name, value, want, wanted) in pieces
        value isa want || print(msg, name, " must be ", wanted, "; got ", typeof(value), ". ")
    end
    throw(ArgumentError("Time: invalid time input types. " * String(take!(msg)) *
                        "For example: Time(2451545.0, TT(), JD())."))
end

Time(val1, val2, scale, format) = _invalid_time_input(
    ("the first part", val1, Real, "Real"), ("the second part", val2, Real, "Real"),
    ("scale", scale, AbstractTimeScale, "a time scale tag such as TT()"),
    ("format", format, AbstractTimeFormat, "a time format tag such as JD()"))

Time(value, scale, format) = _invalid_time_input(
    ("value", value, Union{Real, AbstractString}, "Real, or a string for ISOT()"),
    ("scale", scale, AbstractTimeScale, "a time scale tag such as TT()"),
    ("format", format, AbstractTimeFormat, "a time format tag such as JD()"))

function _time_value(value::Real, scale::Symbol, format::Symbol)
    _validate_inputcoupling(value,format)
    return _from_format(value, scale, format)
end

function _time_value(isostr::AbstractString, scale::Symbol, format::Symbol)
    _validate_inputcoupling(isostr,format)

    y, m, d, h, mi, s = _isot_to_date(isostr)
    jd1, jd2 = dtf2d(scale, y, m, d, h, mi, s)          # ERFA, as Astropy parses it

    jd1, jd2 = _rebalance(jd1, jd2)
    return _time_jd(jd1, jd2, scale, format)
end

"""
    Time(part1, part2, s, f) -> Time

A time from two parts of a date: Julian-date parts for `JD()`, Modified Julian Date parts for
`MJD()`. See [`Time`](@ref AstroEpochs.Time).
"""
Time(jd1::Real, jd2::Real, s::AbstractTimeScale, f::AbstractTimeFormat) =
    _time_parts(jd1, jd2, _scale_symbol(s), _format_symbol(f))


"""
    Time(value, s, f) -> Time

A time from one number: a Julian date for `JD()`, a Modified Julian Date for `MJD()`. See
[`Time`](@ref AstroEpochs.Time).
"""
Time(value::Real, s::AbstractTimeScale, f::AbstractTimeFormat) =
    _time_value(value, _scale_symbol(s), _format_symbol(f))

# The tag as it prints, for messages. A symbol that names no scale or format prints as itself.
@inline _scale_tag_str(s::Symbol) =
    s in TIME_SCALES ? _typename_paren(_scale_tag(s)) : ":" * String(s)

@inline _format_tag_str(f::Symbol) =
    f in TIME_FORMATS ? _typename_paren(_format_tag(f)) : ":" * String(f)

"""
    Time(string, s, ISOT()) -> Time

A time from a string of the form `YYYY-MM-DDTHH:MM:SS`, with optional decimal seconds. The time of day
is required and there is no time zone suffix; the scale comes from `s`. See
[`Time`](@ref AstroEpochs.Time).
"""
Time(isostr::AbstractString, s::AbstractTimeScale, f::AbstractTimeFormat) =
    _time_value(isostr, _scale_symbol(s), _format_symbol(f))

"""
    Time(year, month, day, hour, minute, second, s::AbstractTimeScale, f::AbstractTimeFormat) -> Time

Construct a time from calendar fields in the scale `s`, displayed in the format `f`.

# Arguments
- `year`, `month`, `day`, `hour`, `minute`: Integers.
- `second`: Real, and it can carry a fraction. In UTC, 60 and above is valid only in the last minute
  of a day that ends with a leap second.
- `s`: The time scale of the fields, such as `UTC()`.
- `f`: The format the time displays in, such as `ISOT()`.

# Returns
A `Time` at that calendar date.

# Examples
```julia
t = Time(2016, 12, 31, 23, 59, 60.5, UTC(), ISOT())
t.isot     # "2016-12-31T23:59:60.500"
```
"""
function Time(year::Integer, month::Integer, day::Integer, hour::Integer, minute::Integer,
              second::Real, s::AbstractTimeScale, f::AbstractTimeFormat)
    scale = _scale_symbol(s)
    _check_calendar_fields("$year-$month-$day", year, month, day, hour, minute)
    jd1, jd2 = _rebalance(dtf2d(scale, year, month, day, hour, minute, second)...)
    return _time_jd(jd1, jd2, scale, _format_symbol(f))
end

"""
    Time(t::Time, s::AbstractTimeScale) -> Time

The same instant as `t` in the scale `s`, keeping the format. `t.tdb` says the same thing when the
scale is known as you write the script; this form is for a scale held in a variable.

# Arguments
- `t`: The time to convert.
- `s`: The time scale to convert to, such as `TDB()`.

# Returns
A `Time` in the scale `s`.

# Examples
```jldoctest
a = Time(2451545.0, TT(), JD())
b = Time(2451545.5, TDB(), JD())

println(Time(b, a.scale).scale)

# output

TT()
```
"""
Time(t::Time, s::AbstractTimeScale) = getproperty(t, _scale_symbol(s))

"""
    Time(t::Time, f::AbstractTimeFormat) -> Time

`t` displayed in the format `f`. The format only sets how a time displays, so the date is held
exactly as it was, which rebuilding the time from `t.mjd` or `t.isot` does not do.

# Arguments
- `t`: The time.
- `f`: The format to display in, such as `ISOT()`.

# Returns
A `Time` at the same date, in the same scale.

# Examples
```jldoctest
t = Time(2451545.25, TT(), JD())

println(Time(t, ISOT()).format)

# output

ISOT()
```
"""
Time(t::Time, f::AbstractTimeFormat) = _time_jd(t.jd1, t.jd2, _scale_of(t), _format_symbol(f))

"""
    Time{T}(t::Time) -> Time{T}

`t` with its date held as the numeric type `T`, in the same scale and format.

# Arguments
- `t`: The time.

# Returns
A `Time{T}`. `convert(Time{T}, t)` does the same.

# Examples
```jldoctest
t = Time(2451545.0, 0.25, TT(), JD())

println(typeof(Time{BigFloat}(t)))

# output

Time{BigFloat}
```
"""
Time{T}(t::Time) where {T<:Real} = _time_jd(T(t.jd1), T(t.jd2), _scale_of(t), _format_of(t))

Base.convert(::Type{Time{T}}, t::Time) where {T<:Real} = Time{T}(t)
Base.convert(::Type{Time{T}}, t::Time{T}) where {T<:Real} = t

"""
    tai_minus_utc(t::Time) -> Float64

TAI − UTC at `t`, in seconds, from the leap-second list.

# Arguments
- `t`: The time, in any scale.

# Notes
- Through a leap second the old value holds, until the following UTC midnight.
- After the list's last entry this is the last value, since future leap seconds are unknown.
- Before 1972 the list has no value, so this returns 0 and warns.
- A date that is not finite returns NaN.

# Returns
TAI − UTC in seconds.

# Examples
```julia
tai_minus_utc(Time("2024-01-01T00:00:00", UTC(), ISOT()))     # 37.0
```
"""
function tai_minus_utc(t::Time)
    u = t.utc
    # The UTC midnight that starts the day, from the two parts: the table changes only at midnight,
    # and one number at the size of a Julian date would round the last moments of a day into the next.
    midnight = (u.jd1 - J2000_EPOCH) + floor(u.jd2 + 0.5) - 0.5
    isfinite(midnight) || return NaN
    return _tai_minus_utc_days(midnight)
end

"""
    tdb_minus_tt(t::Time) -> seconds
    tdb_minus_tt(t::Time, position) -> seconds

TDB − TT at `t`, in seconds, from ERFA's `dtdb`: at the geocentre, or at an observer on the Earth.

# Arguments
- `t`: The time, in any scale.
- `position`: The observer's position from the Earth's centre, as three numbers in km, in Earth-fixed
  axes.

# Notes
- This returns a number. It does not change `t`, and `t.tdb` stays at the geocentre whatever you pass.
- The observer's terms take UT1, and UTC stands in for it. For a station on the Earth that changes
  the result by at most about 0.13 ns.

# Returns
TDB − TT in seconds. The calculation promotes to at least `Float64`.

# Examples
```julia
t = Time("2024-01-01T00:00:00", TT(), ISOT())
tdb_minus_tt(t)
tdb_minus_tt(t, [-2353.6, -4641.3, 3677.1])
```
"""
function tdb_minus_tt(t::Time)
    tt = t.tt
    return _tdb_minus_tt_parts(tt.jd1, tt.jd2)
end

function tdb_minus_tt(t::Time, position::AbstractVector)
    length(position) == 3 || throw(ArgumentError(
        "tdb_minus_tt: position is three numbers in km, in Earth-fixed axes; got $(length(position))."))
    tt, ut = t.tt, t.utc
    day = (ut.jd1 - 0.5) + ut.jd2
    fraction_of_day = day - floor(day)      # `mod` of a dual number has no derivative at a whole day
    return dtdb(tt.jd1, tt.jd2, fraction_of_day, atan(position[2], position[1]),
                hypot(position[1], position[2]), position[3])
end

_fixed_jd(jd::Union{Float16, Float32, Float64}) = @sprintf("%.9f", jd)
_fixed_jd(jd) = string(jd)

# Printing a time never fails, even one with no calendar date: an optimizer can step to a NaN time
# and the message about what went wrong still has to print.
function _isot_for_display(t::Time)
    try
        return t.isot
    catch e
        e isa InterruptException && rethrow()
        return string("no calendar date (Julian date ", t.jd, ")")
    end
end

"""
    Base.show(io::IO, ::MIME"text/plain", t::Time)

The full display at the REPL: the value in the time's format, its scale and its format. A Julian date
is written in fixed notation, to the nine decimal places one number holds.
"""
function Base.show(io::IO, ::MIME"text/plain", t::Time)
    format = _format_of(t)
    value = format === :jd ? _fixed_jd(t.jd) : format === :mjd ? t.mjd : _isot_for_display(t)

    println(io, "AstroEpochs.Time")
    println(io, "  value  = ", value)
    println(io, "  scale  = ", _scale_tag_str(_scale_of(t)))
    println(io, "  format = ", _format_tag_str(format))
end

"""
    Base.show(io::IO, t::Time)

One line, for a time printed inside something else: a vector, a tuple, a string. It reads as the call
that makes the time, and for `JD()` and `MJD()` as closely as one or two numbers can say it.
"""
function Base.show(io::IO, t::Time)
    format, scale = _format_of(t), _scale_tag_str(_scale_of(t))
    if format === :isot && isfinite(t.jd1 + t.jd2)
        text = _isot_for_display(t)
        startswith(text, "no calendar") || return print(io, "Time(\"", text, "\", ", scale, ", ISOT())")
    end
    format === :mjd && return print(io, "Time(", t.mjd, ", ", scale, ", MJD())")
    return print(io, "Time(", t.jd1, ", ", t.jd2, ", ", scale, ", JD())")
end

"""
    _validate_format(format::Symbol)

Validate time format against supported formats
"""
function _validate_format(format::Symbol)
    if format ∉ TIME_FORMATS
        supported = join(map(_format_tag_str, collect(TIME_FORMATS)), ", ")
        throw(ArgumentError("Time: invalid time format $( _format_tag_str(format) ). Supported: [$supported]"))
    end
end

"""
    function _validate_scale(scale::Symbol)

Validate time scale against supported scales
"""        
function _validate_scale(scale::Symbol)
    if scale ∉ TIME_SCALES
        supported = join(map(_scale_tag_str, collect(TIME_SCALES)), ", ")
        throw(ArgumentError("Time: invalid time scale $( _scale_tag_str(scale) ). Supported: [$supported]"))
    end
end

"""
    _validate_inputcoupling(value::Any, format::Symbol)

Validate that `value` is compatible with the specified time `format`.
"""
function _validate_inputcoupling(value::Any, format::Symbol)
    if format === :isot
        if !(value isa AbstractString)
            throw(ArgumentError("Time: format $( _format_tag_str(format) ) requires AbstractString input. Got $(typeof(value))"))
        end
    elseif format === :jd || format === :mjd
        if !(value isa Real)
            throw(ArgumentError("Time: format $( _format_tag_str(format) ) requires Real input. Got $(typeof(value))"))
        end
    else
        throw(ArgumentError("Time: unsupported format $( _format_tag_str(format) ) in this constructor."))
    end
end

"""
    Base.getproperty(t::Time, name::Symbol)

Provide field access and on-demand conversions for Time via property syntax.

Arguments
- t::Time — the epoch object.
- name::Symbol — one of:
  - The date's parts: :jd1, :jd2
  - The scale and format tags: :scale, :format
  - Format accessors: :jd (Real), :mjd (Real), :isot (String)
  - Scale conversions: any of `:tt, :tai, :tdb, :utc, :tcg, :tcb` to return a new Time in that scale.

Returns
- Real for :jd, :mjd, :jd1, :jd2
- A tag for :scale and :format, such as `TT()` and `JD()`
- String for :isot
- Time for scale symbols (converted, preserving current `format`)

Notes
- Scale conversions are performed through apply_transforms following 
  the configured conversion graph.
- When converting scales, `jd1/jd2` are rebalanced to keep invariants.
- Unknown properties throw an error.
- Scale conversions preserve format and return a new Time object;
  the original is unchanged.
- Format computations return a time value (as opposed to a new Time struct).

# Examples
```julia
using AstroEpochs

t = Time("2024-02-29T12:34:56.123", TAI(), ISOT())

# The date's parts, and the scale and format
t.jd1; t.jd2
t.scale == TAI()      # true
t.format == ISOT()    # true

# Format accessors
t.jd        # numeric JD
t.mjd       # numeric MJD
t.isot      # ISO 8601 string

# Scale conversion (preserves the current format, ISOT)
tt_time = t.tt
tai_time = tt_time.tai
```
"""
function Base.getproperty(t::Time, name::Symbol)

    # Julian date elements jd1 and jd2
    if name === :jd1
        return getfield(t, :_jd1)
    elseif name === :jd2
        return getfield(t, :_jd2)
    elseif name === :scale
        return _scale_tag(_scale_of(t))
    elseif name === :format
        return _format_tag(_format_of(t))
    end

    # Time format conversion
    if name === :isot
        return _to_isot(t)
    elseif name === :jd
        return getfield(t, :_jd1) + getfield(t, :_jd2)
    elseif name === :mjd
        return _to_mjd(t)
    end

    # Time scale conversion 
    if name in TIME_SCALES
        newjd1, newjd2 = apply_transforms(t.jd1, t.jd2, _scale_of(t), name)
        return _time_jd(newjd1, newjd2, name, _format_of(t))
    end

    # Fall back to normal field access
    error("Unknown property `$name` for Time. Available: " * join(propertynames(t), ", ") * ".")

end


"""
    _two_sum(a, b) -> (sum, error)

`a + b` exactly, as the rounded sum and its rounding error (Shewchuk 1997), as Astropy's
`astropy.time.utils.two_sum`.
"""
@inline function _two_sum(a, b)
    x  = a + b
    eb = x - a
    ea = x - eb
    eb = b - eb
    ea = a - ea
    return x, ea + eb
end

"""
    _rebalance(jd1, jd2) -> (day, fraction)

The two-part date `jd1 + jd2` as a whole number of days and a fraction in [-0.5, 0.5], the sum
formed exactly. This is Astropy's `astropy.time.utils.day_frac` (BSD 3-Clause, Copyright (c)
2011-2026 Astropy Developers; see THIRD_PARTY_NOTICES.md), the normalisation every Astropy Time
holds its date in, so a Time here splits its date as Astropy's does. Halves round to even, as
NumPy's `round` does and Julia's does by default.
"""
function _rebalance(jd1::T, jd2::T) where {T<:Real}
    sum12, err12 = _two_sum(jd1, jd2)
    day = round(sum12)
    # The remainder can land just outside ±0.5, costing a bit; correct for that, taking care at
    # exactly ±0.5 where the rounding error decides which way the day should have gone.
    frac, check = _two_sum(sum12 - day, err12)
    excess = frac * sign(check) != 0.5 ? round(frac) : round(frac + 2check)
    day += excess
    frac = sum12 - day
    frac += err12
    return T(day), T(frac)
end

"""
    _rebalance(a::Real, b::Real)

Promote to a common type and normalise, as the typed method.
"""
@inline _rebalance(a::Real, b::Real) = _rebalance(promote(a, b)...)

# Comparing two times means measuring them on one clock, so times on different clocks are refused
# with the conversion that makes them comparable.
function _require_one_scale(a::Time, b::Time, action::AbstractString, form::AbstractString)
    _scale_of(a) === _scale_of(b) || throw(ArgumentError(
        "Time: cannot $action times in different scales ($(_scale_tag_str(_scale_of(a))) and " *
        "$(_scale_tag_str(_scale_of(b)))); convert one first, for example `$form.$(_scale_of(a))`."))
    return nothing
end

"""
    -(t2::Time, t1::Time) -> days

The difference between two times in the same scale, in days of that scale. Their formats may differ.
Times in different scales throw; convert one first.

In UTC this is the difference between the scale's own dates, not elapsed time. For elapsed time, use
`t2.tai - t1.tai`.

# Examples
```jldoctest
a = Time(2451545.0, TT(), JD())

println((a + 0.25) - a)

# output

0.25
```
"""
function -(t2::Time, t1::Time)::Real
    _require_one_scale(t2, t1, "subtract", "t2 - t1")
    return (t2.jd1 - t1.jd1) + (t2.jd2 - t1.jd2)
end

"""
    -(t::Time, dt::Real) -> Time

`t`, `dt` days earlier, in the same scale and format. See [`+`](@ref +(::Time, ::Real)).
"""
Base.:-(t::Time, dt::Real) = t + (-dt)
Base.:-(t::Time, dt::Unsigned) = t + (-float(dt))     # negating an unsigned number wraps

"""
    +(t::Time, dt::Real) -> Time

`t`, `dt` days later, in the same scale and format. The numeric type promotes, so `dt` may carry
derivatives, and the operands may come in either order.

In UTC this moves along the scale's own dates. On a day that ends with a leap second a date's fraction
runs over 86401 s, so a result can differ from clock time by up to a second. To move by elapsed time,
use `(t.tai + dt).utc`.

# Examples
```jldoctest
t = Time(2451545.0, TT(), JD())

println((t + 0.25).isot)

# output

2000-01-01T18:00:00.000
```
"""
function +(t::Time, dt::Real)
    # Keep AD types by pairing dt with zero(dt)
    δ1, δ2 = _rebalance(dt, zero(dt))
    jd1 = t.jd1 + δ1
    jd2 = t.jd2 + δ2
    return _time_jd(jd1, jd2, _scale_of(t), _format_of(t))
end

function +(dt::Real, t::Time)
    return t + dt
end

"""
    ==(a::Time, b::Time)

Whether two times are in the same scale and hold the same date, to the last bit. The format plays no
part, and neither does the numeric type. `hash` agrees.

Two times written separately are rarely equal to the last bit, even when they mean the same instant:
`Time(2460000.123, TT(), JD())` and `Time(59999.623, TT(), MJD())` differ by 12 µs, because each is
one number rounded on its own. Use [`isapprox`](@ref isapprox(::Time, ::Time)) for those. Times in
different scales are never equal; convert one first.
"""
function Base.:(==)(a::Time, b::Time)
    _scale_of(a) === _scale_of(b) || return false
    return getfield(a, :_jd1) == getfield(b, :_jd1) && getfield(a, :_jd2) == getfield(b, :_jd2)
end

# Equal times hash alike, across numeric types. Adding zero turns -0.0, which equals 0.0 but hashes
# differently, into 0.0.
Base.hash(t::Time, h::UInt) =
    hash(t.jd1 + zero(t.jd1), hash(t.jd2 + zero(t.jd2), hash(_scale_of(t), hash(:Time, h))))

# Order. `<`, `sort`, `min` and `max` come from `isless`.
function Base.isless(a::Time, b::Time)
    _require_one_scale(a, b, "compare", "a < b")
    return (a.jd1 - b.jd1) + (a.jd2 - b.jd2) < 0
end

"""
    isapprox(a::Time, b::Time; atol = 1e-9) -> Bool

Whether two times in the same scale are within `atol` seconds of each other, a nanosecond by default.
Times written separately rarely satisfy `==`, since each number is rounded on its own; this is the
comparison for them. Times in different scales throw; convert one first.

# Arguments
- `a`, `b`: The times, in the same scale.
- `atol`: The tolerance, in seconds.

# Returns
`true` when the two differ by no more than `atol` seconds.

# Examples
```jldoctest
a = Time(2460000.123, TT(), JD())
b = Time(59999.623, TT(), MJD())

println(a == b)
println(isapprox(a, b; atol = 1e-3))

# output

false
true
```
"""
function Base.isapprox(a::Time, b::Time; atol::Real = 1e-9)
    _require_one_scale(a, b, "compare", "isapprox(a, b")
    return abs(a - b) * 86400.0 <= atol
end

# A time is one thing in a broadcast, so `t .+ offsets` gives one time per offset.
Base.broadcastable(t::Time) = Ref(t)

# What `t.` offers: the properties that work, not the fields that hold them.
Base.propertynames(::Time, private::Bool = false) =
    private ? (:jd1, :jd2, :scale, :format, :jd, :mjd, :isot, :tt, :tai, :tdb, :utc, :tcg, :tcb,
               :_jd1, :_jd2, :_scale, :_format) :
              (:jd1, :jd2, :scale, :format, :jd, :mjd, :isot, :tt, :tai, :tdb, :utc, :tcg, :tcb)

# Format transforms, which read and write a Time, so they come after it.
include("formats/formats.jl")

end
