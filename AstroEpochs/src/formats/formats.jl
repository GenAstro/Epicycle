# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# ──────────────────────────────── formats ──────────────────────────────────
#
# Format transforms: JD, MJD and ISOT to and from the split Julian date a Time holds. A format
# never changes the instant, only how it is written; the calendar arithmetic under ISOT is in
# calendar.jl.
#
# Adding a format means: its symbol in TIME_FORMATS, its tag type in AstroEpochs.jl with the
# others, and a reader and a writer here, dispatched from `_from_format` and `getproperty`.

const TIME_FORMATS = Set([:jd, :mjd, :isot])

"""
    _to_mjd(t::Time)

Return the date as a single Modified Julian Date value.
"""
@inline function _to_mjd(t::Time)
    return (t.jd1 - MJD_EPOCH) + t.jd2
end

"""
    isot(t::Time; digits = 3) -> String

Format `t` as a calendar timestamp in its current scale, with `digits` decimal places in the
seconds field. `t.isot` is equivalent to `isot(t; digits = 3)`.

The output has the form `YYYY-MM-DDTHH:MM:SS`, followed by fractional seconds when `digits > 0`.
It does not include a scale or time zone suffix. Rounding carries into the next minute, hour,
or date as needed. UTC output can include `23:59:60` during a positive leap second.

# Notes
- `digits` must be between 0 and 12; other values throw an error.
- For `Time{Float64}`, fractional-day resolution is on the order of 10 ps. Extra displayed digits
  do not guarantee extra accuracy; precision depends on the input and numeric type.
- Nonfinite dates and dates outside the supported calendar range throw an error.

# Returns
The timestamp as a `String`.

# Examples
```jldoctest
t = Time(2451545.0, 0.25, TT(), JD())

println(isot(t; digits = 6))

# output

2000-01-01T18:00:00.000000
```
"""
function isot(t::Time; digits::Integer = 3)
    0 <= digits <= 12 || throw(ArgumentError("isot: digits runs from 0 to 12; got $digits."))
    isfinite(t.jd1 + t.jd2) || throw(ArgumentError(
        "Time: a time whose Julian date is $(t.jd1 + t.jd2) has no calendar date to write."))
    y, m, d, h, mi, s, f = d2dtf(_scale_of(t), digits, t.jd1, t.jd2)
    whole = @sprintf("%04d-%02d-%02dT%02d:%02d:%02d", y, m, d, h, mi, s)
    return digits == 0 ? whole : whole * "." * lpad(string(f), digits, '0')
end

_to_isot(t::Time) = isot(t)

"""
    _from_format(value::Real, scale::Symbol, format::Symbol) -> Time{T}

Create a `Time{T}` from a JD or MJD value, preserving `T = typeof(value)`.
Integer input is promoted to floating point so the date can include a fractional day.
"""
_from_format(value::Integer, scale::Symbol, format::Symbol) = _from_format(float(value), scale, format)

function _from_format(value::Real, scale::Symbol, format::Symbol)
    T = typeof(value)
    if format == :jd
        jd1 = floor(T, value)
        jd2 = T(value - jd1)
        jd1, jd2 = _rebalance(jd1, jd2)
        return _time_jd(jd1, jd2, scale, :jd)
    elseif format == :mjd
        jd1 = floor(T, value) + T(MJD_EPOCH)
        jd2 = T(value - floor(T, value))
        jd1, jd2 = _rebalance(jd1, jd2)
        return _time_jd(jd1, jd2, scale, :mjd)
    else
        throw(ArgumentError("Internal Error: Unsupported time format: $format")) # COV_EXCL_LINE
    end
end

# The checks a string and calendar fields share. `date` is how the input reads in the message.
function _check_calendar_fields(date, year, month, day, hour, minute)
    1 <= month <= 12 || throw(ArgumentError("Time: the month must be between 1 and 12; got $month."))
    1 <= day <= _days_in_month(year, month) || throw(ArgumentError(
        "Time: $date is not a date; $(_MONTH_NAMES[month]) $year has $(_days_in_month(year, month)) days."))
    0 <= hour < 24 || throw(ArgumentError("Time: the hour must be between 0 and 23; got $hour."))
    0 <= minute < 60 || throw(ArgumentError("Time: the minute must be between 0 and 59; got $minute."))
    return nothing
end

"""
    _isot_to_date(isostr::AbstractString) → Tuple{Int, Int, Int, Int, Int, Real}

Parse and validate an ISOT timestamp, returning `(year, month, day, hour, minute, second)`.
"""
function _isot_to_date(isostr::AbstractString)
    # The form this package reads, with optional fractional seconds
    m = match(r"\A([0-9]{4})-([0-9]{2})-([0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2}(?:\.[0-9]+)?)\z", isostr)
    if m === nothing
        throw(ArgumentError("Time: \"$isostr\" is not in the form YYYY-MM-DDTHH:MM:SS with optional " *
                            "decimal seconds. The time of day is required, and the scale comes from " *
                            "the tag, so there is no time zone suffix."))
    end

    y  = parse(Int, m.captures[1])
    mth = parse(Int, m.captures[2])
    d  = parse(Int, m.captures[3])
    h  = parse(Int, m.captures[4])
    mi = parse(Int, m.captures[5])
    s  = parse(Float64, m.captures[6])

    _check_calendar_fields(isostr, y, mth, d, h, mi)

    return (y, mth, d, h, mi, s)
end

