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

Return time value in Modified Julian Date (MJD) format.
"""
@inline function _to_mjd(t::Time)
    return (t.jd1 - MJD_EPOCH) + t.jd2
end

"""
    isot(t::Time; digits = 3) -> String

`t` as an ISO 8601 string with `digits` decimal places of seconds. `t.isot` is this at the default
three places.

The string is written by ERFA's `d2dtf`, as Astropy formats it: the rounding carries into the minute,
hour and date, and in UTC a time inside a leap second reads 23:59:60.

# Arguments
- `t`: The time, written in its own scale.
- `digits`: Decimal places of seconds. With 0 the string ends at the whole second.

# Returns
The string `YYYY-MM-DDThh:mm:ss`, with `digits` decimal places after a point when `digits` is above 0.
It names no scale.

# Notes
- `digits` runs from 0 to 12, and anything else throws. A `Time` resolves about 10 ps, so more than 11
  places carry nothing.
- A date that is not finite, or too far off for the calendar, has no string and throws.

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

Construct a Time{T} from a JD/MJD numeric value, preserving T = typeof(value). An integer becomes
floating point, since a date in integers has no fraction of a day to carry.
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

Validate and parse an ISO 8601 string of the form `"YYYY-MM-DDTHH:MM:SS.sss"`.
Returning calendar fields: `(year, month, day, hour, minute, second)`.
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

