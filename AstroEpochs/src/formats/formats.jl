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
    _to_isot(t::Time)

Return time value in ISOT format rounding to millisecond.
"""
function _to_isot(t::Time)
    jd = t.jd1 + t.jd2
    y, m, d, fd = jd2cal(t.jd1,t.jd2)
    h, mi, s = fd2hms(fd)
    return @sprintf("%04d-%02d-%02dT%02d:%02d:%06.3f", y, m, d, h, mi, s)
end

"""
    _from_format(value::Real, scale::Symbol, format::Symbol) -> Time{T}

Construct a Time{T} from a JD/MJD numeric value, preserving T = typeof(value).
"""
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

"""
    _isot_to_date(isostr::String) → Tuple{Int, Int, Int, Int, Int, Real}

Validate and parse an ISO 8601 string of the form `"YYYY-MM-DDTHH:MM:SS.sss"`.
Returning calendar fields: `(year, month, day, hour, minute, second)`.
"""
function _isot_to_date(isostr::String)
    # Match ISO 8601 format with optional fractional seconds
    m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$", isostr)
    if m === nothing
        throw(ArgumentError("Time: Invalid ISO 8601 format: $isostr"))
    end

    y  = parse(Int, m.captures[1])
    mth = parse(Int, m.captures[2])
    d  = parse(Int, m.captures[3])
    h  = parse(Int, m.captures[4])
    mi = parse(Int, m.captures[5])
    s  = parse(Float64, m.captures[6])

    # Semantic range checks
    if !(1 <= mth <= 12)
        throw(ArgumentError("Month must be between 1 and 12. Got: $mth"))
    end
    if !(1 <= d <= _days_in_month(y, mth))
        throw(ArgumentError("Time: $isostr is not a date; " *
                            "$(_MONTH_NAMES[mth]) $y has $(_days_in_month(y, mth)) days."))
    end
    if !(0 <= h < 24)
        throw(ArgumentError("Hour must be between 0 and 23. Got: $h"))
    end
    if !(0 <= mi < 60)
        throw(ArgumentError("Minute must be between 0 and 59. Got: $mi"))
    end
    if !(0.0 <= s < 60.0)
        throw(ArgumentError("Seconds must be >= 0.0 and < 60.0. Got: $s"))
    end
    
    return (y, mth, d, h, mi, s) 
end

"""
    calhms2jd_prec(Y::I, M::I, D::I, h::I, m::I, sec::N) where {I <: Integer, N <: Number}

Convert calendar date and time to Julian Date
"""
function calhms2jd_prec(Y::I, M::I, D::I, h::I, m::I, sec::N) where {I <: Integer, N <: Number}
    j2000_epoch, daysfrom_j2000 = cal2jd(Y, M, D)
    frac_of_day = hms2fd(h, m, sec)

    return Float64(j2000_epoch + daysfrom_j2000), frac_of_day - 0.5
end
