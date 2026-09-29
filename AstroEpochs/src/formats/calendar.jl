# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT
#
# Adapted from Tempo.jl v1.3.1, src/convert.jl (MIT License).
# Copyright (c) 2022 Andrea Pasquale and Michele Ceresoli. See THIRD_PARTY_NOTICES.md.
# Tempo's routines follow the ERFA library (BSD 3-Clause), itself derived from the IAU SOFA
# library: `cal2jd` from ERFA cal2jd.c and `jd2cal` from ERFA jd2cal.c.
#
# Changes from Tempo: the leap-year test is AstroEpochs' own `_is_leap_year`, which Tempo's
# `isleapyear` duplicated; the arithmetic is unchanged, and test_correctness_tempo_parity.jl
# holds it to Tempo's own outputs.

# ─────────────────────────── Gregorian calendar ────────────────────────────
#
# Format transforms only: these relabel an instant, and never change it.

const _MONTH_NAMES = ("January", "February", "March", "April", "May", "June", "July",
                      "August", "September", "October", "November", "December")

const _MONTH_DAYS               = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
const _PREVIOUS_MONTH_END       = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
const _PREVIOUS_MONTH_END_LEAP  = (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

_is_leap_year(y) = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
_days_in_month(y, m) = m == 2 ? (_is_leap_year(y) ? 29 : 28) : (m in (4, 6, 9, 11) ? 30 : 31)

"Day of the year, 1 on January 1."
_day_in_year(month, day, isleap::Bool) =
    day + (isleap ? _PREVIOUS_MONTH_END_LEAP : _PREVIOUS_MONTH_END)[month]

"""
    hms2fd(hour, minute, second) -> fraction of a day

Hours, minutes and seconds to a day fraction. Throws `DomainError` outside 0-23 h, 0-59 min and
[0, 60) s.
"""
function hms2fd(h::Integer, m::Integer, s::Number)
    if h < 0 || h > 23
        throw(DomainError(h, "the hour shall be between 0 and 23."))
    elseif m < 0 || m > 59
        throw(DomainError(m, "the minutes must be between 0 and 59."))
    elseif s < 0 || s >= 60
        throw(DomainError(s, "the seconds must be between 0.0 and 59.99999999999999."))
    end
    return ((60 * (60 * h + m)) + s) / 86400
end

"""
    fd2hms(fd) -> (hour, minute, second)

A day fraction to hours, minutes and seconds. Throws `DomainError` outside [0, 1] day.
"""
function fd2hms(fd::Number)
    secinday = fd * 86400
    if secinday < 0 || secinday > 86400
        throw(DomainError(secinday,
            "seconds are out of range: they must be between 0 and 86400."))
    end
    hours = Int(secinday ÷ 3600)
    secinday -= 3600 * hours
    mins = Int(secinday ÷ 60)
    secinday -= 60 * mins
    return hours, mins, secinday
end

"""
    cal2jd(year, month, day) -> (2451545, days since J2000)

A Gregorian date to a two-part Julian date: J2000's Julian date, 2451545, and the whole number of
days from J2000 to noon on the date. The sum is the Julian date at noon on the date, so 0h is half
a day earlier. Years before 1583, months outside 1-12 and days outside the month throw
`DomainError`.

Reference: Seidelmann (1992), Explanatory Supplement to the Astronomical Almanac, §12.92.
"""
function cal2jd(Y::Integer, M::Integer, D::Integer)
    if Y < 1583
        throw(DomainError(Y, "the year shall be greater than 1583."))
    elseif M < 1 || M > 12
        throw(DomainError(M, "the month shall be between 1 and 12."))
    end

    isleap = _is_leap_year(Y)
    ly = (M == 2) && isleap
    if (D < 1) || (D > (_MONTH_DAYS[M] + ly))
        throw(DomainError(D, "the day shall be between 1 and $(_MONTH_DAYS[M] + ly)."))
    end

    Y = Y - 1
    d1 = 365 * Y + Y ÷ 4 - Y ÷ 100 + Y ÷ 400 - 730120     # J2000 day of the year's start
    d2 = _day_in_year(M, D, isleap)
    return 2451545, d1 + d2
end

"""
    jd2cal(dj1, dj2) -> (year, month, day, fraction of day)

A two-part Julian date to a Gregorian date and day fraction. The date may be split between the
parts in any way; the fraction is summed with compensation, so the split costs no precision.
Julian dates below -68569.5 (4713 BC January 1) or above 1e9 throw `DomainError`.

References: Seidelmann (1992), §12.92; Klein (2006), A Generalized Kahan-Babuska-Summation
Algorithm, Computing 76, 279-293, §3.
"""
function jd2cal(dj1::Number, dj2::Number)
    dj = dj1 + dj2
    if dj < -68569.5 || dj > 1e9
        throw(DomainError(dj, "the Julian Date shall be between -68569.5 and 1e9."))
    end

    d1 = round(Int, dj1)
    d2 = round(Int, dj2)
    jd = d1 + d2

    # Separate day and fraction, and form f1 + f2 + 0.5 by compensated summation.
    f1, f2 = promote(dj1 - d1, dj2 - d2)
    s  = 0.5 * one(f1)
    cs = zero(f1)
    for x in (f1, f2)
        t = s + x
        cs += abs(s) >= abs(x) ? (s - t) + x : (x - t) + s
        s = t
        if s >= 1
            jd += 1
            s -= 1
        end
    end
    f  = s + cs
    cs = f - s

    # A negative fraction borrows a day.
    if f < 0
        f = s + 1
        cs += (1 - f) + s
        s = f
        f = s + cs
        cs = f - s
        jd -= 1
    end

    # A fraction that rounds to one carries a day.
    if (f - 1) >= -eps(f) / 4
        t = s - 1
        cs += (s - t) - 1
        s = t
        f = s + cs
        if -eps() / 2 < f
            jd += 1
            f = max(f, zero(f))
        end
    end

    # The day number as a Gregorian date.
    l = jd + 68569
    n = (4l) ÷ 146097
    l -= (146097n + 3) ÷ 4
    i = (4000 * (l + 1)) ÷ 1461001
    l -= (1461i) ÷ 4 - 31
    k = (80l) ÷ 2447
    D = (l - (2447k) ÷ 80)
    l = k ÷ 11
    M = k + 2 - 12l
    Y = 100 * (n - 49) + i + l
    return Y, M, D, f
end
