# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT
#
# The Gregorian calendar and the time of day, ported from ERFA 2.0.1 (BSD 3-Clause; derived from
# the IAU SOFA library): `cal2jd` from src/cal2jd.c, `jd2cal` from src/jd2cal.c, `d2tf` from
# src/d2tf.c, and `d2dtf` and `dtf2d` from src/d2dtf.c and src/dtf2d.c. See THIRD_PARTY_NOTICES.md.
# These are the routines Astropy formats and parses dates with, so an ISOT string here is
# Astropy's; test_correctness_erfa_parity.jl holds them to pyerfa's own output.
#
# Two changes from the C: failures throw (the C returns a status), and a date is returned as a
# tuple rather than through pointers. Day numbers come from `round` and `trunc`, which return plain
# numbers for dual numbers, so under automatic differentiation the fraction of a day carries the
# derivative and the calendar date does not.

# ─────────────────────────── Gregorian calendar ────────────────────────────
#
# Format transforms only: these relabel an instant, and never change it.

const _MONTH_NAMES = ("January", "February", "March", "April", "May", "June", "July",
                      "August", "September", "October", "November", "December")

const _MONTH_DAYS = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

_is_leap_year(y) = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
_days_in_month(y, m) = m == 2 ? (_is_leap_year(y) ? 29 : 28) : (m in (4, 6, 9, 11) ? 30 : 31)

# ERFA_DNINT: nearest integer, halves away from zero. Julia's default `round` sends halves to the
# even neighbour, which differs from ERFA exactly at half-days.
_dnint(x) = round(x, RoundNearestTiesAway)

"""
    cal2jd(year, month, day) -> (2400000.5, MJD at 0h)

ERFA's `eraCal2jd`: a Gregorian date to a two-part Julian date for 0h on that date, split as the
MJD zero point and the Modified Julian Date. Throws `DomainError` for a year before 4800 BC, a
month outside 1-12, or a day outside the month.
"""
function cal2jd(iy::Integer, im::Integer, id::Integer)
    iy < -4799 && throw(DomainError(iy, "the year must be 4800 BC (-4799) or later."))
    (1 <= im <= 12) || throw(DomainError(im, "the month must be between 1 and 12."))
    ly = (im == 2) && _is_leap_year(iy)
    (1 <= id <= _MONTH_DAYS[im] + ly) ||
        throw(DomainError(id, "the day must be between 1 and $(_MONTH_DAYS[im] + ly)."))
    my    = div(im - 14, 12)
    iypmy = iy + my
    djm = div(1461 * (iypmy + 4800), 4) + div(367 * (im - 2 - 12 * my), 12) -
          div(3 * div(iypmy + 4900, 100), 4) + id - 2432076
    return MJD_EPOCH, Float64(djm)
end

"""
    jd2cal(dj1, dj2) -> (year, month, day, fraction of day)

ERFA's `eraJd2cal`: a two-part Julian date to a Gregorian date and fraction of a day. The date may
be split between the parts in any way; the fraction is summed with compensation (Klein 2006), so
the split costs no precision. Throws `DomainError` below JD -68569.5 (4713 BC January 1) or above
1e9.
"""
function jd2cal(dj1, dj2)
    dj = dj1 + dj2
    (dj < -68569.5 || dj > 1e9) &&
        throw(DomainError(dj, "the Julian Date must be between -68569.5 and 1e9."))

    # Separate day and fraction, -0.5 <= fraction < 0.5.
    d  = _dnint(dj1)
    f1 = dj1 - d
    jd = Int(d)
    d  = _dnint(dj2)
    f2 = dj2 - d
    jd += Int(d)

    # f1 + f2 + 0.5 by compensated summation.
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
    if (f - 1) >= -eps(Float64) / 4
        t = s - 1
        cs += (s - t) - 1
        s = t
        f = s + cs
        if -eps(Float64) / 2 < f
            jd += 1
            f = max(f, zero(f))
        end
    end

    # The day number as a Gregorian date.
    l = jd + 68569
    n = div(4l, 146097)
    l -= div(146097n + 3, 4)
    i = div(4000 * (l + 1), 1461001)
    l -= div(1461i, 4) - 31
    k = div(80l, 2447)
    id = l - div(2447k, 80)
    l = div(k, 11)
    im = k + 2 - 12l
    iy = 100 * (n - 49) + i + l
    return iy, im, id, f
end

"""
    d2tf(ndp, days) -> (sign, hours, minutes, seconds, fraction)

ERFA's `eraD2tf`: an interval in days to hours, minutes, seconds and a fraction of a second in
units of 10^-ndp, rounded to that resolution. `sign` is `'+'` or `'-'`.
"""
function d2tf(ndp::Integer, days)
    sign = days >= 0 ? '+' : '-'
    a = SECONDS_IN_DAY * abs(days)
    if ndp < 0                                  # pre-round if coarser than a second
        nrs = 1
        for n in 1:(-ndp)
            nrs *= (n == 2 || n == 4) ? 6 : 10
        end
        rs = Float64(nrs)
        a = rs * _dnint(a / rs)
    end
    nrs = 1
    for _ in 1:ndp
        nrs *= 10
    end
    rs = Float64(nrs)
    rm = rs * 60.0
    rh = rm * 60.0
    a = _dnint(rs * a)
    ah = trunc(a / rh); a -= ah * rh
    am = trunc(a / rm); a -= am * rm
    as = trunc(a / rs)
    af = a - as * rs
    return sign, Int(ah), Int(am), Int(as), Int(af)
end

# ────────────────────── date and time of day, with leap seconds ─────────────────────
#
# On a UTC day that ends with a leap second the day is 86401 s long, and the final minute 61 s.
# `_leap_seconds_in_day` finds that from the leap-second table: the change in TAI − UTC between 0h
# today and 0h tomorrow, less any pre-1972 drift (which AstroEpochs does not model, so it is zero).
# A day before the table has none: the step to its first value is not a leap second.

function _leap_seconds_in_day(iy, im, id)
    _in_leap_table(iy, im, id) || return 0.0
    dat0  = _dat(iy, im, id)
    iy2, im2, id2, _ = jd2cal(sum(cal2jd(iy, im, id)), 1.5)   # tomorrow, from its noon, as ERFA
    dat24 = _dat(iy2, im2, id2)
    return dat24 - dat0
end

"""
    d2dtf(scale, ndp, d1, d2) -> (year, month, day, hour, minute, second, fraction)

ERFA's `eraD2dtf`: a two-part Julian date in `scale` (`:utc` or another scale symbol) to a
calendar date and time of day, rounded to `ndp` decimal places of a second. In UTC a time inside a
leap second is 23:59:60.
"""
function d2dtf(scale::Symbol, ndp::Integer, d1, d2)
    a1, b1 = d1, d2
    iy1, im1, id1, fd = jd2cal(a1, b1)

    leap = false
    if scale === :utc
        dleap = _leap_seconds_in_day(iy1, im1, id1)
        leap = abs(dleap) > 0.5
        leap && (fd += fd * dleap / SECONDS_IN_DAY)
    end

    _, h, m, s, f = d2tf(ndp, fd)

    # Rounded past 24 h: tomorrow's date, or 23:59:60 on a leap-second day.
    if h > 23
        iy2, im2, id2, _ = jd2cal(a1 + 1.5, b1 - fd)
        if !leap
            iy1, im1, id1, h, m, s = iy2, im2, id2, 0, 0, 0
        else
            if s > 0
                iy1, im1, id1, h, m, s = iy2, im2, id2, 0, 0, 0
            else
                h, m, s = 23, 59, 60
            end
            if ndp < 0 && s == 60
                iy1, im1, id1, h, m, s = iy2, im2, id2, 0, 0, 0
            end
        end
    end
    return iy1, im1, id1, h, m, s, f
end

"""
    dtf2d(scale, year, month, day, hour, minute, second) -> (d1, d2)

ERFA's `eraDtf2d`: a calendar date and time of day in `scale` to a two-part Julian date, `d1` the
JD at 0h and `d2` the fraction of the day. In UTC, on a day that ends with a leap second, the day
has 86401 s and 23:59 has 61, so 23:59:60.5 is a valid time. Throws `ArgumentError` for an hour,
minute or second outside the day, which ERFA reports as a status.
"""
function dtf2d(scale::Symbol, iy::Integer, im::Integer, id::Integer, ihr::Integer, imn::Integer, sec)
    dj, w = cal2jd(iy, im, id)
    dj += w
    day    = SECONDS_IN_DAY
    seclim = 60.0
    if scale === :utc
        dleap = _leap_seconds_in_day(iy, im, id)
        day += dleap
        (ihr == 23 && imn == 59) && (seclim += dleap)
    end
    (0 <= ihr <= 23) || throw(ArgumentError("Hour must be between 0 and 23. Got: $ihr"))
    (0 <= imn <= 59) || throw(ArgumentError("Minute must be between 0 and 59. Got: $imn"))
    (0 <= sec < seclim) || throw(ArgumentError(
        "Seconds must be >= 0.0 and < $(seclim == 60.0 ? "60.0" : "$(seclim), this being a " *
        "leap-second day in UTC"). Got: $sec"))
    time = (60.0 * (60 * ihr + imn) + sec) / day
    return dj, time
end
