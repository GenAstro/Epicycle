# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Epochs in ephemeris files: the CCSDS forms YYYY-MM-DDThh:mm:ss[.d…][Z] and
# YYYY-DDDThh:mm:ss[.d…][Z] (CCSDS 502.0-B-3 §7.5.10), in any supported time system.

# Time systems an ephemeris file may use, by the name files give them. GPS is not an AstroEpochs
# scale; it is TAI less 19 s exactly, so it is read and written through TAI.
const _TIME_SYSTEMS = Dict("UTC" => :utc, "TAI" => :tai, "TT" => :tt, "TDB" => :tdb,
                           "TCB" => :tcb, "TCG" => :tcg, "GPS" => :gps)
const _GPS_MINUS_TAI = -19.0                                       # s
# Registered systems Epicycle cannot convert, so a file using one is named clearly rather than
# misread.
const _UNSUPPORTED_TIME_SYSTEMS = ("UT1", "GMST", "MET", "MRT", "SCLK")

const _CAL = r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$"
const _DOY = r"^(\d{4})-(\d{3})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$"

_is_leap_year(y) = (y % 4 == 0 && y % 100 != 0) || y % 400 == 0
_days_in(y, m) = m == 2 ? (_is_leap_year(y) ? 29 : 28) : (m in (4, 6, 9, 11) ? 30 : 31)

"""
    _parse_epoch(text, system) -> (time, notes, problem)

Read one epoch in the time system named by `system` (a file's TIME_SYSTEM value, upper case).
`time` is `nothing` when the text cannot be read, and `problem` then says why. `notes` lists
corrections made to read it, such as a space in place of the `T`.
"""
function _parse_epoch(text::AbstractString, system::AbstractString)
    notes = String[]
    s = strip(text)
    if endswith(s, "Z") || endswith(s, "z")
        s = s[1:end-1]
        system == "UTC" || push!(notes, "the Z (UTC) suffix ignored; times are $system, as TIME_SYSTEM says")
    end
    sep = match(r"^(\d{4}-(?:\d{2}-\d{2}|\d{3}))[ t](\d.*)$", s)
    if sep !== nothing
        s = sep[1] * "T" * sep[2]
        push!(notes, "date and time separated by a space or lower-case t; read as T")
    end
    m = match(_CAL, s)
    if m !== nothing
        y, mo, d = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
        (1 <= mo <= 12) || return nothing, notes, "month $mo is not 1 to 12"
        (1 <= d <= _days_in(y, mo)) || return nothing, notes, "day $d does not exist in month $mo of $y"
        h, mi, sec = parse(Int, m[4]), parse(Int, m[5]), parse(Float64, m[6])
    else
        m = match(_DOY, s)
        m === nothing && return nothing, notes,
            "not a CCSDS epoch; expected YYYY-MM-DDThh:mm:ss[.d…] or YYYY-DDDThh:mm:ss[.d…]"
        y, doy = parse(Int, m[1]), parse(Int, m[2])
        ndays = _is_leap_year(y) ? 366 : 365
        (1 <= doy <= ndays) || return nothing, notes, "day of year $doy does not exist in $y"
        mo, d = 1, doy
        while d > _days_in(y, mo)
            d -= _days_in(y, mo)
            mo += 1
        end
        h, mi, sec = parse(Int, m[3]), parse(Int, m[4]), parse(Float64, m[5])
    end
    (0 <= h <= 23) || return nothing, notes, "hour $h is not 0 to 23"
    (0 <= mi <= 59) || return nothing, notes, "minute $mi is not 0 to 59"
    # 60 is a leap second, and only UTC has them (§7.5.10, note).
    (sec < 60 || (system == "UTC" && sec < 61)) ||
        return nothing, notes, "seconds $sec out of range" * (system == "UTC" ? "" : "; only UTC has leap seconds")
    scale = _TIME_SYSTEMS[system]
    if scale === :gps
        jd1, jd2 = AstroEpochs.dtf2d(:tai, y, mo, d, h, mi, sec)
        return Time(jd1, jd2 - _GPS_MINUS_TAI / 86400, :tai, :jd), notes, ""
    end
    jd1, jd2 = AstroEpochs.dtf2d(scale, y, mo, d, h, mi, sec)
    return Time(jd1, jd2, scale, :jd), notes, ""
end

"""
    _format_epoch(t, system, digits) -> String

`t` in the time system `system` as `YYYY-MM-DDThh:mm:ss.d…`, with `digits` decimal places of
seconds. Rounding carries into the minute, hour and date, and a UTC time inside a leap second
reads `:60`.
"""
function _format_epoch(t::Time, system::AbstractString, digits::Integer)
    scale = _TIME_SYSTEMS[system]
    if scale === :gps
        tai = t.tai
        jd1, jd2 = tai.jd1, tai.jd2 + _GPS_MINUS_TAI / 86400
        scale = :tai
    else
        u = getproperty(t, scale)
        jd1, jd2 = u.jd1, u.jd2
    end
    y, mo, d, h, mi, s, f = AstroEpochs.d2dtf(scale, digits, jd1, jd2)
    base = @sprintf("%04d-%02d-%02dT%02d:%02d:%02d", y, mo, d, h, mi, s)
    return digits == 0 ? base : base * "." * lpad(string(f), digits, '0')
end
