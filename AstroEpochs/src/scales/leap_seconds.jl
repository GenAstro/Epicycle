# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Leap seconds: TAI − UTC.
#
# The table comes from the IERS leap-second list as IANA distributes it,
# `leap-seconds.list`, which carries its own expiry date. It is downloaded the first time a
# conversion needs it, stored in a Scratch space, and downloaded again once it has expired, so
# a leap second announced after this release is picked up without a new release. Offline,
# an expired copy is used with a warning, and with no copy at all the built-in table below is
# used; that table ends at the 2017-01-01 leap second.
#
# A list is checked before it is used or stored, so a truncated or altered file is refused rather
# than giving a date the wrong offset until the list expires. Getting a table never stops a
# conversion: a depot that cannot be written and processes downloading at once both still end in
# a table. Each process downloads to its own temporary file.
#
# TAI − UTC changes at 0h UTC on the date in the table. Dates are held as days since J2000,
# JD 2451545.0.
#
# The built-in table is adapted from Tempo.jl v1.3.1, src/leapseconds.jl (MIT License),
# Copyright (c) 2022 Andrea Pasquale and Michele Ceresoli. See THIRD_PARTY_NOTICES.md.
# =============================================================================

using Downloads: download
using Scratch: @get_scratch!
using SHA: sha1

const LEAP_SECONDS_URL = "https://data.iana.org/time-zones/tzdb/leap-seconds.list"

# Seconds to wait for the list. The automatic download happens inside the first UTC conversion,
# which must not hang on a network that drops the request; it falls back to a stored or built-in
# table instead. A refresh the caller asked for can wait longer.
const _LEAP_DOWNLOAD_TIMEOUT = 10.0
const _LEAP_REFRESH_TIMEOUT  = 60.0
const _NTP_EPOCH_JD    = 2415020.5          # 1900-01-01T00:00:00, the list's time origin

struct LeapSecondTable
    starts::Vector{Float64}     # UTC days since J2000 at 0h of each change
    delta::Vector{Float64}      # TAI − UTC from that date, in seconds
    expires::Float64            # UTC days since J2000 after which the list is out of date
end

const _LEAP       = Ref{Union{Nothing, LeapSecondTable}}(nothing)
const _LEAP_LOCK  = ReentrantLock()

"""
    parse_leap_seconds(text::AbstractString) -> LeapSecondTable

Parse an IETF/IANA `leap-seconds.list`: one `NTP-seconds  TAI−UTC` pair per data line, and the
expiry on the line that begins `#@`.

The file is refused with an `ArgumentError` unless its `#h` hash matches and its entries could be
real: every change through 2017 as the built-in table has it, dates in order, and each later change
one second. The hash sits on the last line of the file, so a file cut off anywhere fails it.
"""
function parse_leap_seconds(text::AbstractString)
    updated = expiry = file_hash = nothing
    data = Tuple{SubString{String}, SubString{String}}[]
    for line in eachline(IOBuffer(text))
        s = strip(line)
        isempty(s) && continue
        fields = split(s)
        if fields[1] == "#\$" && length(fields) >= 2
            updated = fields[2]
        elseif fields[1] == "#@" && length(fields) >= 2
            expiry = fields[2]
        elseif fields[1] == "#h" && length(fields) >= 2
            file_hash = lowercase(join(fields[2:end]))
        elseif !startswith(s, "#")
            length(fields) >= 2 || throw(ArgumentError("leap-second list has a line with no offset: $s"))
            push!(data, (fields[1], fields[2]))
        end
    end
    isempty(data) && throw(ArgumentError("leap-second list has no entries"))
    expiry === nothing && throw(ArgumentError("leap-second list has no expiry line (#@)"))
    updated === nothing && throw(ArgumentError("leap-second list has no update line (#\$)"))
    file_hash === nothing && throw(ArgumentError("leap-second list has no hash line (#h); the file may be cut off"))
    computed = bytes2hex(sha1(updated * expiry * join(a * b for (a, b) in data)))
    computed == file_hash || throw(ArgumentError("leap-second list does not match its hash (#h)"))

    starts = [_ntp_to_j2000_days(_list_number(a)) for (a, _) in data]
    delta  = [_list_number(b) for (_, b) in data]
    issorted(starts) && allunique(starts) || throw(ArgumentError("leap-second list is not in date order"))
    _check_leap_entries(starts, delta)
    return LeapSecondTable(starts, delta, _ntp_to_j2000_days(_list_number(expiry)))
end

function _list_number(field)
    value = tryparse(Float64, field)
    (value === nothing || !isfinite(value)) &&
        throw(ArgumentError("leap-second list has $(repr(String(field))) where a number belongs"))
    return value
end

# Every list has the changes the built-in table has, and TAI − UTC has only ever moved a second at a
# time. The hash catches damage; this catches a list that hashes correctly and is still not the list.
function _check_leap_entries(starts, delta)
    builtin = _builtin_leap_table()
    n = length(builtin.starts)
    length(starts) >= n && starts[1:n] == builtin.starts && delta[1:n] == builtin.delta ||
        throw(ArgumentError("leap-second list does not hold the changes through 2017-01-01 " *
                            "that every list has"))
    all(step -> abs(step) == 1.0, diff(delta[n:end])) ||
        throw(ArgumentError("leap-second list changes TAI − UTC by something other than one second"))
    return nothing
end

_ntp_to_j2000_days(ntp_seconds) = _NTP_EPOCH_JD + ntp_seconds / SECONDS_IN_DAY - J2000_EPOCH

_now_j2000_days() = time() / SECONDS_IN_DAY + 2440587.5 - J2000_EPOCH    # Unix epoch as a JD

# The built-in table: each date on which TAI − UTC changed, and its value from that date.
const _BUILTIN_LEAP_SECONDS = (
    (1972, 1, 1, 10.0), (1972, 7, 1, 11.0), (1973, 1, 1, 12.0), (1974, 1, 1, 13.0),
    (1975, 1, 1, 14.0), (1976, 1, 1, 15.0), (1977, 1, 1, 16.0), (1978, 1, 1, 17.0),
    (1979, 1, 1, 18.0), (1980, 1, 1, 19.0), (1981, 7, 1, 20.0), (1982, 7, 1, 21.0),
    (1983, 7, 1, 22.0), (1985, 7, 1, 23.0), (1988, 1, 1, 24.0), (1990, 1, 1, 25.0),
    (1991, 1, 1, 26.0), (1992, 7, 1, 27.0), (1993, 7, 1, 28.0), (1994, 7, 1, 29.0),
    (1996, 1, 1, 30.0), (1997, 7, 1, 31.0), (1999, 1, 1, 32.0), (2006, 1, 1, 33.0),
    (2009, 1, 1, 34.0), (2012, 7, 1, 35.0), (2015, 7, 1, 36.0), (2017, 1, 1, 37.0))

"The built-in table, used when no list has ever been downloaded. It never expires."
function _builtin_leap_table()
    # cal2jd gives 0h on the date as (MJD zero point, MJD); held as days since J2000.
    starts = [(first(cal2jd(y, m, d)) - J2000_EPOCH) + last(cal2jd(y, m, d))
              for (y, m, d, _) in _BUILTIN_LEAP_SECONDS]
    delta  = [Δ for (_, _, _, Δ) in _BUILTIN_LEAP_SECONDS]
    return LeapSecondTable(starts, delta, -Inf)
end

"""
    leap_second_table() -> LeapSecondTable

The table in use, loaded on first call. See the file header for where it comes from.
"""
function leap_second_table()
    t = _LEAP[]
    t === nothing || return t
    return _load_leap_second_table!()
end

function _load_leap_second_table!()
    # Precompiling a package that converts a UTC time must not reach for the network, or an install
    # would depend on it. The table is not kept, so the first real use gets the list.
    ccall(:jl_generating_output, Cint, ()) == 1 && return _builtin_leap_table()
    lock(_LEAP_LOCK) do
        _LEAP[] === nothing || return _LEAP[]
        _LEAP[] = _fetch_leap_second_table()
        return _LEAP[]
    end
end

# Where the downloaded list is kept, across sessions. A depot that cannot be written, as in a
# container or on a shared system, costs only the keeping: the list goes to a temporary folder and
# is downloaded again next session.
function _leap_seconds_path()
    try
        return joinpath(@get_scratch!("leap_seconds"), "leap-seconds.list")
    catch e
        e isa InterruptException && rethrow()
        @warn "AstroEpochs: cannot store the leap-second list in the Julia depot " *
              "($(first(DEPOT_PATH))); it will not be kept for the next session." maxlog = 1
        return joinpath(mktempdir(), "leap-seconds.list")
    end
end

# The stored list, if there is one and it can be read and passes its checks.
function _stored_leap_table(path)
    try
        return isfile(path) ? parse_leap_seconds(read(path, String)) : nothing
    catch e
        e isa InterruptException && rethrow()
        return nothing
    end
end

# The table to use: the stored list while it is current, else a fresh download, else the stored
# list though expired, else the built-in table. The path, URL and timeout are arguments so the
# tests can run each branch against a temporary folder and a URL that fails.
function _fetch_leap_second_table(path = _leap_seconds_path(); url = LEAP_SECONDS_URL,
                                  timeout = _LEAP_DOWNLOAD_TIMEOUT)
    cached = _stored_leap_table(path)
    cached !== nothing && cached.expires > _now_j2000_days() && return cached

    fresh = try
        _download_leap_seconds(path, url, timeout)
    catch e
        e isa InterruptException && rethrow()
        nothing
    end
    fresh !== nothing && return fresh

    # Another process may have stored a current list while this one was failing to.
    latest = _stored_leap_table(path)
    latest !== nothing && latest.expires > _now_j2000_days() && return latest

    if cached !== nothing
        @warn "AstroEpochs: the leap-second list expired and could not be refreshed from " *
              "$(url); using the expired copy. A leap second announced since " *
              "it expired is missing." maxlog = 1
        return cached
    end
    @warn "AstroEpochs: could not download the leap-second list from $(url); " *
          "using the built-in table, which ends at 2017-01-01." maxlog = 1
    return _builtin_leap_table()
end

# Download the list to `path`, parsing it before it replaces the stored copy, and leave no partial
# file behind if either step fails. The temporary file is this process's own and sits beside `path`,
# so the rename replaces the stored copy in one step, and processes downloading at once do not
# touch each other's files.
function _download_leap_seconds(path, url, timeout)
    tmp, io = mktemp(dirname(path))
    close(io)
    try
        download(url, tmp; timeout = timeout)
        table = parse_leap_seconds(read(tmp, String))
        _replace_stored_list(tmp, path)
        return table
    finally
        try
            rm(tmp; force = true)
        catch e
            e isa InterruptException && rethrow()    # a leftover temporary file costs nothing
        end
    end
end

# A rename can lose to another process doing the same one. That is fine if a good list is there.
function _replace_stored_list(tmp, path)
    try
        Base.Filesystem.rename(tmp, path)
    catch e
        _stored_leap_table(path) === nothing && rethrow()
    end
    return nothing
end

"""
    refresh_leap_seconds!() -> LeapSecondTable

Download the IERS leap-second list now and use it for every later conversion.

The list is otherwise downloaded on a session's first UTC conversion, if no current copy is stored,
so this is only needed to pick up a newly announced leap second before the stored list expires, or
in a session that has run past the expiry.

# Returns
The table now in use. Throws when the download fails, rather than falling back as a first
conversion does, because a caller who asks for a refresh wants to know it did not happen.

# Example
```julia
refresh_leap_seconds!()
```
"""
refresh_leap_seconds!() =
    _refresh_leap_seconds!(_leap_seconds_path(), LEAP_SECONDS_URL, _LEAP_REFRESH_TIMEOUT)

# The table in use changes only once the download has parsed, so a failure leaves it as it was.
function _refresh_leap_seconds!(path, url, timeout)
    table = _download_leap_seconds(path, url, timeout)
    lock(_LEAP_LOCK) do
        _LEAP[] = table
    end
    return table
end

"TAI − UTC in seconds at a UTC date given as days since J2000."
function tai_minus_utc(utc_days)
    table = leap_second_table()
    idx = searchsortedlast(table.starts, utc_days)
    idx < 1 && return _before_leap_seconds()
    # The built-in table has no expiry (-Inf); using it already warned, when it was loaded.
    utc_days > table.expires + _DUBIOUS_AFTER_DAYS && isfinite(table.expires) &&
        _past_leap_seconds(table.expires)
    return table.delta[idx]
end

# The warnings are kept out of `tai_minus_utc` because a logging macro contains a try block, which
# automatic differentiation cannot compile, and each shows once per session: a propagation or an
# optimisation converts at every step, and a warning per conversion buries everything else.

# Before 1972 UTC was not an integer offset from TAI, and the list has no value for it.
@noinline function _before_leap_seconds()
    @warn "AstroEpochs: no leap-second data before 1972-01-01; UTC before then uses " *
          "TAI − UTC = 0, which can be wrong by up to ten seconds. Shown once per session." maxlog = 1
    return 0.0
end

# Past the list's expiry a leap second may yet be announced. Every date past it is uncertain, but
# the list expires only months ahead, so warning there would warn on routine mission planning.
# ERFA warns for a year more than five past its release ("dubious year"); this warns five years
# past the list's expiry.
const _DUBIOUS_AFTER_DAYS = 5 * 365.25

@noinline function _past_leap_seconds(expires)
    y, m, d, _ = jd2cal(J2000_EPOCH, expires)
    @warn "AstroEpochs: UTC more than five years after " * @sprintf("%04d-%02d-%02d", y, m, d) *
          ", when the leap-second list expires, assumes no further leap seconds. Shown once " *
          "per session." maxlog = 1
    return nothing
end

# Whether the leap-second table covers a UTC date. Before it, TAI − UTC is taken as 0, and the step
# to the table's first value is a jump at 0h on its first date, not a leap second at the end of the
# day before; treated as a leap second, it made 1971-12-31 a UTC day of 86410 s.
function _in_leap_table(iy, im, id)
    djm0, djm = cal2jd(iy, im, id)
    return (djm0 - J2000_EPOCH) + djm >= first(leap_second_table().starts)
end

"""
    _dat(year, month, day) -> seconds

TAI − UTC at 0h UTC on a Gregorian date, from the leap-second table. It plays the part of ERFA's
`eraDat`, without the pre-1972 drift terms: before 1972 it warns and returns 0.
"""
function _dat(iy, im, id)
    djm0, djm = cal2jd(iy, im, id)
    return tai_minus_utc((djm0 - J2000_EPOCH) + djm)
end

# ─────────────────────────────── UTC ↔ TAI ─────────────────────────────────
#
# Ported from ERFA 2.0.1, src/utctai.c and src/taiutc.c (BSD 3-Clause; derived from the IAU SOFA
# library). See THIRD_PARTY_NOTICES.md.
#
# A UTC two-part date is a quasi-Julian date: on a day that ends with a leap second the fraction of
# the day runs over 86401 SI seconds, so 23:59:60.5 is day fraction 86399.5/86401 and 0h is still a
# whole day. That is ERFA's and Astropy's representation, and so AstroEpochs'.

"""
    utctai(utc1, utc2) -> (tai1, tai2)

ERFA's `eraUtctai`: UTC, as a quasi-Julian date, to TAI. The result keeps the split and order of
the input.
"""
function utctai(utc1, utc2)
    big1 = abs(utc1) >= abs(utc2)
    u1, u2 = big1 ? (utc1, utc2) : (utc2, utc1)

    iy, im, id, fd = jd2cal(u1, u2)
    dat0  = _dat(iy, im, id)                                # TAI − UTC at 0h today
    dat12 = dat0                                            # at 12h: no pre-1972 drift here
    iyt, imt, idt, _ = jd2cal(u1 + 1.5, u2 - fd)
    dat24 = _dat(iyt, imt, idt)                             # at 0h tomorrow

    dlod  = 2.0 * (dat12 - dat0)                            # per-day drift, zero after 1972
    dleap = _in_leap_table(iy, im, id) ?                    # any leap second at the end of today
            dat24 - (dat0 + dlod) : zero(dat24)

    fd *= (SECONDS_IN_DAY + dleap) / SECONDS_IN_DAY         # undo the spread of the leap second
    fd *= (SECONDS_IN_DAY + dlod) / SECONDS_IN_DAY          # pre-1972 UTC seconds to SI seconds

    z1, z2 = cal2jd(iy, im, id)
    a2 = z1 - u1
    a2 += z2
    a2 += fd + dat0 / SECONDS_IN_DAY
    return big1 ? (u1, a2) : (a2, u1)
end

"""
    taiutc(tai1, tai2) -> (utc1, utc2)

ERFA's `eraTaiutc`: TAI to UTC as a quasi-Julian date, by inverting `utctai` in three iterations.
The result keeps the split and order of the input.
"""
function taiutc(tai1, tai2)
    big1 = abs(tai1) >= abs(tai2)
    a1, a2 = big1 ? (tai1, tai2) : (tai2, tai1)
    u1, u2 = a1, a2
    for _ in 1:3
        g1, g2 = utctai(u1, u2)
        u2 += a1 - g1
        u2 += a2 - g2
    end
    return big1 ? (u1, u2) : (u2, u1)
end
