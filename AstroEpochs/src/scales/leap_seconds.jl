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
# TAI − UTC changes at 0h UTC on the date in the table. Dates are held as days since J2000,
# JD 2451545.0.
#
# The built-in table is adapted from Tempo.jl v1.3.1, src/leapseconds.jl (MIT License),
# Copyright (c) 2022 Andrea Pasquale and Michele Ceresoli. See THIRD_PARTY_NOTICES.md.
# =============================================================================

using Downloads: download
using Scratch: @get_scratch!

const LEAP_SECONDS_URL = "https://data.iana.org/time-zones/tzdb/leap-seconds.list"
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
"""
function parse_leap_seconds(text::AbstractString)
    starts, delta = Float64[], Float64[]
    expires = NaN
    for line in eachline(IOBuffer(text))
        s = strip(line)
        isempty(s) && continue
        if startswith(s, "#@")
            expires = _ntp_to_j2000_days(parse(Float64, split(s)[2]))
        elseif !startswith(s, "#")
            fields = split(s)
            push!(starts, _ntp_to_j2000_days(parse(Float64, fields[1])))
            push!(delta,  parse(Float64, fields[2]))
        end
    end
    isempty(starts) && throw(ArgumentError("leap-second list has no entries"))
    issorted(starts) || throw(ArgumentError("leap-second list is not in date order"))
    isnan(expires) && throw(ArgumentError("leap-second list has no expiry line (#@)"))
    return LeapSecondTable(starts, delta, expires)
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
    # cal2jd gives days from J2000 to noon on the date; the change is at 0h, half a day earlier.
    starts = [Float64(last(cal2jd(y, m, d))) - 0.5 for (y, m, d, _) in _BUILTIN_LEAP_SECONDS]
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
    lock(_LEAP_LOCK) do
        _LEAP[] === nothing || return _LEAP[]
        _LEAP[] = _fetch_leap_second_table()
        return _LEAP[]
    end
end

function _fetch_leap_second_table()
    path   = joinpath(@get_scratch!("leap_seconds"), "leap-seconds.list")
    cached = isfile(path) ? _try_parse(read(path, String)) : nothing
    cached !== nothing && cached.expires > _now_j2000_days() && return cached

    fresh = try
        tmp = download(LEAP_SECONDS_URL, path * ".download")
        table = parse_leap_seconds(read(tmp, String))
        mv(tmp, path; force = true)
        table
    catch e
        e isa InterruptException && rethrow()
        nothing
    end
    fresh !== nothing && return fresh

    if cached !== nothing
        @warn "AstroEpochs: the leap-second list expired and could not be refreshed from " *
              "$(LEAP_SECONDS_URL); using the expired copy. A leap second announced since " *
              "it expired is missing." maxlog = 1
        return cached
    end
    @warn "AstroEpochs: could not download the leap-second list from $(LEAP_SECONDS_URL); " *
          "using the built-in table, which ends at 2017-01-01." maxlog = 1
    return _builtin_leap_table()
end

_try_parse(text) = try parse_leap_seconds(text) catch; nothing end

"""
    refresh_leap_seconds!() -> LeapSecondTable

Download the IERS leap-second list now and use it for every later conversion.

The list is otherwise downloaded on the first UTC conversion and again once it expires, so this
is only needed to pick up a newly announced leap second before the stored list expires.

# Returns
The table now in use. Throws when the download fails, rather than falling back as a first
conversion does, because a caller who asks for a refresh wants to know it did not happen.

# Example
```julia
refresh_leap_seconds!()
```
"""
function refresh_leap_seconds!()
    path = joinpath(@get_scratch!("leap_seconds"), "leap-seconds.list")
    tmp  = download(LEAP_SECONDS_URL, path * ".download")
    table = parse_leap_seconds(read(tmp, String))
    mv(tmp, path; force = true)
    lock(_LEAP_LOCK) do
        _LEAP[] = table
    end
    return table
end

"TAI − UTC in seconds at a UTC date given as days since J2000."
function tai_minus_utc(utc_days)
    table = leap_second_table()
    idx = searchsortedlast(table.starts, utc_days)
    idx < 1 && return _before_leap_seconds(utc_days)
    return table.delta[idx]
end

# Before 1972 UTC was not an integer offset from TAI, and the list has no value for it. Warns on
# every call, as Tempo.jl, whose table this was, did. Kept out of `tai_minus_utc` because a logging macro contains a
# try block, which automatic differentiation cannot compile.
@noinline function _before_leap_seconds(utc_days)
    @warn "AstroEpochs: no leap-second data before 1972-01-01; using TAI − UTC = 0 for the " *
          "UTC date $(utc_days) days from J2000, which can be wrong by up to ten seconds."
    return 0.0
end

"UTC → TAI offset in seconds, for a UTC epoch in seconds since J2000."
offset_utc2tai(seconds) = tai_minus_utc(seconds / SECONDS_IN_DAY)

"TAI → UTC offset in seconds, for a TAI epoch in seconds since J2000."
function offset_tai2utc(seconds)
    tai_days = seconds / SECONDS_IN_DAY
    utc_days = tai_days
    for _ in 1:2                    # the step function settles in one iteration; two is safe
        utc_days = tai_days - tai_minus_utc(utc_days) / SECONDS_IN_DAY
    end
    return -tai_minus_utc(utc_days)
end
