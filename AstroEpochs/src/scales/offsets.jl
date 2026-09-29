# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT
#
# Adapted from Tempo.jl v1.3.1, src/offset.jl (MIT License).
# Copyright (c) 2022 Andrea Pasquale and Michele Ceresoli. See THIRD_PARTY_NOTICES.md.
# The functions, constants and arithmetic are Tempo's, unchanged; test_correctness_tempo_parity.jl
# holds them to Tempo's own outputs.

# ─────────────────────────── time-scale offsets ────────────────────────────
#
# Scale transforms: each function takes the epoch as seconds since J2000 in the source scale and
# returns what to add to reach the destination scale, in seconds. UTC ↔ TAI is in leap_seconds.jl;
# which pairs are joined, and how a longer conversion is routed, is in graph.jl.
#
# To add a scale, or a better model for an existing pair, write the offset function here and
# change its entry in OFFSET_TABLE in graph.jl.

# ── TT ↔ TAI: a fixed offset, by definition ─────────────────────────────────

const OFFSET_TAI_TT = 32.184        # s

"TT − TAI, which is 32.184 s at every epoch."
@inline offset_tai2tt(seconds) = OFFSET_TAI_TT

"TAI − TT, which is −32.184 s at every epoch."
@inline offset_tt2tai(seconds) = -OFFSET_TAI_TT

# ── TT ↔ TCG: a linear rate from 1977-01-01T00:00:00 TAI (IAU 2000 Resolution B1.9) ──────

const JD77_SEC = -7.25803167816e8   # 1977-01-01T00:00:32.184 TT, s since J2000
const LG_RATE  = 6.969290134e-10    # L_G, defining constant

"TT − TCG at a TCG epoch, in seconds since J2000."
@inline function offset_tcg2tt(seconds)
    δt = seconds - JD77_SEC
    return -LG_RATE * δt
end

"TCG − TT at a TT epoch, in seconds since J2000."
@inline function offset_tt2tcg(seconds)
    rate = LG_RATE / (1 - LG_RATE)
    δt = seconds - JD77_SEC
    return rate * δt
end

# ── TDB ↔ TCB: a linear rate from the same epoch (IAU 2006 Resolution B3) ─────
#
# The resolution also has a constant term, TDB₀ = −6.55e-5 s, which this form omits.

const LB_RATE = 1.550519768e-8      # L_B, defining constant

"TDB − TCB at a TCB epoch, in seconds since J2000."
@inline function offset_tcb2tdb(seconds)
    δt = seconds - JD77_SEC
    return -LB_RATE * δt
end

"TCB − TDB at a TDB epoch, in seconds since J2000."
@inline function offset_tdb2tcb(seconds)
    rate = LB_RATE / (1 - LB_RATE)
    δt = seconds - JD77_SEC
    return rate * δt
end

# ── TT ↔ TDB: the one-term periodic approximation ──────────────────────────
#
# TDB − TT ≈ k sin(g + e sin g), with g the Earth's mean anomaly. Accurate to about 40 µs over
# 1900-2100; an observer on the Earth's surface sees up to about 4 µs more. The module docstring
# of AstroEpochs lists the more accurate models that could replace it here.
#
# References: https://www.cv.nrao.edu/~rfisher/Ephemerides/times.html#TDB;
# https://github.com/JuliaAstro/AstroTime.jl/issues/26

const _TDB_K  = 1.657e-3            # s, amplitude
const _TDB_EB = 1.671e-2            # orbital eccentricity of the Earth
const _TDB_M0 = 6.239996            # rad, mean anomaly at J2000
const _TDB_M1 = 1.99096871e-7       # rad/s, mean motion

"TDB − TT at a TT epoch, in seconds since J2000."
@inline function offset_tt2tdb(seconds)
    g = _TDB_M0 + _TDB_M1 * seconds
    return _TDB_K * sin(g + _TDB_EB * sin(g))
end

"TT − TDB at a TDB epoch, in seconds since J2000, by three fixed-point iterations."
@inline function offset_tdb2tt(seconds)
    tt = seconds
    offset = 0
    for _ in 1:3
        g = _TDB_M0 + _TDB_M1 * tt
        offset = -_TDB_K * sin(g + _TDB_EB * sin(g))
        tt = seconds + offset
    end
    return offset
end
