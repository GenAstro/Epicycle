# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT
#
# Transforms between adjacent time scales, ported from ERFA 2.0.1 (BSD 3-Clause; derived from the
# IAU SOFA library): src/taitt.c, tttai.c, tttcg.c, tcgtt.c, tttdb.c, tdbtt.c, tdbtcb.c and
# tcbtdb.c. See THIRD_PARTY_NOTICES.md. These are the routines Astropy converts scales with;
# test_correctness_erfa_parity.jl holds them to pyerfa's own output.
#
# ─────────────────────────── time-scale transforms ────────────────────────────
#
# Scale transforms: each takes a two-part Julian date in the source scale and returns it in the
# destination scale. As in ERFA, the correction goes into whichever part is smaller in magnitude,
# so the larger part passes through untouched and no precision is lost to summing the two.
# UTC ↔ TAI is in leap_seconds.jl, TDB − TT in tdb.jl; which pairs are joined, and how a longer
# conversion is routed, is in graph.jl.
#
# To add a scale, write its transforms to and from one existing scale here and add both to
# SCALE_TRANSFORMS in graph.jl.

const TTMTAI = 32.184                      # TT − TAI, s, exact by definition
const _DJM77 = 43144.0                     # MJD of 1977-01-01
const ELG    = 6.969290134e-10             # L_G = 1 − d(TT)/d(TCG)   (IAU 2000 Resolution B1.9)
const ELB    = 1.550519768e-8              # L_B = 1 − d(TDB)/d(TCB)  (IAU 2006 Resolution B3)
const TDB0   = -6.55e-5                    # TDB − TCB at 1977-01-01T00:00:32.184 TT, s

# Add `δ` days to the smaller part of a two-part date.
@inline _add_small(a1, a2, δ) = abs(a1) > abs(a2) ? (a1, a2 + δ) : (a1 + δ, a2)

# ── TAI ↔ TT: a fixed offset ────────────────────────────────────────────────

"eraTaitt: TAI to TT."
taitt(tai1, tai2) = _add_small(tai1, tai2, TTMTAI / SECONDS_IN_DAY)

"eraTttai: TT to TAI."
tttai(tt1, tt2) = _add_small(tt1, tt2, -TTMTAI / SECONDS_IN_DAY)

# ── TT ↔ TCG: a linear rate from 1977-01-01T00:00:32.184 TT ─────────────────

const _T77T = _DJM77 + TTMTAI / SECONDS_IN_DAY     # that epoch as an MJD

"eraTttcg: TT to TCG."
function tttcg(tt1, tt2)
    elgg = ELG / (1.0 - ELG)
    return abs(tt1) > abs(tt2) ?
        (tt1, tt2 + ((tt1 - MJD_EPOCH) + (tt2 - _T77T)) * elgg) :
        (tt1 + ((tt2 - MJD_EPOCH) + (tt1 - _T77T)) * elgg, tt2)
end

"eraTcgtt: TCG to TT."
function tcgtt(tcg1, tcg2)
    return abs(tcg1) > abs(tcg2) ?
        (tcg1, tcg2 - ((tcg1 - MJD_EPOCH) + (tcg2 - _T77T)) * ELG) :
        (tcg1 - ((tcg2 - MJD_EPOCH) + (tcg1 - _T77T)) * ELG, tcg2)
end

# ── TDB ↔ TCB: a linear rate from the same epoch, and the constant TDB0 ─────

const _T77TD = MJD_EPOCH + _DJM77                  # 1977-01-01 as a JD
const _T77TF = TTMTAI / SECONDS_IN_DAY             # and 32.184 s, as a fraction of a day

"eraTdbtcb: TDB to TCB."
function tdbtcb(tdb1, tdb2)
    tdb0 = TDB0 / SECONDS_IN_DAY
    elbb = ELB / (1.0 - ELB)
    if abs(tdb1) > abs(tdb2)
        d = _T77TD - tdb1
        f = tdb2 - tdb0
        return tdb1, f - (d - (f - _T77TF)) * elbb
    else
        d = _T77TD - tdb2
        f = tdb1 - tdb0
        return f - (d - (f - _T77TF)) * elbb, tdb2
    end
end

"eraTcbtdb: TCB to TDB."
function tcbtdb(tcb1, tcb2)
    tdb0 = TDB0 / SECONDS_IN_DAY
    if abs(tcb1) > abs(tcb2)
        d = tcb1 - _T77TD
        return tcb1, tcb2 + tdb0 - (d + (tcb2 - _T77TF)) * ELB
    else
        d = tcb2 - _T77TD
        return tcb1 + tdb0 - (d + (tcb1 - _T77TF)) * ELB, tcb2
    end
end

# ── TT ↔ TDB: TDB − TT from tdb.jl ──────────────────────────────────────────
#
# As Astropy does, TDB − TT is evaluated at the epoch being converted, in whichever of the two
# scales it is given; the two differ by under two milliseconds, which changes TDB − TT by less than
# a picosecond, so neither direction iterates.

"eraTttdb: TT to TDB, with TDB − TT from `tdb_minus_tt`."
tttdb(tt1, tt2) = _add_small(tt1, tt2, tdb_minus_tt(tt1, tt2) / SECONDS_IN_DAY)

"eraTdbtt: TDB to TT, with TDB − TT from `tdb_minus_tt`."
tdbtt(tdb1, tdb2) = _add_small(tdb1, tdb2, -tdb_minus_tt(tdb1, tdb2) / SECONDS_IN_DAY)
