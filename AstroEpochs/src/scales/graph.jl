# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# ───────────────────────────── the scale graph ─────────────────────────────
#
# Which scales exist, which pairs are joined directly by an offset function, and how a conversion
# between any two is routed. The offset functions themselves are in offsets.jl (TT, TAI, TCG, TDB,
# TCB) and leap_seconds.jl (UTC).
#
# Adding a scale means: its offset functions to and from one existing scale, their two entries in
# OFFSET_TABLE, its symbol in TIME_SCALES, and a MULTI_HOPS route to each scale it is not joined to
# directly. Its tag type goes in AstroEpochs.jl with the others.

const TIME_SCALES = Set([:tt, :tai, :tdb, :utc, :tcg, :tcb])

# Predefined multi-hop paths between Time scales
const MULTI_HOPS = Dict{Tuple{Symbol, Symbol}, Vector{Symbol}}(
    (:tai, :tcb) => [:tt, :tdb],
    (:tai, :tcg) => [:tt],
    (:tai, :tdb) => [:tt],
    (:tcb, :tcg) => [:tdb, :tt],
    (:tcb, :tt)  => [:tdb],
    (:tcb, :utc) => [:tdb, :tt, :tai],
    (:tcg, :tdb) => [:tt],
    (:tcg, :utc) => [:tt, :tai],
    (:tdb, :utc) => [:tt, :tai],
    (:tt, :utc)  => [:tai],
)

"""
    const OFFSET_TABLE = Dict{Tuple{Symbol, Symbol}, Function}

Maps adjacent pairs time scales to conversion functions.
"""
const OFFSET_TABLE = Dict{Tuple{Symbol, Symbol}, Function}(
    (:tt, :tai)  => offset_tt2tai,
    (:tai, :tt)  => offset_tai2tt,
    (:tt, :tdb)  => offset_tt2tdb,
    (:tdb, :tt)  => offset_tdb2tt,
    (:tai, :utc) => offset_tai2utc,                # leap_seconds.jl: the IERS list, kept current
    (:utc, :tai) => offset_utc2tai,
    (:tcg, :tt)  => offset_tcg2tt,
    (:tcb,:tdb)  => offset_tcb2tdb,
    (:tt,:tcg)   => offset_tt2tcg,
    (:tdb,:tcb)  => offset_tdb2tcb,
)

"""
    get_conversion_path(from::Symbol, to::Symbol) → Vector{Symbol}

Returns conversion path from `from` to `to` time scales. 
"""
function get_conversion_path(from::Symbol, to::Symbol)::Vector{Symbol}
    if from == to
        return [from]
    elseif haskey(MULTI_HOPS, (from, to))
        return vcat(from, MULTI_HOPS[(from, to)], to)
    elseif haskey(MULTI_HOPS, (to, from))
        # reverse the forward path
        revpath = reverse(MULTI_HOPS[(to, from)])
        return vcat(from, revpath, to)
    elseif haskey(OFFSET_TABLE, (from, to))
        return [from, to]
    elseif haskey(OFFSET_TABLE, (to, from))
        return [from, to]  # still valid if OFFSET_TABLE has both directions
    else
        tag = s -> _scale_tag_str(s)
        error("No known time scale conversion path from $(tag(from)) to $(tag(to))")
    end
end

"""
    apply_offsets(jd1, jd2, from::Symbol, to::Symbol)

Applies sequence of scale conversions from `from` to `to`.
"""
function apply_offsets(jd1::Real, jd2::Real, from::Symbol, to::Symbol)
    from === to && return _rebalance(jd1, jd2)

    path = get_conversion_path(from, to)               # small vector; fine
    T = promote_type(typeof(jd1), typeof(jd2))
    jd1T = T(jd1); jd2T = T(jd2)
    J2000_T = T(J2000_EPOCH)
    inv_SECS_PER_DAY_T = inv(T(SECONDS_IN_DAY))

    # Accumulate in days; rebalance once at the end to reduce churn
    # TODO. Bug, offset is different for TCB etc. 
    @inbounds for i in 1:(length(path)-1)
        src = path[i]; dst = path[i+1]
        offset_fn = OFFSET_TABLE[(src, dst)]
        jd = jd1T + jd2T
        seconds_since_j2000 = (jd - J2000_T) * (1 / inv_SECS_PER_DAY_T)  # == * SECONDS_IN_DAY
        off_sec_T = convert(T, offset_fn(seconds_since_j2000))
        jd2T += off_sec_T * inv_SECS_PER_DAY_T
    end

    return _rebalance(jd1T, jd2T)
end
