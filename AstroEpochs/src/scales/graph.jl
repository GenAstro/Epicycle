# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# ───────────────────────────── the scale graph ─────────────────────────────
#
# Which scales exist, which pairs are joined directly by a transform, and how a conversion between
# any two is routed. The graph and the routes are Astropy's (astropy.time.core MULTI_HOPS), and the
# transforms are ERFA's: offsets.jl (TT, TAI, TCG, TDB, TCB), tdb.jl (TDB − TT) and
# leap_seconds.jl (UTC).
#
# Adding a scale means: its transforms to and from one existing scale, their two entries in
# SCALE_TRANSFORMS, its symbol in TIME_SCALES, and a MULTI_HOPS route to each scale it is not joined
# to directly. Its tag type goes in AstroEpochs.jl with the others.

const TIME_SCALES = Set([:tt, :tai, :tdb, :utc, :tcg, :tcb])

# Routes between scales that are not joined directly, as in Astropy.
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
    const SCALE_TRANSFORMS = Dict{Tuple{Symbol, Symbol}, Function}

The transform for each pair of directly joined scales. Each takes a two-part Julian date in the
first scale and returns it in the second.
"""
const SCALE_TRANSFORMS = Dict{Tuple{Symbol, Symbol}, Function}(
    (:tai, :tt)  => taitt,
    (:tt, :tai)  => tttai,
    (:tt, :tdb)  => tttdb,
    (:tdb, :tt)  => tdbtt,
    (:tt, :tcg)  => tttcg,
    (:tcg, :tt)  => tcgtt,
    (:tdb, :tcb) => tdbtcb,
    (:tcb, :tdb) => tcbtdb,
    (:utc, :tai) => utctai,                 # leap_seconds.jl: the IERS list, kept current
    (:tai, :utc) => taiutc,
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
    elseif haskey(SCALE_TRANSFORMS, (from, to))
        return [from, to]
    else
        tag = s -> _scale_tag_str(s)
        error("No known time scale conversion path from $(tag(from)) to $(tag(to))")
    end
end

"""
    apply_transforms(jd1, jd2, from::Symbol, to::Symbol) -> (jd1, jd2)

Convert a two-part Julian date from scale `from` to scale `to`, one transform per hop of the
route, then rebalance the parts as a Time holds them.
"""
function apply_transforms(jd1::Real, jd2::Real, from::Symbol, to::Symbol)
    from === to && return _rebalance(jd1, jd2)
    T = promote_type(typeof(jd1), typeof(jd2))
    a1, a2 = T(jd1), T(jd2)
    path = get_conversion_path(from, to)
    for i in 1:(length(path) - 1)
        a1, a2 = SCALE_TRANSFORMS[(path[i], path[i + 1])](a1, a2)
    end
    return _rebalance(promote(a1, a2)...)
end
