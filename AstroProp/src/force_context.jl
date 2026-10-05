# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0
#
# What the forces in one evaluation of the equations of motion share.
#
# Every force in a ForceModel is evaluated at the same epoch, and several need the same things
# there: the epoch in TDB, TT and UTC, the central body's rotation into its fixed axes, and the
# positions of the Sun and Moon. Each force used to compute its own. Converting TT to TDB runs the
# full TDB−TT series, and with gravity, drag, relativity, third bodies, SRP and tides each doing
# it, time conversion was a third of the right-hand side in the one-day LEO benchmark
# (EpicycleEnterprise/benchmark/full_force), and the Earth's rotation, computed three times,
# another quarter.
#
# `_build_odes!` hands the forces a `ForceContext` in `params`, reset at each evaluation. A force
# asks it through `force_epoch`, `force_rotation`, `force_position` and `force_state`, which
# compute a value on first use in that evaluation and return the same value after that. With no
# context, or one for another epoch, they compute the value directly, so a force gives the same
# result inside a propagation and out of it.
#
# A context belongs to one propagation and is not safe to share between threads.

mutable struct ForceContext
    time::Union{Nothing, Time{Float64}}
    epoch::Union{Nothing, EpochScales{Float64}}
    rotation_keys::Vector{Tuple{UInt, Int}}
    rotations::Vector{SMatrix{6,6,Float64,36}}
    position_keys::Vector{Tuple{Int, Int}}
    positions::Vector{SVector{3,Float64}}
    state_keys::Vector{Tuple{Int, Int}}
    states::Vector{SVector{6,Float64}}
end

ForceContext() = ForceContext(nothing, nothing,
                              Tuple{UInt, Int}[], SMatrix{6,6,Float64,36}[],
                              Tuple{Int, Int}[], SVector{3,Float64}[],
                              Tuple{Int, Int}[], SVector{6,Float64}[])

# Start an evaluation at `t`. A time that is not Float64, as when the epoch itself is
# differentiated, leaves the context unused.
function _reset!(ctx::ForceContext, t)
    ctx.time  = t isa Time{Float64} ? t : nothing
    ctx.epoch = nothing
    empty!(ctx.rotation_keys); empty!(ctx.rotations)
    empty!(ctx.position_keys); empty!(ctx.positions)
    empty!(ctx.state_keys);    empty!(ctx.states)
    return ctx
end
_reset!(::Nothing, t) = nothing

# The context in `params` if it is for epoch `t`, otherwise `nothing`.
_context(params::NamedTuple, t::Time) =
    haskey(params, :context) && params.context isa ForceContext && params.context.time === t ?
        params.context : nothing
_context(params, t) = nothing

_epoch_scales(t::Time) = EpochScales(t.tdb.jd, t.tt.jd, t.utc.jd)

"""
    force_epoch(params, t::Time) -> EpochScales

The epoch `t` as TDB, TT and UTC Julian dates (`.tdb`, `.tt`, `.utc`), converted once per
evaluation of the equations of motion. A force reads its time scales here rather than from
`t.tdb` and `t.utc`, which convert again on every call.
"""
function force_epoch(params, t::Time)
    ctx = _context(params, t)
    ctx === nothing && return _epoch_scales(t)
    e = ctx.epoch
    if e === nothing
        e = _epoch_scales(t)
        ctx.epoch = e
    end
    return e
end

"""
    force_rotation(params, orientation, naifid, t::Time) -> SMatrix{6,6}

`body_fixed_rotation(orientation, naifid, t)`, computed once per evaluation for each orientation
and body.
"""
function force_rotation(params, orientation::AbstractOrientationModel, naifid::Integer, t::Time)
    ctx = _context(params, t)
    ctx === nothing && return body_fixed_rotation(orientation, naifid, t)
    key = (objectid(orientation), Int(naifid))
    i = findfirst(==(key), ctx.rotation_keys)
    i === nothing || return ctx.rotations[i]
    M = SMatrix{6,6,Float64,36}(body_fixed_rotation(orientation, naifid, force_epoch(params, t)))
    push!(ctx.rotation_keys, key)
    push!(ctx.rotations, M)
    return M
end

"""
    force_position(params, from::CelestialBody, to::CelestialBody, t::Time) -> SVector{3}

`translate(from, to, jd_tdb)`, km in ICRF axes, looked up once per evaluation for each pair.
"""
function force_position(params, from::CelestialBody, to::CelestialBody, t::Time)
    ctx = _context(params, t)
    ctx === nothing && return SVector{3}(translate(from, to, t.tdb.jd))
    key = (from.naifid, to.naifid)
    i = findfirst(==(key), ctx.position_keys)
    i === nothing || return ctx.positions[i]
    r = SVector{3,Float64}(translate(from, to, force_epoch(params, t).tdb))
    push!(ctx.position_keys, key)
    push!(ctx.positions, r)
    return r
end

"""
    force_state(params, from::CelestialBody, to::CelestialBody, t::Time) -> SVector{6}

`translate_state(from, to, jd_tdb)`, km and km/s in ICRF axes, looked up once per evaluation for
each pair.
"""
function force_state(params, from::CelestialBody, to::CelestialBody, t::Time)
    ctx = _context(params, t)
    ctx === nothing && return SVector{6}(translate_state(from, to, t.tdb.jd))
    key = (from.naifid, to.naifid)
    i = findfirst(==(key), ctx.state_keys)
    i === nothing || return ctx.states[i]
    s = SVector{6,Float64}(translate_state(from, to, force_epoch(params, t).tdb))
    push!(ctx.state_keys, key)
    push!(ctx.states, s)
    return s
end
