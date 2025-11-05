# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation

__precompile__()
"""
Module containing time system implementations for astronomical times.
Supports high-precision time representations using dual-float Julian Date
storage and conversions between scales such as TT, TAI, UTC, TDB, TCB, TCG.

The API is inspired by AstroPy.Time. The numerics are built on Tempo.jl

Notes
- This module does not handle leap seconds directly. Those are handled in Tempo.jl. 
- This module currently mixes symbols and instances for time scales and formats.
  Future versions will standardize on typed tags (e.g., TT(), TDB(), JD()) to avoid ambiguity.

 Current TDB↔TT conversion model
 - Microsecond-level approximation intended for design/trade studies.
 - Not suitable for navigation-grade timing; 

 Future work (https://github.com/JuliaAstro/AstroTime.jl/issues/26)
 - Fast (ERFA harmonic series):
   Implement Δ_tt2tdb_erfa(TT) using the SOFA/ERFA dtdb trigonometric series
   with identical fundamental arguments and coefficient tables for TT→TDB,
   and invert TDB→TT via a short fixed-point iteration using the same Δ.
 - High-precision (ephemeris-based):
   Compute Δ from relativistic terms using Earth barycentric position/velocity
   (r·v/c^2) and gravitational potential (U/c^2) from ephemeris, with
   optional topocentric corrections; use the same Δ for forward/inverse to
   ensure consistency.
"""
module AstroEpochs

using Printf
using Tempo

using AstroBase 

export Time
export TAI, TT, TDB, UTC, TCB, TCG
export JD, MJD, ISOT

import Base: +,-

const SECONDS_IN_DAY = 86400.0
const MJD_EPOCH = Tempo.DJM0 
const J2000_EPOCH = 2451545.0

abstract type AbstractTimeScale end
abstract type AbstractTimeFormat end

"""
    TT <: AbstractTimeScale

Terrestrial Time scale.
"""
struct TT   <: AbstractTimeScale end

"""
    TDB <: AbstractTimeScale

Barycentric Dynamical Time scale.
"""
struct TDB  <: AbstractTimeScale end

"""
    UTC <: AbstractTimeScale    

Coordinated Universal Time scale.
"""
struct UTC  <: AbstractTimeScale end

"""
    TCB <: AbstractTimeScale

Barycentric Coordinate Time scale.
"""
struct TCB  <: AbstractTimeScale end

"""
    TCG <: AbstractTimeScale

Geocentric Coordinate Time scale.
"""
struct TCG  <: AbstractTimeScale end

"""
    TAI <: AbstractTimeScale

International Atomic Time scale.
"""
struct TAI  <: AbstractTimeScale end

"""
    JD <: AbstractTimeFormat

Julian Date format.
"""
struct JD     <: AbstractTimeFormat end

"""
    MJD <: AbstractTimeFormat

Modified Julian Date format.
"""
struct MJD    <: AbstractTimeFormat end   

"""
    ISOT <: AbstractTimeFormat

ISO 8601 Time format.
"""
struct ISOT   <: AbstractTimeFormat end

# Map tags -> existing Symbol API
@inline _scale_symbol(::TT)  = :tt       # COV_EXCL_LINE
@inline _scale_symbol(::TDB) = :tdb      # COV_EXCL_LINE
@inline _scale_symbol(::UTC) = :utc      # COV_EXCL_LINE
@inline _scale_symbol(::TCB) = :tcb      # COV_EXCL_LINE
@inline _scale_symbol(::TCG) = :tcg      # COV_EXCL_LINE
@inline _scale_symbol(::TAI) = :tai      # COV_EXCL_LINE

@inline _format_symbol(::JD)     = :jd    # COV_EXCL_LINE
@inline _format_symbol(::MJD)    = :mjd   # COV_EXCL_LINE
@inline _format_symbol(::ISOT)   = :isot  # COV_EXCL_LINE

# Map  Symbol -> tags
@inline _scale_tag(::Val{:tt})  = TT()
@inline _scale_tag(::Val{:tdb}) = TDB()
@inline _scale_tag(::Val{:utc}) = UTC()
@inline _scale_tag(::Val{:tcb}) = TCB()
@inline _scale_tag(::Val{:tcg}) = TCG()
@inline _scale_tag(::Val{:tai}) = TAI()

@inline _format_tag(::Val{:jd})   = JD()
@inline _format_tag(::Val{:mjd})  = MJD()
@inline _format_tag(::Val{:isot}) = ISOT()

""" 
    _typename_paren(x) = string(nameof(typeof(x)), "()")

Print a type name with parentheses for display.
"""
_typename_paren(x) = string(nameof(typeof(x)), "()")

# Define supported time scales and formats
const TIME_FORMATS = Set([:jd, :mjd, :isot])
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
    (:tt, :tai)  => Tempo.offset_tt2tai,
    (:tai, :tt)  => Tempo.offset_tai2tt,
    (:tt, :tdb)  => Tempo.offset_tt2tdb,
    (:tdb, :tt)  => Tempo.offset_tdb2tt,
    (:tai, :utc) => Tempo.offset_tai2utc,
    (:utc, :tai) => Tempo.offset_utc2tai,
    (:tcg, :tt)  => Tempo.offset_tcg2tt,
    (:tcb,:tdb)  => Tempo.offset_tcb2tdb,
    (:tt,:tcg)   => Tempo.offset_tt2tcg,
    (:tdb,:tcb)  => Tempo.offset_tdb2tcb,
)

"""
    Time{T<:Real}

High-precision astronomical epoch represented as a split Julian Date with an associated time scale and format.

Fields
- _jd1::Real — first component of the split Julian Date (access via `t.jd1`).
- _jd2::Real — second component of the split Julian Date (access via `t.jd2`).
- scale::Symbol — time scale tag, one of: :tt, :tai, :tdb, :utc, :tcg, :tcb.
- format::Symbol — time format tag, one of: :jd, :mjd, :isot.

# Notes:
- Invariant: `jd1 + jd2` equals the epoch’s Julian Date. Internally, `_rebalance` keeps `_jd2 ∈ [-0.5, 0.5)`.
- Property access performs on-demand conversions:
  - `t.tt`, `t.tdb`, `t.utc`, … return a new Time converted to that scale.
  - `t.jd` and `t.mjd` return numeric date values; `t.isot` returns an ISO 8601 string.
- Constructors accept Symbols (shown) and also typed tags (e.g., `TT(), TDB(), JD(), MJD(), ISOT()`).
- Time arithmetic uses days: `t + 1.0` advances by one day; `t2 - t1` returns a Real (days).
- Validation is strict: scales and formats must be supported; ISO input strings must match YYYY-MM-DDTHH:MM:SS.sss.
- Time currently mixes symbols and instances for time scales and formats.
  Future versions will standardize on typed tags (e.g., TT(), TDB(), JD()) to avoid ambiguity.
- To maintain precision, jd2 should remain small (~< 1.0); use jd1 for large offsets.

  # Examples
```julia
using AstroEpochs

# Construct from JD whole (TT scale)
t1 = Time(2451545.25, TT(), JD())

# Construct from JD parts (TT scale)
t1 = Time(2451545.0, 0.25, TT(), JD())

# Convert scale and format
t2 = t1.tdb
println(typeof(t2))
println(t2.mjd)      # numeric MJD
println(t2.isot)     # ISO 8601 string

# Construct from MJD value (TDB scale)
t3 = Time(58000.0, TDB(), MJD())
println(t3.jd)       # numeric JD

# Construct from ISO string (UTC scale), tagged as MJD format
t4 = Time("2017-01-01T00:00:00.000", UTC(), ISOT())
println(t4.jd)       # numeric JD

# Arithmetic (days)
t5 = t3 + 2.0
println(t5.jd - t3.jd)

# Difference (days), same scale/format required
dt = t5 - t3
println(dt)
```
"""
struct Time{T<:Real}
    _jd1::T
    _jd2::T
    scale::Symbol
    format::Symbol

    # Core typed constructor (validating), infers T from jd1/jd2
    function Time(jd1::T, jd2::T, scale::Symbol, format::Symbol) where {T<:Real}
        jd1, jd2 = _rebalance(jd1, jd2)
        _validate_scale(scale)
        _validate_format(format)
        new{T}(jd1, jd2, scale, format)
    end
end

"""
    Time(jd1::Real, jd2::Real, scale::Symbol, format::Symbol)

Promotes to `T = promote_type(typeof(jd1), typeof(jd2))` and constructs `Time{T}`.
"""
function Time(jd1::Real, jd2::Real, scale::Symbol, format::Symbol)
    T = promote_type(typeof(jd1), typeof(jd2))
    return Time(T(jd1), T(jd2), scale, format)
end

"""
    Time(jd1, jd2, scale, format)

Construct time given single-value Julian Date with validation.
"""
function Time(jd1, jd2, scale, format)
    msg = IOBuffer()
    bad = false
    if !(jd1 isa Real); print(msg, "jd1 must be Real; got ", typeof(jd1), ". "); 
        bad = true end
    if !(jd2 isa Real); print(msg, "jd2 must be Real; got ", typeof(jd2), ". "); 
        bad = true end
    if !(scale isa Symbol); print(msg, "scale must be Symbol; got ", typeof(scale), ". "); 
        bad = true end
    if !(format isa Symbol); print(msg, "format must be Symbol; got ", 
        typeof(format), ". "); bad = true end
    if !bad
        # This should be an unreachable error, throw internal error if encountered
        error("Internal error: invalid time input types.") # COV_EXCL_LINE
    end
    throw(ArgumentError("Time: invalid time input types. " * String(take!(msg))))
end

"""
    Time(value::Real, scale::Symbol, format::Symbol)

Construct time given single-value Julian Date.
"""
function Time(value::Real, scale::Symbol, format::Symbol)
    _validate_scale(scale)
    _validate_format(format)
    _validate_inputcoupling(value,format)
    return _from_format(value, scale, format)
end

"""
    Time(isostr::String, scale::Symbol, format::Symbol) → Time

Construct time given time in ISOT format.
""" 
function Time(isostr::String, scale::Symbol, format::Symbol)

    _validate_scale(scale)
    _validate_format(format)
    _validate_inputcoupling(isostr,format)
    
    y, m, d, h, mi, s = _isot_to_date(isostr)
    jd1, jd2 = calhms2jd_prec(y, m, d, h, mi, s)

    jd1, jd2 = _rebalance(jd1, jd2)
    return Time(jd1, jd2, scale, format)
end

""" 
    Time(jd1::Real, jd2::Real, s::AbstractTimeScale, f::AbstractTimeFormat) -> Time

Construct time given JD parts and typed scale/format tags.
"""
Time(jd1::Real, jd2::Real, s::AbstractTimeScale, f::AbstractTimeFormat) =
    Time(jd1, jd2, _scale_symbol(s), _format_symbol(f))


"""
    Time(value::Real, s::AbstractTimeScale, f::AbstractTimeFormat) -> Time

Construct time given single-value JD and typed scale/format tags.
"""
Time(value::Real, s::AbstractTimeScale, f::AbstractTimeFormat) =
    Time(value, _scale_symbol(s), _format_symbol(f))

""" 
    _typename_paren(x) = string(nameof(typeof(x)), "()")

Helper function to print type name with parentheses in show(). 
"""
@inline function _scale_tag_str(s::Symbol)
    try
        return _typename_paren(_scale_tag(Val(s)))
    catch
        return ":" * String(s)
    end
end

""" 
    _typename_paren(x) = string(nameof(typeof(x)), "()")

Helper function to print type name with parentheses in show(). 
"""
@inline function _format_tag_str(f::Symbol)
    try
        return _typename_paren(_format_tag(Val(f)))
    catch
        return ":" * String(f)
    end
end

"""
    Time(isostr::String, s::AbstractTimeScale, f::AbstractTimeFormat) -> Time

Construct time given ISOT string and typed scale/format tags.
"""
Time(isostr::String, s::AbstractTimeScale, f::AbstractTimeFormat) =
    Time(isostr, _scale_symbol(s), _format_symbol(f))

"""
    Base.show(io::IO, ::MIME"text/plain", t::Time)

Pretty, stable text/plain rendering for Time.
"""
function Base.show(io::IO, ::MIME"text/plain", t::Time)
    value = if t.format == :jd
        t.jd
    elseif t.format == :mjd
        t.mjd
    elseif t.format == :isot
        t.isot
    else
        "Unknown format: $(t.format)"
    end

    scale_str = try
        _typename_paren(_scale_tag(Val(t.scale)))
    catch
        string(t.scale)
    end
    format_str = try
        _typename_paren(_format_tag(Val(t.format)))
    catch
        string(t.format)
    end

    println(io, "AstroEpochs.Time")
    println(io, "  value  = ", value)
    println(io, "  scale  = ", scale_str)
    println(io, "  format = ", format_str)
end

"""
    Base.show(io::IO, t::Time)

Delegate to the text/plain renderer so print/println/sprint(show, t) use the same output.
"""
function Base.show(io::IO, t::Time)
    show(io, MIME"text/plain"(), t)
end

"""
    _validate_format(format::Symbol)

Validate time format against supported formats
"""
function _validate_format(format::Symbol)
    if format ∉ TIME_FORMATS
        supported = join(map(_format_tag_str, collect(TIME_FORMATS)), ", ")
        throw(ArgumentError("Time: invalid time format $( _format_tag_str(format) ). Supported: [$supported]"))
    end
end

"""
    function _validate_scale(scale::Symbol)

Validate time scale against supported scales
"""        
function _validate_scale(scale::Symbol)
    if scale ∉ TIME_SCALES
        supported = join(map(_scale_tag_str, collect(TIME_SCALES)), ", ")
        throw(ArgumentError("Time: invalid time scale $( _scale_tag_str(scale) ). Supported: [$supported]"))
    end
end

"""
    _validate_inputcoupling(value::Any, format::Symbol)

Validate that `value` is compatible with the specified time `format`.
"""
function _validate_inputcoupling(value::Any, format::Symbol)
    if format === :isot
        if !(value isa AbstractString)
            throw(ArgumentError("Time: format $( _format_tag_str(format) ) requires AbstractString input. Got $(typeof(value))"))
        end
    elseif format === :jd || format === :mjd
        if !(value isa Real)
            throw(ArgumentError("Time: format $( _format_tag_str(format) ) requires Real input. Got $(typeof(value))"))
        end
    else
        throw(ArgumentError("Time: unsupported format $( _format_tag_str(format) ) in this constructor."))
    end
end

"""
    Base.getproperty(t::Time, name::Symbol)

Provide field access and on-demand conversions for Time via property syntax.

Arguments
- t::Time — the epoch object.
- name::Symbol — one of:
  - Raw fields: :jd1, :jd2, :scale, :format
  - Format accessors: :jd (Real), :mjd (Real), :isot (String)
  - Scale conversions: any of `:tt, :tai, :tdb, :utc, :tcg, :tcb` to return a new Time in that scale.

Returns
- Real for :jd, :mjd, :jd1, :jd2
- Symbol for :scale, :format
- String for :isot
- Time for scale symbols (converted, preserving current `format`)

Notes
- Scale conversions are performed through apply_offsets following 
  the configured conversion graph.
- When converting scales, `jd1/jd2` are rebalanced to keep invariants.
- Unknown properties throw an error.
- Scale conversions preserve format and return a new Time object;
  the original is unchanged.
- Format computations return a time value (as opposed to a new Time struct).

# Examples
```julia
using AstroEpochs

t = Time("2024-02-29T12:34:56.123", TAI(), ISOT())

# Raw fields
t.jd1; t.jd2; t.scale; t.format

# Format accessors
t.jd        # numeric JD
t.mjd       # numeric MJD
t.isot      # ISO 8601 string

# Scale conversion (preserves current format :isot)
tt_time = t.tt
tai_time = tt_time.tai
```
"""
function Base.getproperty(t::Time, name::Symbol)

    # Julian date elements jd1 and jd2
    if name === :jd1
        return getfield(t, :_jd1)
    elseif name === :jd2
        return getfield(t, :_jd2)
    elseif name === :scale || name === :format
        return getfield(t, name)
    end

    # Time format conversion
    if name === :isot
        return _to_isot(t)
    elseif name === :jd
        return t.jd1 + t.jd2
    elseif name === :mjd 
        return _to_mjd(t)
    end

    # Time scale conversion 
    if name in TIME_SCALES
        newjd1, newjd2 = apply_offsets(t.jd1, t.jd2, t.scale, name)
        return Time(newjd1, newjd2, name, t.format)
    end

    # Fall back to normal field access
    error("Unknown property `$name` for Time object.")

end

"""
    _to_mjd(t::Time)

Return time value in Modified Julian Date (MJD) format.
"""
@inline function _to_mjd(t::Time)
    return (t.jd1 - MJD_EPOCH) + t.jd2
end

"""
    _to_isot(t::Time)

Return time value in ISOT format rounding to millisecond.
"""
function _to_isot(t::Time)
    jd = t.jd1 + t.jd2
    y, m, d, fd = Tempo.jd2cal(t.jd1,t.jd2)
    h, mi, s = Tempo.fd2hms(fd)
    return @sprintf("%04d-%02d-%02dT%02d:%02d:%06.3f", y, m, d, h, mi, s)
end

"""
    _from_format(value::Real, scale::Symbol, format::Symbol) -> Time{T}

Construct a Time{T} from a JD/MJD numeric value, preserving T = typeof(value).
"""
function _from_format(value::Real, scale::Symbol, format::Symbol)
    T = typeof(value)
    if format == :jd
        jd1 = floor(T, value)
        jd2 = T(value - jd1)
        jd1, jd2 = _rebalance(jd1, jd2)
        return Time(jd1, jd2, scale, :jd)
    elseif format == :mjd
        jd1 = floor(T, value) + T(MJD_EPOCH)
        jd2 = T(value - floor(T, value))
        jd1, jd2 = _rebalance(jd1, jd2)
        return Time(jd1, jd2, scale, :mjd)
    else
        throw(ArgumentError("Internal Error: Unsupported time format: $format")) # COV_EXCL_LINE
    end
end

"""
    function _rebalance(jd1::Real, jd2::Real)

Rebalance `jd1` and `jd2` so that `-0.5 ≤ jd2 ≤ 0.5`, 
with the remainder placed in `jd1`.
"""
function _rebalance(jd1::T, jd2::T) where {T<:Real}
    half = T(0.5)
    oneT = one(T)

    # Shift up/down by whole days until jd2 ∈ [-0.5, 0.5)
    # TODO. Optimize this for large offsets?
    while jd2 >= half
        jd1 += oneT
        jd2 -= oneT
    end
    while jd2 < -half
        jd1 -= oneT
        jd2 += oneT
    end

    return jd1, jd2
end

""" 
    function _rebalance(a::Real, b::Real)

# Promote to a common T and reuse the typed method
"""
@inline _rebalance(a::Real, b::Real) = _rebalance(promote(a, b)...)

"""
    function -(t2::Time, t1::Time)

  Subtract two `Time` objects, returning the difference in days.
"""
function -(t2::Time, t1::Time)::Real
    if t1.scale != t2.scale
        error("Cannot subtract Times with different scales (", t1.scale, " vs ", t2.scale, ")")
    end
    if t1.format != t2.format
        error("Cannot subtract Times with different formats (", t1.format, " vs ", t2.format, ")")
    end
    return (t2.jd1 - t1.jd1) + (t2.jd2 - t1.jd2)
end

"""
    -(t::Time{T}, dt::Real) -> Time{PT}

Subtract days `dt` from `t` returning new time object.
"""
function Base.:-(t::Time{T}, dt::Real) where {T<:Real}
    δ1, δ2 = _rebalance(-T(dt), zero(T))
    return Time(getfield(t, :_jd1) + δ1, getfield(t, :_jd2) + δ2, t.scale, t.format)
end

"""
    function +(t::Time, dt::Real)

Add a scalar (in days) to a Time object, returning a new Time struct.
"""
function +(t::Time, dt::Real)
    # Keep AD types by pairing dt with zero(dt)
    δ1, δ2 = _rebalance(dt, zero(dt))
    jd1 = t.jd1 + δ1
    jd2 = t.jd2 + δ2
    return Time(jd1, jd2, t.scale, t.format)
end

"""
    function +(t::Time, dt::Real)

Add a scalar (in days) to a Time object, commutative version.
"""
function +(dt::Real, t::Time)
    return t + dt
end

"""
    ==(a::Time, b::Time)

Equality operator for Time.  Same scale and format, and identical jd1/jd2 values.
"""
function Base.:(==)(a::Time, b::Time)
    a.scale === b.scale || return false
    a.format === b.format || return false
    Ta = typeof(getfield(a, :_jd1))
    Tb = typeof(getfield(b, :_jd1))
    PT = promote_type(Ta, Tb)
    va = PT(getfield(a, :_jd1)) + PT(getfield(a, :_jd2))
    vb = PT(getfield(b, :_jd1)) + PT(getfield(b, :_jd2))
    return va == vb
end

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
    
"""
    _isot_to_date(isostr::String) → Tuple{Int, Int, Int, Int, Int, Real}

Validate and parse an ISO 8601 string of the form `"YYYY-MM-DDTHH:MM:SS.sss"`.
Returning calendar fields: `(year, month, day, hour, minute, second)`.
"""
function _isot_to_date(isostr::String)
    # Match ISO 8601 format with optional fractional seconds
    m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)$", isostr)
    if m === nothing
        throw(ArgumentError("Time: Invalid ISO 8601 format: $isostr"))
    end

    y  = parse(Int, m.captures[1])
    mth = parse(Int, m.captures[2])
    d  = parse(Int, m.captures[3])
    h  = parse(Int, m.captures[4])
    mi = parse(Int, m.captures[5])
    s  = parse(Float64, m.captures[6])

    # Semantic range checks
    if !(1 <= mth <= 12)
        throw(ArgumentError("Month must be between 1 and 12. Got: $mth"))
    end
    if !(1 <= d <= 31)
        throw(ArgumentError("Day must be between 1 and 31. Got: $d"))
    end
    if !(0 <= h < 24)
        throw(ArgumentError("Hour must be between 0 and 23. Got: $h"))
    end
    if !(0 <= mi < 60)
        throw(ArgumentError("Minute must be between 0 and 59. Got: $mi"))
    end
    if !(0.0 <= s < 60.0)
        throw(ArgumentError("Seconds must be >= 0.0 and < 60.0. Got: $s"))
    end
    
    return (y, mth, d, h, mi, s) 
end

"""
    calhms2jd_prec(Y::I, M::I, D::I, h::I, m::I, sec::N) where {I <: Integer, N <: Number}

Convert calendar date and time to Julian Date
"""
function calhms2jd_prec(Y::I, M::I, D::I, h::I, m::I, sec::N) where {I <: Integer, N <: Number}
    # TODO: Submit Patch to Tempo.jl
    j2000_epoch, daysfrom_j2000 = Tempo.cal2jd(Y, M, D)
    frac_of_day = Tempo.hms2fd(h, m, sec)

    return Float64(j2000_epoch + daysfrom_j2000), frac_of_day - 0.5
end

end 


