# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

__precompile__()
"""
Module containing time system implementations for astronomical times.
Supports high-precision time representations using dual-float Julian Date
storage and conversions between scales such as TT, TAI, UTC, TDB, TCB, TCG.

The API is inspired by AstroPy.Time. The time-scale offsets and the calendar arithmetic are
adapted from Tempo.jl (MIT); see THIRD_PARTY_NOTICES.md.

Organization. Every algorithm here is one of two kinds, and each has its own folder:
- `scales/` changes which clock an instant is measured on: `offsets.jl` (TT, TAI, TCG, TDB, TCB),
  `leap_seconds.jl` (UTC) and `graph.jl` (which pairs are joined, and how a conversion is routed).
- `formats/` changes how an instant is written: `calendar.jl` (Gregorian ↔ Julian date) and
  `formats.jl` (JD, MJD, ISOT).
This file holds the types, the `Time` struct, its construction, accessors and arithmetic.

Notes
- Leap seconds come from the IERS leap-second list, downloaded on first use and refreshed when
  it expires; see `scales/leap_seconds.jl`.
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

using EpicycleBase

export Time
export TAI, TT, TDB, UTC, TCB, TCG
export JD, MJD, ISOT
export refresh_leap_seconds!

import Base: +,-

include("constants.jl")

# Format transforms that need no Time: the calendar, which the leap-second table also uses.
include("formats/calendar.jl")

# Scale transforms: the offset functions, then the graph that joins them.
include("scales/offsets.jl")
include("scales/leap_seconds.jl")
include("scales/graph.jl")

abstract type AbstractTimeScale end
abstract type AbstractTimeFormat end

"""
    TT()

Terrestrial Time scale.

# Example
```jldoctest
Time(2451545.0, TT(), JD()).scale

# output

:tt
```
"""
struct TT   <: AbstractTimeScale end

"""
    TDB()

Barycentric Dynamical Time scale.

# Example
```jldoctest
Time(2451545.0, TDB(), JD()).scale

# output

:tdb
```
"""
struct TDB  <: AbstractTimeScale end

"""
    UTC()

Coordinated Universal Time scale.

# Example
```jldoctest
Time(2451545.0, UTC(), JD()).scale

# output

:utc
```
"""
struct UTC  <: AbstractTimeScale end

"""
    TCB()

Barycentric Coordinate Time scale.

# Example
```jldoctest
Time(2451545.0, TCB(), JD()).scale

# output

:tcb
```
"""
struct TCB  <: AbstractTimeScale end

"""
    TCG()

Geocentric Coordinate Time scale.

# Example
```jldoctest
Time(2451545.0, TCG(), JD()).scale

# output

:tcg
```
"""
struct TCG  <: AbstractTimeScale end

"""
    TAI()

International Atomic Time scale.

# Example
```jldoctest
Time(2451545.0, TAI(), JD()).scale

# output

:tai
```
"""
struct TAI  <: AbstractTimeScale end

"""
    JD()

Julian Date format.

# Example
```jldoctest
Time(2451545.0, TT(), JD()).format

# output

:jd
```
"""
struct JD     <: AbstractTimeFormat end

"""
    MJD()

Modified Julian Date format.

# Example
```jldoctest
Time(51544.5, TT(), MJD()).format

# output

:mjd
```
"""
struct MJD    <: AbstractTimeFormat end   

"""
    ISOT()

ISO 8601 Time format.

# Example
```jldoctest
Time("2000-01-01T12:00:00.000", TT(), ISOT()).format

# output

:isot
```
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


# Marks the constructor that takes Julian-date parts whatever the format; see `_time_jd`.
struct _JDParts end

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

    # Julian-date parts, whatever the format; `format` only sets how the time displays.
    # Reached through `_time_jd`, never directly by a caller.
    function Time{T}(::_JDParts, jd1::T, jd2::T, scale::Symbol, format::Symbol) where {T<:Real}
        jd1, jd2 = _rebalance(jd1, jd2)
        _validate_scale(scale)
        _validate_format(format)
        new{T}(jd1, jd2, scale, format)
    end
end


"""
    _time_jd(jd1, jd2, scale::Symbol, format::Symbol) -> Time

A `Time` from the two parts of a Julian date, displayed in `format`. The internal way to copy or
shift a time while keeping its format; the public two-part constructor reads its parts in the
format it names.
"""
function _time_jd(jd1::Real, jd2::Real, scale::Symbol, format::Symbol)
    T = promote_type(typeof(jd1), typeof(jd2))
    return Time{T}(_JDParts(), T(jd1), T(jd2), scale, format)
end

"""
    Time(val1::Real, val2::Real, scale::Symbol, format::Symbol)

A time from two parts of a date in `format`, summed: Julian-date parts for `:jd`, Modified
Julian Date parts for `:mjd`, as in AstroPy. Splitting a date into a whole and a fractional part
keeps full precision. An ISOT time is a string; two numbers with `:isot` raise an `ArgumentError`.
The numeric type is `promote_type(typeof(val1), typeof(val2))`.
"""
function Time(val1::Real, val2::Real, scale::Symbol, format::Symbol)
    _validate_format(format)
    format === :isot && throw(ArgumentError(
        "Time: an ISOT time is a string, such as Time(\"2024-01-01T00:00:00\", UTC(), ISOT()); " *
        "two numbers are the parts of a Julian date or a Modified Julian Date."))
    format === :mjd && return _time_jd(val1 + MJD_EPOCH, val2, scale, :mjd)
    return _time_jd(val1, val2, scale, format)
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
    return _time_jd(jd1, jd2, scale, format)
end

"""
    Time(val1::Real, val2::Real, s::AbstractTimeScale, f::AbstractTimeFormat) -> Time

Construct a time from two parts of a date in the format `f`, with typed scale and format tags:
Julian-date parts for `JD()`, Modified Julian Date parts for `MJD()`.
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
        return _time_jd(newjd1, newjd2, name, t.format)
    end

    # Fall back to normal field access
    error("Unknown property `$name` for Time object.")

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
    t1.scale === t2.scale || throw(ArgumentError(
        "Time: cannot subtract times in different scales ($(_scale_tag_str(t2.scale)) and " *
        "$(_scale_tag_str(t1.scale))); convert one first, for example `t2 - t1.$(t2.scale)`."))
    return (t2.jd1 - t1.jd1) + (t2.jd2 - t1.jd2)
end

"""
    -(t::Time{T}, dt::Real) -> Time{PT}

Subtract days `dt` from `t` returning new time object.
"""
function Base.:-(t::Time{T}, dt::Real) where {T<:Real}
    δ1, δ2 = _rebalance(-T(dt), zero(T))
    return _time_jd(getfield(t, :_jd1) + δ1, getfield(t, :_jd2) + δ2, t.scale, t.format)
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
    return _time_jd(jd1, jd2, t.scale, t.format)
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

# Format transforms, which read and write a Time, so they come after it.
include("formats/formats.jl")

end
