# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# A state that knows where and when it is.
#
# `axes_rotation` and `origin_translation` make you supply the epoch and the
# source frame on every call. That is right for the hot path and wrong for
# everything else — most of the time the thing being converted already knows
# both, and repeating them is how they end up mismatched.
#
# `Coordinate` is a state, an epoch, and a coordinate system composed into one
# object, so converting it takes a target and nothing else. It is what makes
# conversion work with no `Spacecraft` anywhere: an astronomer with a vector, a
# script, a test.
#
# The mechanism is the subject interface — `state_of`, `frame_of`, `epoch_of`.
# `Coordinate` implements it here; `Spacecraft` implements it in `AstroModels`.
# Every user-facing conversion takes either, through one code path, so the same
# call reads the same whether or not a spacecraft is involved.
#
# Acceleration is deliberately absent. A coordinate is kinematics. Acceleration
# is not a property of the subject at all — two identical spacecraft at
# identical states have different accelerations under different force models —
# so it arrives through `params`, pushed down by whoever computed it.
# =============================================================================

# --- The subject interface --------------------------------------------------
#
# Three questions any convertible thing must answer. Declared here because
# AstroFrames is what consumes them; implemented by whoever owns the type.

"""
    state_of(subject) -> AbstractOrbitState

The subject's orbital state, in whatever representation it holds.

# Arguments
- `subject`: a [`Coordinate`](@ref), a `Spacecraft`, or another type that
  implements the conversion contract.

# Returns
An `AstroStates` orbit state in the representation stored by the subject, such
as Cartesian, Keplerian, or equinoctial. AstroFrames converts the state to
Cartesian form when required.

# Notes
Types that support AstroFrames conversions implement this method together with
[`frame_of`](@ref) and [`epoch_of`](@ref).

# Example
```julia
AstroFrames.state_of(sc::MySatellite) = sc.orbit
AstroFrames.frame_of(sc::MySatellite) = sc.frame
AstroFrames.epoch_of(sc::MySatellite) = sc.epoch
```
"""
function state_of end

"""
    frame_of(subject) -> AbstractCoordinateSystem

The coordinate system the subject's state is expressed in.

# Arguments
- `subject` — see [`state_of`](@ref).

# Returns
The subject's current `CoordinateSystem` — origin and axes.

# Notes
This is what makes a conversion need only a target: the source comes from the
subject rather than from the caller, so the two cannot disagree.

# Example

```jldoctest
using AstroEpochs: Time, UTC, ISOT
frame_of(Coordinate([7000.0, 0, 1300, 0, 7.35, 1], EarthMJ2000Eq,
                    Time("2020-01-01T00:00:00.000", UTC(), ISOT())))

# output
CoordinateSystem:
  origin = Earth
  axes   = MJ2000Eq
```
"""
function frame_of end

"""
    epoch_of(subject) -> Time

The epoch the subject's state is valid at.

# Arguments
- `subject` — see [`state_of`](@ref).

# Returns
An `AstroEpochs.Time`, carrying its own scale.

# Notes
Every frame transform needs an epoch, and taking it from the subject is what
stops a state being converted at the wrong instant.

# Example

```jldoctest
using AstroEpochs: Time, UTC, ISOT
epoch_of(Coordinate([7000.0, 0, 1300, 0, 7.35, 1], EarthMJ2000Eq,
                    Time("2020-01-01T00:00:00.000", UTC(), ISOT()))) isa Time

# output
true
```
"""
function epoch_of end

# Guard the boundary: a type that has not implemented the interface must say so
# by name, not fail as a `MethodError` several calls down.
for f in (:state_of, :frame_of, :epoch_of)
    @eval function $f(subject)
        throw(ArgumentError(
            "$(typeof(subject)) cannot be converted between coordinate systems: it does " *
            "not implement the subject interface. A convertible type needs `state_of`, " *
            "`frame_of` and `epoch_of`. Wrap a bare state in a `Coordinate` " *
            "(`Coordinate(state, coord_sys, epoch)`), or add the three methods."))
    end
end

# That fallback matches every type, so `hasmethod` cannot answer "does this
# type implement the interface" — it is always true. Comparing against the
# fallback's own method can, and reading it from the fallback rather than
# keeping a second list means a type that implements `state_of` is recognised
# without also having to register itself somewhere.
const _SUBJECT_FALLBACK = which(state_of, Tuple{Any})

"""
    carries_own_state(x) -> Bool

Whether `x` carries its own state, as a spacecraft or `Coordinate` does, rather
than obtaining its position from an ephemeris.

Implementing `state_of` is what makes this true; there is nothing else to
register.
"""
carries_own_state(x) = which(state_of, Tuple{typeof(x)}) !== _SUBJECT_FALLBACK

"""
    _require_epoch_match(origin, jd_tdb, what)

An origin that carries its own state carries it at one epoch. Using it at any
other epoch is a silent error of whatever the origin moved in between, so it is
made loud instead.
"""
function _require_epoch_match(origin, jd_tdb::Real, what::AbstractString)
    origin_jd = _scales(epoch_of(origin)).tdb
    isapprox(origin_jd, jd_tdb; atol = 1e-9) && return nothing
    throw(ArgumentError(
        "$(what), but the origin's epoch (JD $(origin_jd) TDB) is not the epoch " *
        "being converted (JD $(jd_tdb) TDB) — a difference of " *
        "$((jd_tdb - origin_jd) * 86_400) s.
" *
        "
Propagate the origin to that epoch, or pass the reference orbit " *
        "yourself as `reference_state` in params."))
end

# --- Coordinate --------------------------------------------------------------

"""
    Coordinate(state, coord_sys, epoch)

`Coordinate` associates a state with the coordinate system in which it is
expressed and the epoch at which it is valid. A conversion therefore requires
only the target coordinate system.

# Fields
- `state` — an `AstroStates` orbit state. Any representation; Cartesian is the
  usual one. Position in km, velocity in km/s.
- `coord_sys` — the `CoordinateSystem` the state is expressed in.
- `time` — an `AstroEpochs.Time`, carrying its own scale.

# Notes
It composes rather than duplicating: it holds an `AstroStates` state, not raw
vectors, so Keplerian, equinoctial and every other representation come along
for free and there is no second state type in the codebase.

The frame is a field, not a type parameter. Putting it on the state types would
make one struct do two jobs, and small single-focus structs composed together
are the rule throughout Epicycle.

`Coordinate` carries no acceleration. Frames that require acceleration receive
it through `params`.

# Example
```julia
epoch = Time(2458849.5, 0.0, TDB(), JD())
c = Coordinate([7000.0, 0.0, 0.0, 0.0, 7.546, 0.0],
               CoordinateSystem(earth, MJ2000Eq()), epoch)

c_ec = Coordinate(c, CoordinateSystem(earth, MJ2000Ec()))   # re-express
st   = CartesianState(c, CoordinateSystem(earth, ITRF()))   # straight to a state
```
"""
struct Coordinate{S, CS<:AbstractCoordinateSystem, TT}
    state::S
    coord_sys::CS
    time::TT
end

"""
    Coordinate(v::AbstractVector, coord_sys, epoch)

Constructs a `Coordinate` from a six-element Cartesian state vector with
position in km followed by velocity in km/s. The constructor stores the vector
as a `CartesianState`.
"""
Coordinate(v::AbstractVector, cs::AbstractCoordinateSystem, epoch) =
    Coordinate(CartesianState(v), cs, epoch)

state_of(c::Coordinate) = c.state
frame_of(c::Coordinate) = c.coord_sys
epoch_of(c::Coordinate) = c.time

"""
    Base.show(io::IO, ::MIME"text/plain", c::Coordinate)

Show a coordinate as its frame, epoch, and state.
"""
function Base.show(io::IO, ::MIME"text/plain", c::Coordinate)
    # Kept to four lines. The nested `Time` and state have multi-line displays
    # of their own, and letting them through buries the three things a reader
    # actually wants: where, when, and what.
    println(io, "Coordinate:")
    println(io, "  origin = ", _origin_display(c.coord_sys.origin))
    println(io, "  axes   = ", nameof(typeof(c.coord_sys.axes)))
    println(io, "  epoch  = ", _epoch_summary(c.time))
    print(io,   "  state  = ", nameof(typeof(c.state)), " ",
                _state_summary(c.state))
end

# Always shown in TDB, whatever scale the epoch was built in. One scale in the
# display means two coordinates can be compared by eye; showing each in its own
# scale would make equal epochs look different.
_epoch_summary(t) = string(t.tdb.jd, " JD TDB")
_epoch_summary(t::Real) = string(t, " JD TDB")

_state_summary(s::AbstractOrbitState) = _format_elements(to_vector(s))
_state_summary(s::OrbitState)         = _format_elements(s.state)
_state_summary(s)                     = ""

_format_elements(v) =
    string("[", join((string(round(x; sigdigits = 8)) for x in v), ", "), "]")

Base.show(io::IO, c::Coordinate) = show(io, MIME"text/plain"(), c)

# --- Converting ---------------------------------------------------------------

"""
    _cartesian_vector(subject) -> AbstractVector

The subject's state as a plain Cartesian six-vector in its own frame.

Any representation is accepted; the gravitational parameter needed to leave a
Keplerian-like representation comes from the subject's own origin, which is
what makes this need no extra argument.
"""
_cartesian_vector(subject) = _as_cartesian(state_of(subject), frame_of(subject).origin)

# A subject may hold a concrete state or the tagged `OrbitState` container, and
# both have to work — a `Spacecraft` holds the latter.
_as_cartesian(s::CartesianState, origin) = to_vector(s)

_as_cartesian(s::AbstractOrbitState, origin) =
    to_vector(CartesianState(s, get_gravparam(origin)))

_as_cartesian(s::OrbitState, origin) =
    _as_cartesian(state_tag_to_type(s.statetype)(copy(s.state)), origin)

"""
    _reexpress(subject, target::AbstractCoordinateSystem, params) -> SVector{6}

The subject's state as a Cartesian six-vector in `target`.

Rotates the state into the target axes, then adds the origin offset expressed in
those axes. AstroFrames skips the translation when the origins match, which
keeps a pure rotation differentiable.
"""
function _reexpress(subject, target::AbstractCoordinateSystem, params::NamedTuple)
    source = frame_of(subject)
    epoch  = epoch_of(subject)
    v      = _cartesian_vector(subject)

    params = _reference_orbit_from_origin(target, epoch, params)

    out = axes_rotation(source.axes, target.axes, epoch, params) * SVector{6}(v)

    if source.origin !== target.origin
        out = out + origin_translation(source.origin, target.origin, target.axes, epoch, params)
    end
    return out
end

_reexpress(subject, target::AbstractCoordinateSystem) =
    _reexpress(subject, target, NamedTuple())

"""
    _reference_orbit_from_origin(target, epoch, params) -> NamedTuple

Supply `reference_state` from the target's origin when the origin *is* the
orbit the frame is defined by.

`CoordinateSystem(sc, RIC())` names the spacecraft as the origin, so its state
also defines the reference orbit. AstroFrames supplies that state automatically
unless `params` already contains `reference_state`.

Only the user tier does this. `axes_rotation` is the operator tier and takes
`params` exactly as given, so the raw call stays predictable and allocation
free.

An origin that is a body carries no state, so nothing is filled in and the
frame reports what it needs, as before.
"""
function _reference_orbit_from_origin(target::AbstractCoordinateSystem,
                                      epoch, params::NamedTuple)
    needs_reference_orbit(target.axes) || return params
    haskey(params, :reference_state)   && return params

    origin = target.origin
    carries_own_state(origin) || return params

    _require_epoch_match(origin, _scales(epoch).tdb,
        "$(nameof(typeof(target.axes))) axes take their reference orbit from the origin")

    reference = axes_rotation(frame_of(origin).axes, ICRF(), epoch) *
                SVector{6}(_cartesian_vector(origin))
    return merge(params, (; reference_state = reference))
end

"""
    Coordinate(subject, target::AbstractCoordinateSystem)
    Coordinate(subject, target::AbstractCoordinateSystem, params::NamedTuple)

The subject's coordinate, re-expressed in `target`.

`subject` may be a `Coordinate`, a `Spacecraft`, or another type that implements
`state_of`, `frame_of`, and `epoch_of`. The subject supplies its epoch and
source frame.

# Arguments
- `subject` — what to convert. See [`state_of`](@ref).
- `target` — the coordinate system to express it in.
- `params` — only for frames that need data this package cannot obtain on its
  own, such as the `reference_state` an orbit-relative frame is defined by.

# Returns
A `Coordinate` holding a `CartesianState` in `target`, at the subject's epoch.

# Notes
The conversion applies both the axes rotation and, when the origin differs, the
translation. Omitting the translation produces an error equal to the origin
separation, about 384,000 km from Earth to the Moon.

An origin change reads an ephemeris and is therefore not differentiable. A pure
rotation is.

# Example
```julia
c    = Coordinate(state, CoordinateSystem(earth, MJ2000Eq()), epoch)
c_ec = Coordinate(c, CoordinateSystem(earth, MJ2000Ec()))
c_mo = Coordinate(c, CoordinateSystem(moon,  MoonME()))

# a frame defined by another orbit needs that orbit
c_ric = Coordinate(c, CoordinateSystem(earth, RIC()), (; reference_state = chief))
```
"""
Coordinate(subject, target::AbstractCoordinateSystem, params::NamedTuple) =
    Coordinate(CartesianState(_reexpress(subject, target, params)), target, epoch_of(subject))

Coordinate(subject, target::AbstractCoordinateSystem) =
    Coordinate(subject, target, NamedTuple())

# --- The user tier, for callers with no Spacecraft ----------------------------

"""
    CartesianState(c::Coordinate)

The coordinate's state as a `CartesianState`, converting the representation if
it holds another one.
"""
CartesianState(c::Coordinate) = CartesianState(_cartesian_vector(c))

"""
    CartesianState(c::Coordinate, target::AbstractCoordinateSystem[, params])

The coordinate's state expressed in `target`, as a `CartesianState`.

The same call `AstroModels` provides for a `Spacecraft`, so a script that has a
state rather than a vehicle reads identically. See [`Coordinate`](@ref) for
what a conversion does.

# Example
```julia
CartesianState(c, CoordinateSystem(earth, ITRF()))
```
"""
CartesianState(c::Coordinate, target::AbstractCoordinateSystem, params::NamedTuple) =
    CartesianState(_reexpress(c, target, params))

CartesianState(c::Coordinate, target::AbstractCoordinateSystem) =
    CartesianState(c, target, NamedTuple())
