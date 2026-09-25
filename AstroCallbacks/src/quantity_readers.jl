# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Orbital quantities, as plain functions.
#
# Each takes a subject and, where the answer depends on one, a coordinate
# system: `inclination(sat, EarthMJ2000Ec)`. The dependency is positional,
# mapping one-for-one onto what a GMAT user writes as
# `sat.EarthMJ2000Ec.INC`. Omit it and the subject's own frame is used. A
# `Spacecraft` or `Coordinate` declares its frame through `frame_of`, so that
# reads a stated fact rather than guessing one.
#
# Where a quantity is an element of a state representation, the conversion
# happens in the function body. Nothing declares it: `element` takes the
# representation as an argument, so `semi_major_axis` names `KeplerianState`
# the same way it would name any other. That is what `AbstractCalc` carried as
# `calc_input_statetag(::SMA) = Keplerian()` and no longer needs to.
# =============================================================================

"""
    element(StateT, subject, cs, field) -> value
    element(StateT, subject, cs, field, params) -> value

One element of representation `StateT`, read in `cs`.

The gravitational parameter comes from the coordinate system's origin, so it is
not a separate argument.

# Examples
```julia
element(KeplerianState,         sat, EarthMJ2000Eq, :sma)
element(OutGoingAsymptoteState, sat, EarthMJ2000Eq, :c3)
```

# Returns
The requested field from `StateT`. Units follow the selected state field.
"""
function element(StateT, subject, cs::AbstractCoordinateSystem, field::Symbol,
                 params::NamedTuple)
    cart = CartesianState(Coordinate(subject, cs, params))
    return getfield(StateT(cart, get_gravparam(cs.origin)), field)
end

element(StateT, subject, cs::AbstractCoordinateSystem, field::Symbol) =
    element(StateT, subject, cs, field, NamedTuple())

element(StateT, subject, field::Symbol) =
    element(StateT, subject, frame_of(subject), field, NamedTuple())

"""
    set_element!(StateT, subject, cs, field[, params]; to) -> subject

Write one element of representation `StateT`, measured in `cs`, holding the rest of the
representation. The write counterpart of [`element`](@ref), and the one place a setter
converts: the subject's state is read into `StateT` in `cs`, `field` is replaced by `to`, and the
result is converted back and stored through [`set_state!`](@ref) in the subject's own frame.

# Arguments
- `StateT`: the representation the element belongs to, such as `KeplerianState`.
- `subject`: a subject with a `set_state!` method.
- `cs`: the coordinate system the element is measured in.
- `field`: the element's field name in `StateT`, such as `:sma`.
- `params`: data for an orbit-relative coordinate system, as the readers take it.
- `to`: the new value, in the units of the field.

# Returns
`subject`, mutated. Throws `ArgumentError` when the new elements do not define an orbit, and
leaves the subject unchanged.
"""
function set_element!(StateT, subject, cs::AbstractCoordinateSystem, field::Symbol,
                      params::NamedTuple = NamedTuple(); to)
    mu   = get_gravparam(cs.origin)
    held = StateT(CartesianState(Coordinate(subject, cs, params)), mu)
    vals = [f === field ? to : getfield(held, f) for f in fieldnames(typeof(held))]
    new  = StateT(promote(vals...)...)
    _check_elements(new)
    cart = CartesianState(new, mu)
    all(isfinite, to_vector(cart)) || throw(ArgumentError(
        "$(nameof(StateT)) with $field = $to does not define an orbit; the state was not " *
        "changed."))
    return set_state!(subject, cart, cs, params)
end

# The conditions under which Keplerian elements describe an orbit. Checked before converting,
# because the conversion does not refuse them: some give NaN with a logged warning, and the rest
# give a state that reads back as different elements from the ones written.
_check_elements(::Any) = nothing

function _check_elements(b::BrouwerMeanLongState)
    0 <= b.ecc < 0.99 || throw(ArgumentError(
        "mean eccentricity must be at least 0 and below 0.99 for the Brouwer-Lyddane theory; got " *
        "$(b.ecc). The state was not changed."))
    b.sma * (1 - b.ecc) > 3000 || throw(ArgumentError(
        "mean periapsis radius $(b.sma * (1 - b.ecc)) km is below the 3000 km the " *
        "Brouwer-Lyddane theory holds for. The state was not changed."))
    return nothing
end

function _check_elements(s::SphericalRADECState)
    s.r > 0 || throw(ArgumentError(
        "a position magnitude must be positive; got $(s.r) km. The state was not changed."))
    s.v >= 0 || throw(ArgumentError(
        "a velocity magnitude must be zero or more; got $(s.v) km/s. The state was not changed."))
    return nothing
end

function _check_elements(k::KeplerianState)
    a, e, i, ta = k.sma, k.ecc, k.inc, k.ta
    e >= 0 || throw(ArgumentError(
        "eccentricity must be zero or more; got $e. The state was not changed."))
    a * (1 - e) > 0 || throw(ArgumentError(
        "semi-major axis $a km and eccentricity $e do not describe an orbit: an ellipse " *
        "(e < 1) needs a positive semi-major axis and a hyperbola (e > 1) a negative one, and " *
        "e = 1 has no semi-major axis. The state was not changed."))
    0 <= i <= π || throw(ArgumentError(
        "inclination must be between 0 and π rad; got $i. The state was not changed."))
    if e > 1
        limit = acos(-1 / e)
        abs(rem2pi(ta, RoundNearest)) < limit || throw(ArgumentError(
            "true anomaly $ta rad is beyond the asymptote of a hyperbola with eccentricity " *
            "$e, which reaches ±$(round(limit, digits = 4)) rad. The state was not changed."))
    end
    return nothing
end

"""
    set_state!(subject, cart::CartesianState, cs[, params]) -> subject

Store a Cartesian state given in `cs` on `subject`, in the frame and form the subject keeps. A
subject becomes writable through every element setter by having a method of this function.

# Returns
`subject`, mutated.
"""
function set_state!(sc::Spacecraft, cart::CartesianState, cs::AbstractCoordinateSystem,
                    params::NamedTuple = NamedTuple())
    own = Coordinate(Coordinate(to_vector(cart), cs, epoch_of(sc)), frame_of(sc), params)
    set_posvel!(sc, to_vector(CartesianState(own)))
    return sc
end

# --- Keplerian elements ------------------------------------------------------
#
# Written out rather than generated, so each carries its own docstring and
# `?semi_major_axis` says something useful.

for (fn, field, what) in
        ((:semi_major_axis,      :sma,  "Semi-major axis, in km."),
         (:eccentricity,         :ecc,  "Eccentricity."),
         (:inclination,          :inc,  "Inclination, in radians."),
         (:argument_of_periapsis,:aop,  "Argument of periapsis, in radians."),
         (:true_anomaly,         :ta,   "True anomaly, in radians."))
    @eval begin
        $fn(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
            element(KeplerianState, subject, cs, $(QuoteNode(field)), params)
        $fn(subject, cs::AbstractCoordinateSystem) =
            element(KeplerianState, subject, cs, $(QuoteNode(field)))
        $fn(subject) =
            element(KeplerianState, subject, frame_of(subject), $(QuoteNode(field)))
    end
end

# --- Keplerian setters -------------------------------------------------------
#
# Each holds the other five elements, in the coordinate system it is given. On an equatorial orbit
# the RAAN is undefined, and on a circular one so is the argument of periapsis; setting either
# there leaves the orbit where it was, because the state has no place to put it.

for (fn, field, what, unit) in
        ((:semi_major_axis,       :sma,  "semi-major axis",          " km"),
         (:eccentricity,          :ecc,  "eccentricity",             ""),
         (:inclination,           :inc,  "inclination",              " rad"),
         (:raan,                  :raan, "right ascension of the ascending node", " rad"),
         (:argument_of_periapsis, :aop,  "argument of periapsis",    " rad"),
         (:true_anomaly,          :ta,   "true anomaly",             " rad"))
    fn! = Symbol(fn, "!")
    doc = """
        $(fn!)(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to`$unit, measured in `cs` or the spacecraft's own frame. The other
    Keplerian elements are held in that frame, and the state is stored in the spacecraft's own
    frame as before.

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when the elements that
    result do not describe an orbit.

    # Example
    ```julia
    $(fn!)(sat; to = $(fn === :semi_major_axis ? "7200.0" : fn === :eccentricity ? "0.01" : "0.5"))
    $(fn!)(sat, EarthMJ2000Ec; to = $(fn === :semi_major_axis ? "7200.0" : fn === :eccentricity ? "0.01" : "0.5"))
    ```
    """
    @eval begin
        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::Real)
            return set_element!(KeplerianState, sat, cs, $(QuoteNode(field)), params; to = to)
        end
        @doc $doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)
    end
end

"""
    semi_major_axis(subject[, cs][, params])

Semi-major axis in km, measured in `cs` or the subject's frame.

# Returns
The semi-major axis in km.

# Example
```julia
semi_major_axis(sat, EarthMJ2000Eq)
```
""" semi_major_axis

"""
    eccentricity(subject[, cs][, params])

Eccentricity measured in `cs` or the subject's frame.

# Returns
The dimensionless eccentricity.

# Example
```julia
eccentricity(sat)
```
""" eccentricity

"""
    inclination(subject[, cs][, params])

Inclination in radians, measured in `cs` or the subject's frame.

# Returns
The inclination in radians.

# Example
```julia
inclination(sat, EarthMJ2000Eq)
```
""" inclination

"""
    argument_of_periapsis(subject[, cs][, params])

Argument of periapsis in radians, measured in `cs` or the subject's frame.

# Returns
The argument of periapsis in radians.

# Example
```julia
argument_of_periapsis(sat)
```
""" argument_of_periapsis

"""
    true_anomaly(subject[, cs][, params])

True anomaly in radians, measured in `cs` or the subject's frame.

# Returns
The true anomaly in radians.

# Example
```julia
true_anomaly(sat)
```
""" true_anomaly

# `raan` already exists in frame_aware_quantities.jl and is left alone.

"""
    mean_long_sma(subject[, cs][, params])

Brouwer-Lyddane long-period mean semi-major axis, in km. Settable with
[`mean_long_sma!`](@ref).

The mean element averages out the short-period variation an osculating
semi-major axis carries, so a station-keeping box stated in mean elements does
not chase the once-per-orbit oscillation. The conversion runs the mean-element
theory and carries its own convergence tolerance, so a round trip agrees to
about 1e-4 km rather than to machine precision.

Measured in `cs`; without one, in the subject's own frame. See [`element`](@ref).

# Returns
The mean semi-major axis in km.

# Example
```julia
mean_long_sma(sat)
mean_long_sma(sat, EarthMJ2000Eq)
```
"""
mean_long_sma(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    element(BrouwerMeanLongState, subject, cs, :sma, params)
mean_long_sma(subject, cs::AbstractCoordinateSystem) =
    element(BrouwerMeanLongState, subject, cs, :sma)
mean_long_sma(subject) =
    element(BrouwerMeanLongState, subject, frame_of(subject), :sma)

"""
    mean_long_sma!(sat[, cs][, params]; to) -> sat

Set the Brouwer-Lyddane long-period mean semi-major axis of `sat` to `to` km, measured in `cs`
or the spacecraft's own frame. The other mean elements are held in that frame, the mean anomaly
among them, and the state is stored in the spacecraft's own frame as before.

# Notes
The conversion back from mean elements iterates to a tolerance, so the value read back agrees
with `to` to about 1e-4 km. The theory holds for `0 ≤ e < 0.99` and a mean periapsis radius above
3000 km; outside that the setter throws.

# Returns
`sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when the mean elements that
result are outside the theory's domain.

# Example
```julia
mean_long_sma!(sat; to = 7000.0)
```
"""
mean_long_sma!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
               params::NamedTuple = NamedTuple(); to::Real) =
    set_element!(BrouwerMeanLongState, sat, cs, :sma, params; to = to)

EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof(mean_long_sma),
                           cs::AbstractCoordinateSystem = frame_of(sat),
                           params::NamedTuple = NamedTuple(); to) =
    mean_long_sma!(sat, cs, params; to = to)

"""
    outgoing_rla(subject[, cs][, params])

Right ascension of the outgoing hyperbolic asymptote, in radians.

Defined on a hyperbolic orbit, where it names the direction the spacecraft
departs along. Reading it from a bound orbit converts through a state that has
no outgoing asymptote and the result means nothing.

Measured in `cs`; without one, in the subject's own frame. See [`element`](@ref).

# Returns
The outgoing asymptote right ascension in radians.

# Example
```julia
outgoing_rla(sat, EarthMJ2000Eq)
```
"""
outgoing_rla(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    element(OutGoingAsymptoteState, subject, cs, :rla, params)
outgoing_rla(subject, cs::AbstractCoordinateSystem) =
    element(OutGoingAsymptoteState, subject, cs, :rla)
outgoing_rla(subject) =
    element(OutGoingAsymptoteState, subject, frame_of(subject), :rla)

# --- position and velocity components ---------------------------------------

for (fn, idx, what) in ((:position_x, 1, "X"), (:position_y, 2, "Y"), (:position_z, 3, "Z"))
    @eval begin
        $fn(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
            position_vector(subject, cs, params)[$idx]
        $fn(subject, cs::AbstractCoordinateSystem) = position_vector(subject, cs)[$idx]
        $fn(subject) = position_vector(subject)[$idx]
    end
end

"""
    position_x(subject[, cs][, params])

X position component in km, expressed in `cs` or the subject's frame.

# Returns
The X position component in km.

# Example
```julia
position_x(sat, EarthMJ2000Eq)
```
""" position_x

"""
    position_y(subject[, cs][, params])

Y position component in km, expressed in `cs` or the subject's frame.

# Returns
The Y position component in km.

# Example
```julia
position_y(sat, EarthMJ2000Eq)
```
""" position_y

"""
    position_z(subject[, cs][, params])

Z position component in km, expressed in `cs` or the subject's frame.

# Returns
The Z position component in km.

# Example
```julia
position_z(sat, EarthMJ2000Eq)
```
""" position_z

# --- Cartesian setters -------------------------------------------------------
#
# A Cartesian component or vector is written in the coordinate system given, holding the other
# components there, and the state is stored in the spacecraft's own frame as before.

# Replace components `idx` of the six-vector [x, y, z, vx, vy, vz] as read in `cs`.
function _set_cartesian!(sat::Spacecraft, cs::AbstractCoordinateSystem, params::NamedTuple,
                         idx, to)
    x = collect(to_vector(CartesianState(Coordinate(sat, cs, params))))
    length(to) == length(idx) || throw(ArgumentError(
        "expected $(length(idx)) component$(length(idx) == 1 ? "" : "s"); got $(length(to))."))
    x = convert(Vector{promote_type(eltype(x), eltype(to))}, x)
    x[idx] .= to
    all(isfinite, x) || throw(ArgumentError(
        "a Cartesian state must be finite; got $x. The state was not changed."))
    return set_state!(sat, CartesianState(x), cs, params)
end

for (fn, idx, what, unit) in
        ((:position_x,      1:1, "X component of position", "km"),
         (:position_y,      2:2, "Y component of position", "km"),
         (:position_z,      3:3, "Z component of position", "km"),
         (:position_vector, 1:3, "position vector",         "km"),
         (:velocity_vector, 4:6, "velocity vector",         "km/s"))
    fn! = Symbol(fn, "!")
    scalar = length(idx) == 1
    example = scalar ? "7000.0" : (fn === :position_vector ? "[7000.0, 0.0, 0.0]" :
                                                             "[0.0, 7.5, 0.0]")
    doc = """
        $(fn!)(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to`, in $unit, measured in `cs` or the spacecraft's own frame. The
    other Cartesian components are held in that frame, and the state is stored in the
    spacecraft's own frame as before.

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when `to` has the wrong
    number of components or is not finite.

    # Example
    ```julia
    $(fn!)(sat; to = $example)
    $(fn!)(sat, EarthMJ2000Ec; to = $example)
    ```
    """
    T = scalar ? :Real : :(AbstractVector{<:Real})
    wrap = scalar ? :([to]) : :to
    @eval begin
        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::$T)
            return _set_cartesian!(sat, cs, params, $idx, $wrap)
        end
        @doc $doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)
    end
end

"""
    position_magnitude(subject[, cs][, params])

Distance from the coordinate system's origin, in km.

# Returns
The position magnitude in km.

# Example
```julia
position_magnitude(sat)
position_magnitude(sat, EarthMJ2000Eq)
```
"""
position_magnitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    norm(position_vector(subject, cs, params))
position_magnitude(subject, cs::AbstractCoordinateSystem) = norm(position_vector(subject, cs))
position_magnitude(subject) = norm(position_vector(subject))

"""
    velocity_magnitude(subject[, cs][, params])

Speed relative to the coordinate system, in km/s.

Frame-dependent in a stronger sense than position: a rotating target frame
contributes a term of size `ω × r`, about 0.47 km/s at Earth's equator.

# Returns
The velocity magnitude in km/s.

# Example
```julia
velocity_magnitude(sat)
velocity_magnitude(sat, EarthFixed)
```
"""
velocity_magnitude(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    norm(velocity_vector(subject, cs, params))
velocity_magnitude(subject, cs::AbstractCoordinateSystem) = norm(velocity_vector(subject, cs))
velocity_magnitude(subject) = norm(velocity_vector(subject))

# --- magnitude setters --------------------------------------------------------
#
# A magnitude is the r or v of the spherical state, so setting it holds that vector's direction
# and leaves the other vector alone. A zero vector has no direction to hold, so it is refused.

for (fn, field, what, unit) in
        ((:position_magnitude, :r, "distance from the coordinate system's origin", "km"),
         (:velocity_magnitude, :v, "speed",                                        "km/s"))
    fn! = Symbol(fn, "!")
    doc = """
        $(fn!)(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to` $unit, measured in `cs` or the spacecraft's own frame. The
    vector's direction is held, and so is the other vector; the state is stored in the
    spacecraft's own frame as before.

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when `to` is negative or
    the vector is zero, which has no direction to hold.

    # Example
    ```julia
    $(fn!)(sat; to = $(fn === :position_magnitude ? "7200.0" : "7.6"))
    ```
    """
    @eval begin
        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::Real)
            $fn(sat, cs, params) > 0 || throw(ArgumentError(
                $("the $what is zero, so there is no direction to hold; set the vector instead.")))
            return set_element!(SphericalRADECState, sat, cs, $(QuoteNode(field)), params; to = to)
        end
        @doc $doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)
    end
end

"""
    position_dot_velocity(subject[, cs][, params])

The dot product of position and velocity, in km²/s.

Zero at periapsis and apoapsis, and the standard way to stop on either: it
rises through zero at periapsis and falls through zero at apoapsis, so
`direction` picks which. Unlike true anomaly it does not wrap, so a
root-finder brackets it directly.

# Returns
The position-velocity dot product in km^2/s.

# Example
```julia
position_dot_velocity(sat)
```
"""
position_dot_velocity(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
    dot(position_vector(subject, cs, params), velocity_vector(subject, cs, params))
position_dot_velocity(subject, cs::AbstractCoordinateSystem) =
    dot(position_vector(subject, cs), velocity_vector(subject, cs))
position_dot_velocity(subject) = dot(position_vector(subject), velocity_vector(subject))

"""
    epoch(subject) -> Time

When the subject's state is valid. Takes no coordinate system.

# Returns
The subject's [`Time`](@ref) epoch.

# Example
```julia
epoch(sat)
```
"""
epoch(subject) = epoch_of(subject)

"""
    epoch!(spacecraft; to) -> spacecraft

Move the spacecraft's epoch to `to`, keeping its state's numbers in its own frame. It does not
propagate: the spacecraft is at the same position and velocity at the new epoch. This is what lets
a solver vary a date, such as a launch or arrival epoch.

# Returns
`spacecraft`, mutated.

# Example
```julia
epoch!(sat; to = Time("2024-02-01T00:00:00", UTC(), ISOT()))
```
"""
epoch!(sc::Spacecraft; to::Time) = (sc.time = to; sc)

EpicycleBase.set_quantity!(sc::Spacecraft, ::typeof(epoch); to) = epoch!(sc; to = to)

# --- maneuver quantities -----------------------------------------------------

"""
    delta_v(maneuver) -> [dv1, dv2, dv3]

The maneuver's delta-V, in km/s, in the axes it was built with.

Settable with [`delta_v!`](@ref), which is also what lets a solver vary it:
`Vary(delta_v, toi; ...)` reads through this function and writes through that one. The
spacecraft is not an argument: a delta-V belongs to the maneuver, and applying it to a
spacecraft is what `maneuver!` does.

# Returns
A three-element delta-V vector in km/s.

# Example
```julia
toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1)
delta_v(toi)
```
"""
delta_v(man) = [man.element1, man.element2, man.element3]

"""
    delta_v!(maneuver; to) -> maneuver

Write a three-component delta-V, in km/s, onto a maneuver, in the axes it was built with.

# Arguments
- `maneuver`: the maneuver to write to.
- `to`: three delta-V components, in km/s, in the maneuver's own axes.

# Returns
`maneuver`, mutated. Throws `ArgumentError` when `to` does not have three components, and leaves
the maneuver unchanged.

# Example
```julia
toi = ImpulsiveManeuver(axes = VNB())
delta_v!(toi; to = [0.1, 0.0, 0.0])
```
"""
function delta_v!(man::ImpulsiveManeuver; to::AbstractVector{<:Real})
    length(to) == 3 || throw(ArgumentError(
        "delta_v must have 3 components; got $(length(to))."))
    setfield!(man, :element1, oftype(man.element1, to[1]))
    setfield!(man, :element2, oftype(man.element2, to[2]))
    setfield!(man, :element3, oftype(man.element3, to[3]))
    return man
end

EpicycleBase.set_quantity!(man::ImpulsiveManeuver, ::typeof(delta_v); to) =
    delta_v!(man; to = to)

"""
    delta_v_magnitude(maneuver)

The size of the maneuver's delta-V, in km/s.

Settable with [`delta_v_magnitude!`](@ref), which scales the delta-V and holds its direction.

# Returns
The delta-V magnitude in km/s.

# Example
```julia
toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1)
delta_v_magnitude(toi)
```
"""
delta_v_magnitude(man) = norm(delta_v(man))

"""
    delta_v_magnitude!(maneuver; to) -> maneuver

Scale the maneuver's delta-V to `to` km/s, holding its direction.

# Returns
`maneuver`, mutated. Throws `ArgumentError`, and leaves the maneuver unchanged, when `to` is
negative or the delta-V is zero, which has no direction to hold.

# Example
```julia
toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1)
delta_v_magnitude!(toi; to = 0.25)
```
"""
function delta_v_magnitude!(man::ImpulsiveManeuver; to::Real)
    to >= 0 || throw(ArgumentError(
        "a delta-V magnitude must be zero or more; got $to km/s."))
    m = delta_v_magnitude(man)
    m > 0 || throw(ArgumentError(
        "the maneuver's delta-V is zero, so there is no direction to hold; set `delta_v!` " *
        "instead."))
    return delta_v!(man; to = delta_v(man) .* (to / m))
end

EpicycleBase.set_quantity!(man::ImpulsiveManeuver, ::typeof(delta_v_magnitude); to) =
    delta_v_magnitude!(man; to = to)

"""
    state(spacecraft)

The spacecraft's position and velocity, as a six-vector in its own coordinate
system.

Settable with [`state!`](@ref), which is also what lets a solver estimate it: `Vary(state,
sat; guess = [...], covariance = ...)` reads through this function and writes through that one.
An orbit determination problem varies this the same way a targeting problem varies
[`delta_v`](@ref).

The name is shared with the state block of a transcribed phase, which is the
same idea about a different subject: what the solver is free to move.

# Returns
A six-element Cartesian state in km and km/s, in the spacecraft's coordinate system.

# Example
```julia
state(sat)
```
"""
state(sc::Spacecraft) = to_posvel(sc)

"""
    state!(spacecraft; to) -> spacecraft

Replace the spacecraft's position and velocity with a six-vector in km and km/s, in its own
coordinate system.

# Returns
`spacecraft`, mutated.

# Example
```julia
state!(sat; to = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
```
"""
state!(sc::Spacecraft; to::AbstractVector{<:Real}) = (set_posvel!(sc, to); sc)

EpicycleBase.set_quantity!(sc::Spacecraft, ::typeof(state); to) = state!(sc; to = to)

# --- body quantities ---------------------------------------------------------

"""
    gravitational_parameter(body)

The body's gravitational parameter mu, in km^3/s^2. Takes no coordinate system.

Settable with [`gravitational_parameter!`](@ref), which is also what lets a solver estimate
it: `Vary(gravitational_parameter, earth; ...)` reads through this function and writes through
that one.

# Returns
The gravitational parameter in km^3/s^2.

# Example
```jldoctest
gravitational_parameter(earth)

# output

398600.4418
```
"""
gravitational_parameter(body) = AstroUniverse.get_gravparam(body)

"""
    gravitational_parameter!(body; to) -> body

Set the body's gravitational parameter mu, in km^3/s^2.

# Returns
`body`, mutated.

# Example
```julia
moon_copy = deepcopy(moon)
gravitational_parameter!(moon_copy; to = 4902.8)
```
"""
gravitational_parameter!(body::CelestialBody; to::Real) =
    (AstroUniverse.set_gravparam!(body, to); body)

EpicycleBase.set_quantity!(body::CelestialBody, ::typeof(gravitational_parameter); to) =
    gravitational_parameter!(body; to = to)

# --- output partials ---------------------------------------------------------
#
# d(quantity)/d(state), one row per component of the quantity and six columns,
# for the quantities that are functions of the Cartesian state and of nothing
# else. Declared here beside the reader each one differentiates.
#
# These are the derivatives worth writing by hand. Each is a few lines, each is
# checked against a finite difference in the tests, and each saves the solver a
# forward-mode pass over a reader that would otherwise be differentiated on
# every iteration.
#
# The Keplerian elements and the asymptote quantities are not here. Their
# closed forms are long enough to be worth getting wrong, so they fall to
# automatic differentiation until someone needs the speed. That is the correct
# answer, only slower, which is what the fallback is for.
#
# Every method here is the subject's-own-frame form. Read in a named frame, a
# quantity is a different function of the state and needs its own method.

_pv(subject) = (position_vector(subject), velocity_vector(subject))

"""d(r)/d(x) — position is the first half of the state."""
EpicycleBase.output_partial(subject, ::typeof(position_vector)) =
    hcat(Matrix{Float64}(I, 3, 3), zeros(3, 3))

"""d(v)/d(x) — velocity is the second half."""
EpicycleBase.output_partial(subject, ::typeof(velocity_vector)) =
    hcat(zeros(3, 3), Matrix{Float64}(I, 3, 3))

for (fn, idx) in ((:position_x, 1), (:position_y, 2), (:position_z, 3))
    @eval EpicycleBase.output_partial(subject, ::typeof($fn)) =
        reshape(Float64[i == $idx for i in 1:6], 1, 6)
end

"""d(|r|)/d(x) = [r' / |r|, 0]."""
function EpicycleBase.output_partial(subject, ::typeof(position_magnitude))
    r = position_vector(subject)
    return hcat(transpose(r ./ norm(r)), zeros(1, 3))
end

"""d(|v|)/d(x) = [0, v' / |v|]."""
function EpicycleBase.output_partial(subject, ::typeof(velocity_magnitude))
    v = velocity_vector(subject)
    return hcat(zeros(1, 3), transpose(v ./ norm(v)))
end

"""d(r.v)/d(x) = [v', r']."""
function EpicycleBase.output_partial(subject, ::typeof(position_dot_velocity))
    r, v = _pv(subject)
    return hcat(transpose(v), transpose(r))
end

"""
d(a)/d(x), from the vis-viva relation a = 1 / (2/|r| - |v|^2 / mu).

Differentiating that gives d(a)/d(r) = 2 a^2 r / |r|^3 and
d(a)/d(v) = 2 a^2 v / mu, so mu is the only thing needed beyond the state.
"""
function EpicycleBase.output_partial(subject, ::typeof(semi_major_axis))
    r, v  = _pv(subject)
    mu    = get_gravparam(frame_of(subject).origin)
    rn    = norm(r)
    a     = 1 / (2 / rn - dot(v, v) / mu)
    twoa2 = 2 * a^2
    return hcat((twoa2 / rn^3) .* transpose(r), (twoa2 / mu) .* transpose(v))
end

# --- traits ------------------------------------------------------------------
#
# What each quantity carries beyond its value. Declaring these is the whole
# opt-in; see `EpicycleBase/src/quantity_traits.jl`.

EpicycleBase.label(::typeof(semi_major_axis))      = "Semi-major axis"
EpicycleBase.label(::typeof(eccentricity))         = "Eccentricity"
EpicycleBase.label(::typeof(inclination))          = "Inclination"
EpicycleBase.label(::typeof(raan))                 = "RAAN"
EpicycleBase.label(::typeof(argument_of_periapsis))= "Argument of periapsis"
EpicycleBase.label(::typeof(true_anomaly))         = "True anomaly"
EpicycleBase.label(::typeof(position_vector))      = "Position"
EpicycleBase.label(::typeof(velocity_vector))      = "Velocity"
EpicycleBase.label(::typeof(position_x))           = "X"
EpicycleBase.label(::typeof(position_y))           = "Y"
EpicycleBase.label(::typeof(position_z))           = "Z"
EpicycleBase.label(::typeof(position_magnitude))   = "Position magnitude"
EpicycleBase.label(::typeof(velocity_magnitude))   = "Velocity magnitude"
EpicycleBase.label(::typeof(epoch))                = "Epoch"
EpicycleBase.label(::typeof(state))                = "State"
EpicycleBase.label(::typeof(delta_v))              = "Delta-V"
EpicycleBase.label(::typeof(position_dot_velocity)) = "r-dot-v"
EpicycleBase.label(::typeof(mean_long_sma))        = "Long-period mean semi-major axis"
EpicycleBase.label(::typeof(outgoing_rla))         = "Outgoing RLA"
EpicycleBase.label(::typeof(delta_v_magnitude))    = "Delta-V magnitude"
EpicycleBase.label(::typeof(gravitational_parameter)) = "Gravitational parameter"

# Angles wrap. A stopping condition has to unwrap before it brackets a root;
# nothing else needs to know.
for f in (:inclination, :raan, :argument_of_periapsis, :true_anomaly, :outgoing_rla)
    @eval EpicycleBase.is_cyclic(::typeof($f)) = true
    @eval EpicycleBase.cycle(::typeof($f))     = 2π
end

# --- finding quantities ------------------------------------------------------

"""
    subject_type(quantity) -> Type

The type of subject `quantity` reads, which is how [`quantities`](@ref) finds the quantities for a
spacecraft, a maneuver or a body. `Any` by default.

A custom quantity declares its subject type alongside its label:
`AstroCallbacks.subject_type(::typeof(drag_area)) = Spacecraft`.

# Returns
A type; a `Union` when the quantity reads more than one kind of subject.

# Example
```julia
AstroCallbacks.subject_type(delta_v)          # ImpulsiveManeuver
```
"""
subject_type(::Any) = Any

for f in (:position_vector, :velocity_vector, :position_x, :position_y, :position_z,
          :position_magnitude, :velocity_magnitude, :position_dot_velocity,
          :semi_major_axis, :eccentricity, :inclination, :raan, :argument_of_periapsis,
          :true_anomaly, :mean_long_sma, :outgoing_rla, :epoch)
    @eval subject_type(::typeof($f)) = Union{Spacecraft, Coordinate}
end
subject_type(::typeof(state))                   = Spacecraft
subject_type(::typeof(delta_v))                 = ImpulsiveManeuver
subject_type(::typeof(delta_v_magnitude))       = ImpulsiveManeuver
subject_type(::typeof(gravitational_parameter)) = CelestialBody

"""
    quantities([T], [mod]) -> Vector{Function}

The quantities that read a subject of type `T`, in alphabetical order.

# Arguments
- `T`: a subject type such as `Spacecraft`, `ImpulsiveManeuver` or `CelestialBody`. Without it,
  every quantity is returned.
- `mod`: the module to look in, `AstroCallbacks` by default. `Main` finds the quantities a script
  has defined.

# Notes
A function is a quantity when it declares a label, and it reads a `T` when `T` is a subtype of its
[`subject_type`](@ref). A quantity that declares no subject type appears only when `T` is omitted.

# Returns
The quantity functions themselves, so a result can be passed to `StopAt` or `Constraint`.

# Example
```julia
quantities(ImpulsiveManeuver)     # [delta_v, delta_v_magnitude]
quantities(Spacecraft)
```
"""
function quantities(T::Type = Any, mod::Module = AstroCallbacks)
    found = Function[]
    for n in names(mod; all = true)
        startswith(String(n), '#') && continue
        isdefined(mod, n) || continue
        f = getfield(mod, n)
        f isa Function && EpicycleBase.label(f) != "quantity" || continue
        S = subject_type(f)
        (T === Any || (S !== Any && T <: S)) && push!(found, f)
    end
    return sort!(unique!(found); by = f -> String(nameof(f)))
end
