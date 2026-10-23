# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# One quantity per element of a state representation (#11, #21).
#
# Every field of every representation in AstroStates is readable by name, and settable by the
# same name with `!`. A field that more than one representation carries with the same meaning,
# such as the semi-major axis or the equinoctial h, is one quantity; a field whose name is reused
# for a different value, such as the h of the modified equinoctial set, gets its own. The Keplerian
# and Cartesian quantities, and the ones another representation stands behind (apsis radii,
# semi-latus rectum, true longitude, C3), are defined elsewhere; `_FIELD_QUANTITY` below names
# them all, so a setter's documentation can say what it holds.
#
# A setter writes its element through `set_element!`: the state is read into the representation
# in the coordinate system given, the one element is replaced, and the result is stored in the
# spacecraft's own frame. The rest of that representation is held.
# =============================================================================

using AstroStates: SphericalRADECState, SphericalAZIFPAState, EquinoctialState,
                   AlternateEquinoctialState, ModifiedEquinoctialState, DelaunayState,
                   BrouwerMeanLongState, BrouwerMeanShortState, OutGoingAsymptoteState,
                   IncomingAsymptoteState

# The quantity that reads each field of each representation. A setter's documentation lists the
# fields it holds by these names.
const _FIELD_QUANTITY = Dict{Tuple{Symbol, Symbol}, Symbol}(
    (:CartesianState, :position) => :position_vector,
    (:CartesianState, :velocity) => :velocity_vector,
    (:KeplerianState, :sma) => :semi_major_axis, (:KeplerianState, :ecc) => :eccentricity,
    (:KeplerianState, :inc) => :inclination, (:KeplerianState, :raan) => :raan,
    (:KeplerianState, :aop) => :argument_of_periapsis, (:KeplerianState, :ta) => :true_anomaly,
    (:ModifiedKeplerianState, :rp) => :periapsis_radius,
    (:ModifiedKeplerianState, :ra) => :apoapsis_radius,
    (:ModifiedKeplerianState, :inc) => :inclination, (:ModifiedKeplerianState, :raan) => :raan,
    (:ModifiedKeplerianState, :aop) => :argument_of_periapsis,
    (:ModifiedKeplerianState, :ta) => :true_anomaly,
    (:ModifiedEquinoctialState, :p) => :semilatus_rectum,
    (:ModifiedEquinoctialState, :L) => :true_longitude)
for S in (:SphericalRADECState, :SphericalAZIFPAState)
    _FIELD_QUANTITY[(S, :r)]   = :position_magnitude
    _FIELD_QUANTITY[(S, :v)]   = :velocity_magnitude
    _FIELD_QUANTITY[(S, :ra)]  = :right_ascension
    _FIELD_QUANTITY[(S, :dec)] = :declination
end
for S in (:EquinoctialState, :AlternateEquinoctialState)
    _FIELD_QUANTITY[(S, :a)]     = :semi_major_axis
    _FIELD_QUANTITY[(S, :h)]     = :equinoctial_h
    _FIELD_QUANTITY[(S, :k)]     = :equinoctial_k
    _FIELD_QUANTITY[(S, :mlong)] = :mean_longitude
end
for S in (:OutGoingAsymptoteState, :IncomingAsymptoteState)
    _FIELD_QUANTITY[(S, :rp)] = :periapsis_radius
    _FIELD_QUANTITY[(S, :c3)] = :c3
    _FIELD_QUANTITY[(S, :ta)] = :true_anomaly
end

"""
    _read_field(subject, S, cs, params, field, what)

Field `field` of representation `S` of `subject` in `cs`. Throws an `ArgumentError` naming `what`
when the conversion to `S` fails for this orbit, rather than returning the `NaN` it gives.
"""
function _read_field(subject, S, cs::AbstractCoordinateSystem, params::NamedTuple, field::Symbol,
                     what::String)
    value = getfield(_state(S, subject, cs, params), field)
    isnan(value) && throw(ArgumentError(
        "$what is undefined for this orbit: it has no $(nameof(S)) (see the warning logged by the " *
        "conversion)."))
    return value
end

# fn, representation, field, label, unit, cyclic, what, definition, notes, example value, hyperbolic
const _STATE_ELEMENT_QUANTITIES = (
    # --- spherical ------------------------------------------------------------------------------
    (:right_ascension, :SphericalAZIFPAState, :ra, "Right ascension", "rad", true,
     "right ascension of the position", "",
     "Measured in the equatorial plane of `cs` from its X axis, in [0, 2π).", "1.0", false),
    (:declination, :SphericalAZIFPAState, :dec, "Declination", "rad", false,
     "declination of the position", "",
     "Measured from the equatorial plane of `cs`, in [−π/2, π/2].", "0.2", false),
    (:velocity_right_ascension, :SphericalRADECState, :rav, "Velocity right ascension", "rad", true,
     "right ascension of the velocity", "",
     "The direction of the velocity vector, measured as `right_ascension` measures the position.",
     "1.5", false),
    (:velocity_declination, :SphericalRADECState, :decv, "Velocity declination", "rad", false,
     "declination of the velocity", "",
     "The direction of the velocity vector, measured as `declination` measures the position.",
     "0.1", false),
    (:velocity_azimuth, :SphericalAZIFPAState, :vazi, "Velocity azimuth", "rad", true,
     "azimuth of the velocity", "",
     "Measured east of north in the local horizontal plane, in [0, 2π). In a body-fixed `cs` it is " *
     "the azimuth relative to the rotating body.", "1.2", false),
    (:flight_path_angle, :SphericalAZIFPAState, :fpa, "Flight path angle", "rad", false,
     "flight path angle", "",
     "Measured from the radial direction to the velocity, in [0, π]: π/2 is horizontal flight and " *
     "less than π/2 is climbing. This is GMAT's convention; the angle above the local horizontal " *
     "is π/2 minus this one.", "1.5", false),
    # --- equinoctial ----------------------------------------------------------------------------
    (:equinoctial_h, :EquinoctialState, :h, "Equinoctial h", "", false,
     "equinoctial h", "e sin(ω + Ω)",
     "The same in `EquinoctialState` and `AlternateEquinoctialState`. Defined for closed orbits.",
     "0.01", false),
    (:equinoctial_k, :EquinoctialState, :k, "Equinoctial k", "", false,
     "equinoctial k", "e cos(ω + Ω)",
     "The same in `EquinoctialState` and `AlternateEquinoctialState`. Defined for closed orbits.",
     "0.01", false),
    (:equinoctial_p, :EquinoctialState, :p, "Equinoctial p", "", false,
     "equinoctial p", "tan(i/2) sin Ω",
     "Not the p of `ModifiedEquinoctialState`, which is the semi-latus rectum. Undefined at " *
     "i = π.", "0.1", false),
    (:equinoctial_q, :EquinoctialState, :q, "Equinoctial q", "", false,
     "equinoctial q", "tan(i/2) cos Ω", "Undefined at i = π.", "0.1", false),
    (:alternate_equinoctial_p, :AlternateEquinoctialState, :altp, "Alternate equinoctial p", "",
     false, "alternate equinoctial p", "sin(i/2) sin Ω",
     "Stays finite at i = π, where `equinoctial_p` does not.", "0.1", false),
    (:alternate_equinoctial_q, :AlternateEquinoctialState, :altq, "Alternate equinoctial q", "",
     false, "alternate equinoctial q", "sin(i/2) cos Ω",
     "Stays finite at i = π, where `equinoctial_q` does not.", "0.1", false),
    (:mean_longitude, :EquinoctialState, :mlong, "Mean longitude", "rad", true,
     "mean longitude", "Ω + ω + M",
     "In [0, 2π). The same in `EquinoctialState` and `AlternateEquinoctialState`. Stays defined on " *
     "circular and equatorial orbits, where its three terms separately do not.", "1.0", false),
    # --- modified equinoctial -------------------------------------------------------------------
    (:modified_equinoctial_f, :ModifiedEquinoctialState, :f, "Modified equinoctial f", "", false,
     "modified equinoctial f", "e cos(ω + Ω)", "Equal to `equinoctial_k`.", "0.01", false),
    (:modified_equinoctial_g, :ModifiedEquinoctialState, :g, "Modified equinoctial g", "", false,
     "modified equinoctial g", "e sin(ω + Ω)", "Equal to `equinoctial_h`.", "0.01", false),
    (:modified_equinoctial_h, :ModifiedEquinoctialState, :h, "Modified equinoctial h", "", false,
     "modified equinoctial h", "tan(i/2) cos Ω",
     "Not the h of `EquinoctialState`: it equals `equinoctial_q`.", "0.1", false),
    (:modified_equinoctial_k, :ModifiedEquinoctialState, :k, "Modified equinoctial k", "", false,
     "modified equinoctial k", "tan(i/2) sin Ω",
     "Not the k of `EquinoctialState`: it equals `equinoctial_p`.", "0.1", false),
    # --- Delaunay -------------------------------------------------------------------------------
    (:delaunay_l, :DelaunayState, :l, "Delaunay l", "rad", true,
     "Delaunay l", "the mean anomaly M", "Defined for closed orbits.", "1.0", false),
    (:delaunay_g, :DelaunayState, :g, "Delaunay g", "rad", true,
     "Delaunay g", "the argument of periapsis ω", "Defined for closed orbits.", "1.0", false),
    (:delaunay_h, :DelaunayState, :h, "Delaunay h", "rad", true,
     "Delaunay h", "the right ascension of the ascending node Ω", "Defined for closed orbits.",
     "1.0", false),
    (:delaunay_L, :DelaunayState, :L, "Delaunay L", "km²/s", false,
     "Delaunay L", "√(μa)", "Defined for closed orbits.", "54000.0", false),
    (:delaunay_G, :DelaunayState, :G, "Delaunay G", "km²/s", false,
     "Delaunay G", "L √(1 − e²), the angular momentum magnitude",
     "Defined for closed orbits; at most L.", "52000.0", false),
    (:delaunay_H, :DelaunayState, :H, "Delaunay H", "km²/s", false,
     "Delaunay H", "G cos i, the polar component of the angular momentum",
     "Defined for closed orbits; between −G and G.", "48000.0", false),
    # --- asymptotes -----------------------------------------------------------------------------
    (:outgoing_asymptote_right_ascension, :OutGoingAsymptoteState, :rla,
     "Outgoing asymptote right ascension", "rad", true,
     "right ascension of the outgoing asymptote", "",
     "In [0, 2π). On a closed orbit the apoapsis direction stands in for the asymptote.",
     "1.0", true),
    (:outgoing_asymptote_declination, :OutGoingAsymptoteState, :dla,
     "Outgoing asymptote declination", "rad", false,
     "declination of the outgoing asymptote", "",
     "In [−π/2, π/2]. On a closed orbit the apoapsis direction stands in for the asymptote.",
     "0.2", true),
    (:outgoing_b_vector_azimuth, :OutGoingAsymptoteState, :bpa,
     "Outgoing B-vector azimuth", "rad", true,
     "B-vector azimuth of the outgoing asymptote", "",
     "Counter-clockwise from south, in [0, 2π) (GMAT's OutgoingBVAZI).", "1.0", true),
    (:incoming_asymptote_right_ascension, :IncomingAsymptoteState, :rla,
     "Incoming asymptote right ascension", "rad", true,
     "right ascension of the incoming asymptote", "",
     "In [0, 2π). The asymptote points from the body back along the approach, opposite to the " *
     "incoming velocity at infinity. On a closed orbit the apoapsis direction stands in for it.",
     "1.0", true),
    (:incoming_asymptote_declination, :IncomingAsymptoteState, :dla,
     "Incoming asymptote declination", "rad", false,
     "declination of the incoming asymptote", "",
     "In [−π/2, π/2]. On a closed orbit the apoapsis direction stands in for the asymptote.",
     "0.2", true),
    (:incoming_b_vector_azimuth, :IncomingAsymptoteState, :bpa,
     "Incoming B-vector azimuth", "rad", true,
     "B-vector azimuth of the incoming asymptote", "",
     "Counter-clockwise from south, in [0, 2π) (GMAT's IncomingBVAZI). The B-plane quantities " *
     "`b_dot_t`, `b_dot_r` and `b_plane_angle` come from it.", "1.0", true),
)

# The Brouwer-Lyddane mean elements, long-period and short-period, six each.
const _BROUWER_THEORY_NOTE =
    "The conversion runs the mean-element theory and iterates to a tolerance, so a round trip " *
    "agrees to about 1e-4 km in the semi-major axis rather than to machine precision. The theory " *
    "holds for `0 ≤ e < 0.99` and a mean periapsis radius above 3000 km."
const _BROUWER_QUANTITIES = Tuple(
    (Symbol("brouwer_", kind, "_", name), S, field,
     "Brouwer $kind-period mean $what", unit, cyclic,
     "Brouwer-Lyddane $kind-period mean $what", "", _BROUWER_THEORY_NOTE, example, false)
    for (kind, S) in (("long", :BrouwerMeanLongState), ("short", :BrouwerMeanShortState))
    for (name, field, what, unit, cyclic, example) in (
        ("semi_major_axis",       :sma,  "semi-major axis",       "km",  false, "7000.0"),
        ("eccentricity",          :ecc,  "eccentricity",          "",    false, "0.01"),
        ("inclination",           :inc,  "inclination",           "rad", false, "0.5"),
        ("raan",                  :raan, "RAAN",                  "rad", true,  "1.0"),
        ("argument_of_periapsis", :aop,  "argument of periapsis", "rad", true,  "1.0"),
        ("mean_anomaly",          :ma,   "mean anomaly",          "rad", true,  "1.0")))

for (fn, S, field, _...) in (_STATE_ELEMENT_QUANTITIES..., _BROUWER_QUANTITIES...)
    _FIELD_QUANTITY[(S, field)] = fn
end

"The quantities a setter of `field` in `S` holds, as a list for its docstring."
function _held_list(S::Symbol, field::Symbol)
    others = [f for f in fieldnames(getfield(AstroStates, S)) if f !== field]
    names  = ["`$(_FIELD_QUANTITY[(S, f)])`" for f in others]
    return join(names[1:end-1], ", ") * " and " * names[end]
end

const _CLOSED_EXAMPLE = "sat = Spacecraft(state = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]))"
const _HYPERBOLIC_EXAMPLE = "sat = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 12.5, 1.0]))"

for (fn, S, field, lbl, u, cyclic, what, definition, note, example, hyperbolic) in
        (_STATE_ELEMENT_QUANTITIES..., _BROUWER_QUANTITIES...)
    fn!       = Symbol(fn, "!")
    defined   = isempty(definition) ? "" : ", $definition"
    unit_text = isempty(u) ? "" : " $u"
    returns   = isempty(u) ? "The dimensionless $what." : "The $what in $u."
    make_sat  = hyperbolic ? _HYPERBOLIC_EXAMPLE : _CLOSED_EXAMPLE
    reader_doc = """
        $fn(subject[, cs][, params])

    The $what of `subject`$defined, measured in `cs` or the subject's own frame.

    # Arguments
    See [`position_vector`](@ref).

    # Notes
    - $note
    - The `$field` field of `$S`. Settable with [`$fn!`](@ref).

    # Returns
    $returns

    # Examples
    ```julia
    $make_sat
    $fn(sat)
    $fn(sat, EarthMJ2000Ec)
    ```
    """
    setter_doc = """
        $fn!(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to`$unit_text, measured in `cs` or the spacecraft's own frame. The
    rest of `$S` is held in that frame: $(_held_list(S, field)).

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when the elements that
    result do not describe an orbit.

    # Examples
    ```julia
    $make_sat
    $fn!(sat; to = $example)
    $fn!(sat, EarthMJ2000Ec; to = $example)
    ```
    """
    @eval begin
        $fn(subject, cs::AbstractCoordinateSystem, params::NamedTuple) =
            _read_field(subject, $S, cs, params, $(QuoteNode(field)), $(uppercasefirst(what)))
        $fn(subject, cs::AbstractCoordinateSystem) = $fn(subject, cs, NamedTuple())
        $fn(subject) = $fn(subject, frame_of(subject))
        @doc $reader_doc $fn

        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::Real)
            isfinite(to) || throw(ArgumentError(
                $("$what must be finite; got ") * string(to) * ". The state was not changed."))
            return set_element!(sat, $S, cs, $(QuoteNode(field)), params; to = to)
        end
        @doc $setter_doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)

        EpicycleBase.label(::typeof($fn)) = $lbl
        subject_type(::typeof($fn)) = Union{Spacecraft, Coordinate}
    end
    isempty(u) || @eval EpicycleBase.unit(::typeof($fn)) = $u
    if cyclic
        @eval EpicycleBase.is_cyclic(::typeof($fn)) = true
        @eval EpicycleBase.cycle(::typeof($fn))     = 2π
    end
end

# --- setters for elements read through another representation ------------------------------------
#
# The semi-latus rectum and the true longitude are fields of ModifiedEquinoctialState, read in
# orbit_parameters.jl. C3 has no setter: the outgoing and incoming asymptote sets would hold
# different things, and neither is the obvious one.

for (fn, field, what, unit, example) in
        ((:semilatus_rectum, :p, "semi-latus rectum", " km",  "7000.0"),
         (:true_longitude,   :L, "true longitude",    " rad", "1.0"))
    fn! = Symbol(fn, "!")
    doc = """
        $fn!(sat[, cs][, params]; to) -> sat

    Set the $what of `sat` to `to`$unit, measured in `cs` or the spacecraft's own frame. The
    rest of `ModifiedEquinoctialState` is held in that frame: $(_held_list(:ModifiedEquinoctialState, field)).

    # Returns
    `sat`, mutated. Throws `ArgumentError`, and leaves `sat` unchanged, when the elements that
    result do not describe an orbit.

    # Examples
    ```julia
    $_CLOSED_EXAMPLE
    $fn!(sat; to = $example)
    ```
    """
    @eval begin
        function $fn!(sat::Spacecraft, cs::AbstractCoordinateSystem = frame_of(sat),
                      params::NamedTuple = NamedTuple(); to::Real)
            isfinite(to) || throw(ArgumentError(
                $("$what must be finite; got ") * string(to) * ". The state was not changed."))
            return set_element!(sat, ModifiedEquinoctialState, cs, $(QuoteNode(field)), params;
                                to = to)
        end
        @doc $doc $fn!
        EpicycleBase.set_quantity!(sat::Spacecraft, ::typeof($fn),
                                   cs::AbstractCoordinateSystem = frame_of(sat),
                                   params::NamedTuple = NamedTuple(); to) =
            $fn!(sat, cs, params; to = to)
    end
end

# --- the conditions each representation's elements must meet -------------------------------------
#
# Checked by `set_element!` before converting back, because a conversion does not refuse every
# element it cannot represent: a declination past the pole, for one, converts to a state that
# reads back as a different declination.

_refuse(msg) = throw(ArgumentError(msg * " The state was not changed."))

function _check_elements(s::SphericalAZIFPAState)
    s.r > 0 || _refuse("a position magnitude must be positive; got $(s.r) km.")
    s.v >= 0 || _refuse("a velocity magnitude must be zero or more; got $(s.v) km/s.")
    abs(s.dec) <= π/2 || _refuse("a declination must be between −π/2 and π/2 rad; got $(s.dec).")
    0 <= s.fpa <= π || _refuse("a flight path angle must be between 0 and π rad; got $(s.fpa).")
    return nothing
end

function _check_elements(e::Union{EquinoctialState, AlternateEquinoctialState})
    e.a > 0 || _refuse("the equinoctial semi-major axis must be positive; got $(e.a) km.")
    e.h^2 + e.k^2 < 1 || _refuse(
        "equinoctial h and k give eccentricity $(hypot(e.h, e.k)); the equinoctial elements " *
        "describe closed orbits only.")
    return nothing
end

function _check_elements(m::ModifiedEquinoctialState)
    m.p > 0 || _refuse("a semi-latus rectum must be positive; got $(m.p) km.")
    return nothing
end

function _check_elements(d::DelaunayState)
    d.L > 0 || _refuse("Delaunay L must be positive; got $(d.L) km²/s.")
    0 < d.G <= d.L || _refuse(
        "Delaunay G must be positive and at most L ($(d.L) km²/s); got $(d.G) km²/s.")
    abs(d.H) <= d.G || _refuse(
        "Delaunay H must be between −G and G ($(d.G) km²/s); got $(d.H) km²/s.")
    return nothing
end

function _check_elements(a::Union{OutGoingAsymptoteState, IncomingAsymptoteState})
    a.rp > 0 || _refuse("a periapsis radius must be positive; got $(a.rp) km.")
    abs(a.dla) <= π/2 || _refuse(
        "an asymptote declination must be between −π/2 and π/2 rad; got $(a.dla).")
    return nothing
end
