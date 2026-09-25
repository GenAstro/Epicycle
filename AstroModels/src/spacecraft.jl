# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

"""
    mutable struct Spacecraft{S<:OrbitState, TT<:Time, CS<:AbstractCoordinateSystem, T<:Real}

`Spacecraft` stores an orbital state, epoch, coordinate system, physical
properties, visualization model, and optional trajectory history.

# Fields
- state::S — orbital state as an OrbitState struct
- time::TT — epoch as a Time struct
- mass::T — total mass [kg]. **Private.** `sc.mass` throws; read it with `total_mass(sc)`.
  The field is replaced by a mass model in a later release and the accessor is what keeps
  callers unchanged across that.
- name::String — user label.
- history::SpacecraftHistory — trajectory history organized into segments
- coord_sys::CS — coordinate system (origin and axes) associated with the spacecraft.
- cad_model::CADModel — 3D model for visualization
- drag::Union{AbstractDragGeometry, Nothing} — drag geometry (e.g. `SphericalDrag`), or `nothing`
- srp::Union{AbstractSRPGeometry, Nothing} — SRP geometry (e.g. `SphericalSRP`), or `nothing`
- save_history::Bool — whether propagators append a segment to `history` on each call.
  Default `true`; set `false` to skip the ephemeris save (useful for benchmarking or
  when only the final state is wanted).

  When `false`, `propagate!` still updates `sc.state` and `sc.time` to the final
  integration point exactly as it does with `true`, but does not push a
  `HistorySegment` onto `sc.history`.
- notified::Set{Symbol} — which once-per-spacecraft notices have already been issued.
  Managed by `notify_once!`; a promoted copy shares it so a solver warns once rather
  than once per iterate.

# Notes
- The keyword constructor accepts any subset of fields that differ from the defaults.
- State can be provided two ways as shown in the example below. Assigning `sc.state` accepts
  any state representation and converts it to the representation and numeric type the
  spacecraft already holds, taking μ from the coordinate system origin when needed.
- history stores trajectory data in segments; use `history.segments` to access individual HistorySegments
- Numeric parameter T is chosen by promotion: T = promote_type(eltype(state), typeof(time.jd1), typeof(mass)).

# Examples
```julia
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse

sc = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T00:00:00", TAI(), ISOT()),
    mass  = 1000.0,
    name  = "Demo"
)

# alternative state definition
sc = Spacecraft(state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03],Cartesian()))
```
"""
mutable struct Spacecraft{S<:OrbitState, TT<:Time, CS<:AbstractCoordinateSystem, T<:Real} <: AbstractPoint
    state::S
    time::TT
    mass::T
    name::String
    history::SpacecraftHistory
    coord_sys::CS
    cad_model::CADModel
    drag::Union{AbstractDragGeometry, Nothing}
    srp::Union{AbstractSRPGeometry, Nothing}
    save_history::Bool
    notified::Set{Symbol}
end

"""
    total_mass(sc::Spacecraft) -> Real

Total current spacecraft mass [kg].

# Arguments
- `sc`: The spacecraft.

# Returns
The total mass [kg], in the spacecraft's numeric type.

# Notes
This is the only supported way to read mass. The field it reads is private and
`sc.mass` throws, because the field is replaced by a mass model in a later release
and this accessor is what keeps callers unchanged across that. Today it returns a
single pool; once tanks are modeled it becomes dry mass plus the propellant
remaining in each tank.

Mass changes only as a consequence of applying a maneuver. There is no setter,
because assigning a total would require deciding which tank it came from and the
caller cannot know that.

# Example
```jldoctest
sc = Spacecraft(mass = 1500.0)
total_mass(sc)

# output
1500.0
```
"""
total_mass(sc::Spacecraft) = getfield(sc, :mass)

# The mass field is private. It is replaced by a mass model in a later release, so a
# caller holding `sc.mass` holds a number whose meaning is about to change; throwing
# means nobody carries a stale reading forward silently.
# Every other name falls through to `getfield`, which the compiler resolves away when
# the symbol is a literal, so field access costs what it always did.
@inline function Base.getproperty(sc::Spacecraft, f::Symbol)
    if f === :mass
        throw(ErrorException(
            "`Spacecraft` has no public `mass` field. Use `total_mass(sc)` to read total " *
            "mass [kg]. Mass changes only by applying a maneuver."))
    end
    return getfield(sc, f)
end

@inline function Base.setproperty!(sc::Spacecraft, f::Symbol, v)
    if f === :mass
        throw(ErrorException(
            "`Spacecraft` has no public `mass` field, and mass has no setter. Assigning a " *
            "total would require deciding which tank it came from. Change mass by applying " *
            "a maneuver, or construct the spacecraft with the mass you want."))
    end
    f === :state && return setfield!(sc, :state, _as_orbit_state(sc, v))
    return setfield!(sc, f, convert(fieldtype(typeof(sc), f), v))
end

# `sc.state = CartesianState(...)` is what a user writes, so any state representation is
# accepted. The representation is part of the spacecraft's type, so the value is converted to
# the representation the spacecraft already holds, as `set_posvel!` does, and stored in the
# spacecraft's numeric type. A value already of the field's type goes straight in, which is the
# propagator's case.
function _as_orbit_state(sc::Spacecraft, v)
    v isa fieldtype(typeof(sc), :state) && return v
    v isa Union{AbstractState, OrbitState} || throw(ArgumentError(
        "Spacecraft state must be a state such as `CartesianState([...])` or " *
        "`KeplerianState(...)`, got a $(typeof(v)). To set position and velocity from a " *
        "vector, use `set_posvel!(sc, x)`."))
    os     = v isa OrbitState ? v : OrbitState(v)
    target = getfield(sc, :state).statetype
    T      = typeof(getfield(sc, :mass))
    os.statetype == target && return OrbitState(T.(copy(os.state)), target)

    from    = state_tag_to_type(os.statetype)(copy(os.state))
    tostate = state_tag_to_type(target)
    to      = applicable(tostate, from) ? tostate(from) :
              tostate(from, _origin_mu(sc, "Spacecraft state", os.statetype, target))
    return OrbitState(T.(to_vector(to)), target)
end

"""
    notify_once!(sc::Spacecraft, tag::Symbol, message::AbstractString) -> Bool

Issue `message` as a warning the first time `tag` is raised for `sc`, and stay silent
afterwards.

# Arguments
- `sc`: The spacecraft the notice concerns.
- `tag`: The notice's identity. `:lumped_mass_burn` and `:non_positive_mass` are the
  two in use.
- `message`: The warning text.

# Returns
`true` if the message was issued, `false` if `tag` had already been raised for this
spacecraft.

# Notes
The record lives on the spacecraft rather than in a table held by the caller. A table
keyed on object identity sees every promoted or reconstructed copy as a new
spacecraft, so under a solver it would warn once per iterate instead of once.
`Base.promote` shares this set with the promoted copy for the same reason;
`deepcopy` copies it, which carries the contents without aliasing.
"""
function notify_once!(sc::Spacecraft, tag::Symbol, message::AbstractString)
    tag in getfield(sc, :notified) && return false
    push!(getfield(sc, :notified), tag)
    @warn message
    return true
end

"""
    state_eltype(os::OrbitState) = eltype(os.state)

Returns the type of the elements in the state vector of an OrbitState.
"""
state_eltype(os::OrbitState) = eltype(os.state)

"""
    Spacecraft(state::Union{AbstractState,OrbitState}, time::TT; mass=1000.0, name="unnamed",
               history=nothing, coord_sys=CoordinateSystem(earth, ICRF())) where {TT<:Time}

Outer positional constructor for Spacecraft that promotes numeric types as needed.
"""
function Spacecraft(state::Union{AbstractState,OrbitState}, time::TT;
    mass::Real = 1000.0,
    name::AbstractString = "unnamed",
    history::Union{Nothing,SpacecraftHistory} = nothing,
    coord_sys::CS = CoordinateSystem(earth, ICRF()),
    cad_model::CADModel = CADModel(),
    drag::Union{AbstractDragGeometry, Nothing} = nothing,
    srp::Union{AbstractSRPGeometry, Nothing} = nothing,
    save_history::Bool = true,
    ) where {TT<:Time, CS<:AbstractCoordinateSystem}

    # Normalize to OrbitState
    os = state isa OrbitState ? state : OrbitState(state)

    # Determine numeric promotion type from state elements, time.jd1, and mass
    state_T = state_eltype(os)
    time_T  = typeof(time.jd1)
    mass_T  = typeof(mass)

    # Validate they are Real-typed (support Float, BigFloat, Dual, etc.)
    for (nm, Ty) in (("state elements", state_T), ("time.jd1", time_T), ("mass", mass_T))
        Ty <: Real || throw(ArgumentError("Spacecraft: $nm must be Real-typed, got $Ty"))
    end

    Tnum = promote_type(state_T, time_T, mass_T)

    # Convert state to T if needed
    os_T = state_eltype(os) === Tnum ? os : OrbitState(Tnum.(copy(os.state)), os.statetype)

    # Mass in T
    mass_T = Tnum(mass)

    # Rebuild time with promoted numeric type (preserve scale/format)
    t_T = AstroEpochs._time_jd(Tnum(time.jd1), Tnum(time.jd2),     # Julian-date parts, format kept
                               getfield(time, :scale), getfield(time, :format))
    TTIME = typeof(t_T)

    # History default: empty SpacecraftHistory
    hist_T = history === nothing ? SpacecraftHistory() : history

    return Spacecraft{typeof(os_T), TTIME, CS, Tnum}(os_T, t_T, mass_T, String(name), hist_T,
                                                     coord_sys, cad_model, drag, srp, save_history,
                                                     Set{Symbol}())
end

"""
    Spacecraft(; state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                      time = Time("2015-09-21T12:23:12", UTC(), ISOT()),
                      mass = 1000.0,
                      name = "unnamed",
                      history = nothing,
                      coord_sys = CoordinateSystem(earth, ICRF()))

Kwarg outer constructor for Spacecraft with defaults for all fields.
"""
function Spacecraft(; state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                      time = Time("2015-09-21T12:23:12", UTC(), ISOT()),
                      mass = 1000.0,
                      name = "unnamed",
                      history = nothing,
                      coord_sys = CoordinateSystem(earth, ICRF()),
                      cad_model = CADModel(),
                      drag = nothing,
                      srp = nothing,
                      save_history::Bool = true)
    Spacecraft(state, time; mass=mass, name=name, history=history, coord_sys=coord_sys,
               cad_model=cad_model, drag=drag, srp=srp, save_history=save_history)
end

"""
    Base.show(io::IO, sc::Spacecraft)

Pretty-print a Spacecraft in a human-readable, multi-line summary.
"""
function Base.show(io::IO, sc::Spacecraft)
     println(io, "Spacecraft: ", sc.name)
     _indent_and_print(io, sc.time, "  ")
     _indent_and_print(io, sc.state, "  ")
     _indent_and_print(io, sc.coord_sys, "  ")
     println(io, "  Total Mass = ", total_mass(sc), " kg")
     sc.drag === nothing ? println(io, "  Drag = none") : _indent_and_print(io, sc.drag, "  ")
     sc.srp  === nothing ? println(io, "  SRP  = none") : _indent_and_print(io, sc.srp,  "  ")
     println(io, "  Save History = ", sc.save_history)
     _indent_and_print(io, sc.cad_model, "  ")
 end

"""
    _indent_and_print(io::IO, obj, prefix::AbstractString)  

Indent and print composed objects using their own show methods
"""
 function _indent_and_print(io::IO, obj, prefix::AbstractString)
    # Use the MIME"text/plain" show to support types without 1-arg show
    s = repr(MIME"text/plain"(), obj)
    for line in split(chomp(s), '\n')
        println(io, prefix, line)
    end
end

"""
    Base.deepcopy_internal(sc::Spacecraft, dict::IdDict)

Deep copy a spacecraft to ensure no aliasing of inner mutable fields
"""
function Base.deepcopy_internal(sc::Spacecraft, dict::IdDict)
    # coord_sys, cad_model and the geometries are shared rather than copied. The geometries and
    # the CAD model are immutable. A coordinate system's origin can be another spacecraft, and
    # copying it would give the copy a private duplicate of that spacecraft.
    out = Spacecraft(
        state         = Base.deepcopy_internal(getfield(sc, :state), dict),
        time         = Base.deepcopy_internal(getfield(sc, :time), dict),
        mass         = getfield(sc, :mass),
        name         = getfield(sc, :name),
        history      = Base.deepcopy_internal(getfield(sc, :history), dict),
        coord_sys    = getfield(sc, :coord_sys),
        cad_model    = getfield(sc, :cad_model),
        drag         = getfield(sc, :drag),
        srp          = getfield(sc, :srp),
        save_history = getfield(sc, :save_history),
    )
    union!(getfield(out, :notified), getfield(sc, :notified))
    return out
end

"""
    get_state(sc::Spacecraft, target::AbstractOrbitStateType) -> AbstractOrbitState

Return the spacecraft's orbital state as the concrete type specified by `target`,
converting if needed. Does not mutate `sc`.

# Arguments
- `sc::Spacecraft`: the spacecraft.
- `target::AbstractOrbitStateType`: the representation wanted, such as `Keplerian()`.

# Notes
- If the conversion requires the gravitational parameter μ, it is taken from
  `sc.coord_sys.origin.mu`. An origin without one, such as another spacecraft,
  raises an `ArgumentError`.
- If the current state already has the `target` representation, a copy of it is
  returned.
- State is returned in the coordinate system of the spacecraft; no coordinate
  transformations are performed. Use `CartesianState(sc, cs)` for that.

# Returns
The state as the concrete type for `target`, such as `KeplerianState`.

# Examples
```julia
using AstroModels, AstroStates, AstroEpochs
sc = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
                time  = Time("2015-09-21T12:23:12", TAI(), ISOT()))
kep = get_state(sc, Keplerian())    # μ from the default Earth origin
```
"""
function get_state(sc::Spacecraft, target::AbstractOrbitStateType)::AbstractOrbitState

    # Convert the state tag to a concrete type caller
    tostate = state_tag_to_type(target)  

    # Create a concrete state from the current spacecraft's OrbitState
    concrete_state_caller = state_tag_to_type(sc.state.statetype)                   
    from_state = concrete_state_caller(copy(sc.state.state))

    # Fast path: already in required concrete type, nothing to do. 
    if from_state isa tostate
        return from_state
    end

    # A conversion that needs no μ has a one-argument method; otherwise μ comes from the origin.
    applicable(tostate, from_state) && return tostate(from_state)
    return tostate(from_state, _origin_mu(sc, "get_state", sc.state.statetype, target))
end

# The origin's gravitational parameter, for a conversion that needs one. Checked rather than
# caught, so an error raised inside a conversion reaches the caller as itself.
function _origin_mu(sc::Spacecraft, caller, from, to)
    origin = sc.coord_sys.origin
    hasfield(typeof(origin), :mu) && return getfield(origin, :mu)
    throw(ArgumentError(
        "$caller: μ is required to convert from $from to $to, but the coordinate system " *
        "origin has no gravitational parameter. Use a celestial body such as `earth` as " *
        "the origin, or a state representation that does not need μ."))
end

"""
    to_posvel(sc::Spacecraft) -> Vector{<:Real}

Return the Cartesian position-velocity vector [x, y, z, vx, vy, vz] for `sc` in its
current coordinate system, in km and km/s, converting the stored orbital state if needed.

# Arguments
- `sc::Spacecraft`: the spacecraft.

# Notes
A state stored in another representation is converted as [`get_state`](@ref)
converts it, taking μ from the coordinate system origin when the conversion needs it.
The returned vector is a copy; changing it does not change `sc`.

# Returns
A 6-element vector in the spacecraft's numeric type.

# Examples
```julia
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroUniverse
sc = Spacecraft(
           state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
           time=Time("2015-09-21T12:23:12", TAI(), ISOT()));
to_posvel(sc)

# output
6-element Vector{Float64}:
 7000.0
  300.0
    0.0
    0.0
    7.5
    0.03
```
"""
function to_posvel(sc::Spacecraft)
    
    # If orbit state with type Cartesian(), this is fast just return
    if sc.state isa OrbitState && sc.state.statetype == Cartesian()
        return copy(sc.state.state)
    end

    # For any other state type, convert to Cartesian using get_state
    # (which handles μ extraction from coord_sys.origin when needed)
    cart_state = get_state(sc, Cartesian())
    return cart_state.posvel
end

"""
    set_posvel!(sc::Spacecraft, x::AbstractVector{<:Real})

Set the Cartesian position-velocity vector [x, y, z, vx, vy, vz] for `sc` in-place,
in km and km/s, in the spacecraft's current coordinate system.

# Arguments
- `sc::Spacecraft`: the spacecraft.
- `x::AbstractVector{<:Real}`: the new position and velocity, six elements.

# Notes
The spacecraft keeps its state representation: a spacecraft holding Keplerian
elements is given the elements of `x`, converted with μ from the coordinate system
origin. A vector that is not six elements long, or an origin without μ when the
conversion needs one, raises an `ArgumentError`.

# Returns
`nothing`. The function replaces `sc.state` while preserving the spacecraft's
numeric type.

# Examples
```julia
using AstroModels, AstroStates, AstroEpochs
sc = Spacecraft(
           state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
           time=Time("2015-09-21T12:23:12", TAI(), ISOT())
       );
set_posvel!(sc, [7050.0, 0.0, 0.0, 0.0, 7.6, 0.0]);
to_posvel(sc)

# output
6-element Vector{Float64}:
 7050.0
    0.0
    0.0
    0.0
    7.6
    0.0
```
"""
function set_posvel!(sc::Spacecraft, x::AbstractVector{<:Real})

    # Validate input length 
    if length(x) != 6
        throw(ArgumentError("set_posvel!: expected a length-6 vector, got length=$(length(x))"))
    end

    # If orbit state with type Cartesian(), this is fast just return
    if (sc.state isa OrbitState) && (sc.state.statetype == Cartesian())
        # Preserve the spacecraft's numeric type T
        T = typeof(total_mass(sc))
        sc.state = OrbitState(T.(x), Cartesian())
        return
    end

    # For other state types, convert the Cartesian vector to the current state type
    cart_state = CartesianState(x)
    TargetType = state_tag_to_type(sc.state.statetype)
    target_state = applicable(TargetType, cart_state) ? TargetType(cart_state) :
        TargetType(cart_state, _origin_mu(sc, "set_posvel!", Cartesian(), sc.state.statetype))
    T = typeof(total_mass(sc))
    sc.state = OrbitState(T.(to_vector(target_state)), sc.state.statetype)
    return
end

"""
    Base.promote(sc::Spacecraft{S,TT,CS,T}, ::Type{Tnew}) where {S,TT,CS,T,Tnew<:Real}

Promotes a Spacecraft to a new numeric type `Tnew` for automatic differentiation support.
The state, time, and mass are promoted to `Tnew`, while history remains as Float64 for efficiency.

This enables AD workflows where computation types (e.g., ForwardDiff.Dual) are promoted
while preserving Float64 ephemeris storage.

# Arguments
- `sc::Spacecraft`: The spacecraft to promote
- `::Type{Tnew}`: Target numeric type (e.g., ForwardDiff.Dual{Nothing,Float64,3})

# Returns
- `Spacecraft{S_new, TT_new, CS, Tnew}`: Promoted spacecraft

# Example
```julia
using ForwardDiff
sc = Spacecraft(state=CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]))
sc_dual = promote(sc, ForwardDiff.Dual{Nothing,Float64,3})
```
"""
function Base.promote(sc::Spacecraft{S,TT,CS,T}, ::Type{Tnew}) where {S,TT,CS,T,Tnew<:Real}
    # Promote the state to new type
    state_promoted = OrbitState(Tnew.(copy(sc.state.state)), sc.state.statetype)
    
    # Promote the time to new type (preserve scale/format)
    time_promoted = AstroEpochs._time_jd(Tnew(sc.time.jd1), Tnew(sc.time.jd2),  # Julian-date parts, format kept
                        getfield(sc.time, :scale), getfield(sc.time, :format))
    
    # Promote the mass to new type
    mass_promoted = Tnew(total_mass(sc))
    
    # Keep history as Float64 (no promotion needed - already Float64)
    # This is the key benefit: ephemeris storage remains efficient
    history_preserved = deepcopy(sc.history)
    
    # Create new spacecraft with promoted types
    return Spacecraft{typeof(state_promoted), typeof(time_promoted), CS, Tnew}(
        state_promoted,
        time_promoted,
        mass_promoted,
        sc.name,
        history_preserved,
        sc.coord_sys,
        sc.cad_model,
        getfield(sc, :drag),
        getfield(sc, :srp),
        getfield(sc, :save_history),
        getfield(sc, :notified),
    )
end
