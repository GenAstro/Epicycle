"""
    mutable struct Spacecraft{S<:OrbitState, TT<:Time, CS<:AbstractCoordinateSystem, T<:Real}

Spacecraft struct with state, time, mass, name, history, and coordinate system and other data

Fields
- state::S — orbital state as an OrbitState struct
- time::TT — epoch as a Time struct
- mass::T — total mass 
- name::String — user label.
- history::Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}} — time-tagged history (e.g., [(time, posvel)]) grouped in segments.
- coord_sys::CS — coordinate system (origin and axes) associated with the spacecraft.
- cad_model::CADModel — 3D model for visualization

# Notes:
- Use the keyword constructor to create spacecraft and only define the fields you want to change from the defaults.
- State can be provided two ways as shown in the example below
- history is an internal structure to store ephemeris during integration.  It will be hidden in future releases.
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

# alternative epoch definition
sc = Spacecraft(state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03],Cartesian()))
```
"""
mutable struct Spacecraft{S<:OrbitState, TT<:Time, CS<:AbstractCoordinateSystem, T<:Real} <: AbstractPoint
    state::S
    time::TT
    mass::T
    name::String
    history::Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}}
    coord_sys::CS
    cad_model::CADModel
end

"""
    state_eltype(os::OrbitState) = eltype(os.state)

Returns the type of the elements in the state vector of an OrbitState.
"""
state_eltype(os::OrbitState) = eltype(os.state)

"""
    Spacecraft(state::Union{AbstractState,OrbitState}, time::TT; mass=1000.0, name="unnamed",
               history=nothing, coord_sys=CoordinateSystem(earth, ICRFAxes())) where {TT<:Time}

Outer positional constructor for Spacecraft that promotes numeric types as needed.
"""
function Spacecraft(state::Union{AbstractState,OrbitState}, time::TT;
    mass::Real = 1000.0,
    name::AbstractString = "unnamed",
    history::Union{Nothing,AbstractVector} = nothing,
    coord_sys::CS = CoordinateSystem(earth, ICRFAxes()),
    cad_model::CADModel = CADModel()
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
    t_T = Time(Tnum(time.jd1), Tnum(time.jd2),
               getfield(time, :scale), getfield(time, :format))
    TTIME = typeof(t_T)

    # History default with Float64 types (always, regardless of T)
    hist_T = history === nothing ?
        Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}}() :
        convert(Vector{Vector{Tuple{Time{Float64}, Vector{Float64}}}}, history)

    return Spacecraft{typeof(os_T), TTIME, CS, Tnum}(os_T, t_T, mass_T, String(name), hist_T, coord_sys, cad_model)
end

"""
    Spacecraft(; state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                      time = Time("2015-09-21T12:23:12", UTC(), ISOT()),
                      mass = 1000.0,
                      name = "unnamed",
                      history = nothing,
                      coord_sys = CoordinateSystem(earth, ICRFAxes()),
                      cad_model = CADModel())

Kwarg outer constructor for Spacecraft with defaults for all fields.
"""
function Spacecraft(; state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
                      time = Time("2015-09-21T12:23:12", UTC(), ISOT()),
                      mass = 1000.0,
                      name = "unnamed",
                      history = nothing,
                      coord_sys = CoordinateSystem(earth, ICRFAxes()),
                      cad_model = CADModel())
    Spacecraft(state, time; mass=mass, name=name, history=history, coord_sys=coord_sys, cad_model=cad_model)
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
     println(io, "  Total Mass = ", sc.mass, " kg")
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
    return Spacecraft(
        # TODO.  implement deep copy on composed objects and call here
        state     = Base.deepcopy_internal(getfield(sc, :state), dict),
        time      = Base.deepcopy_internal(getfield(sc, :time), dict),
        mass      = getfield(sc, :mass),
        name      = getfield(sc, :name),
        history   = Base.deepcopy_internal(getfield(sc, :history), dict),
        coord_sys = getfield(sc, :coord_sys),
        cad_model = getfield(sc, :cad_model),
    )
end

"""
    get_state(sc::Spacecraft, target::AbstractOrbitStateType) -> AbstractOrbitState

Return the spacecraft's orbital state as the concrete type specified by `target`,
converting if needed. Does not mutate `sc`.

Notes
- If conversion requires the gravitational parameter μ, it is taken from
  `sc.coord_sys.origin.mu` when available; otherwise an error is thrown.
- If the current state already matches `target`, the existing state is returned
  without conversion to the new type.
- State is returned in the coordinate system of the spacecraft; no coordinate
  transformations are performed.

# Examples

```julia
using AstroModels
sc = Spacecraft(
           state=CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
           time=Time("2015-09-21T12:23:12", TAI(), ISOT())
           );
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

    # Attempt conversion without μ first; if that fails, try with μ from the coord system origin
    try
        return tostate(from_state)
    catch
        origin = sc.coord_sys.origin
        T = typeof(origin)
        if hasfield(T, :mu)
            return tostate(from_state, getfield(origin, :mu))
        else
            typein = sc.state.statetype
            error("get_state: μ is required to convert from $typein to $target but the coordinate
            system does not have a celestial body (with a μ) as its origin.  Change to state types
            that do not require μ or change the coordinate system.")
        end
    end
end

"""
    to_posvel(sc::Spacecraft) -> Vector{<:Real}

Return the Cartesian position-velocity vector [x, y, z, vx, vy, vz] for `sc` in its
current coordinate system, converting the stored orbital state if needed.

Currently supported
- OrbitState with `statetype == Cartesian()` — returns the internal 6-vector (fast path).
- CartesianState — returns `to_vector(state)`.

Limitations
- Other state types currently throw an ArgumentError. Use `get_state(sc, Cartesian())`
  first or extend conversions.

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

    # TODO. Generalize this function, possibly replacing it.
    throw(ArgumentError("Conversion not implemented yet to convert state. TODO."))
end

"""
    set_posvel!(sc::Spacecraft, x::AbstractVector{<:Real})

Set the Cartesian position-velocity vector [x, y, z, vx, vy, vz] for `sc` in-place.

Supported
- `OrbitState` with `statetype == Cartesian()`.

Limitations
- Other state representations currently throw an ArgumentError. 
  Use `get_state(sc, Cartesian())` first.

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
        T = typeof(sc.mass)
        sc.state = OrbitState(T.(x), Cartesian())
        return
    end

    # TODO. Generalize this function, possibly replacing it.
    throw(ArgumentError("Conversion not implemented yet to convert state. TODO."))

end

"""
    push_history_segment!(sc::Spacecraft, segment::Vector{<:Tuple})

Append a new history segment to the spacecraft.
Always stores data as Float64 regardless of input numeric types.

# Arguments
- sc::Spacecraft: The spacecraft to update
- segment::Vector{Tuple{Time, Vector}}: Complete segment with time-state pairs

# Returns
- sc::Spacecraft: The same spacecraft instance (for chaining)

# Notes
- Converts all times and position/velocity vectors to Float64 for efficient storage
- History is for ephemeris logging, not differentiation
- Works for single-point segments (maneuvers) or multi-point segments (propagation)
"""
function push_history_segment!(sc::Spacecraft, segment::Vector{<:Tuple})
    # Helper function to safely convert any Real to Float64 (handles Dual numbers)
    function to_float64(x::Real)
        # For regular numbers
        if x isa AbstractFloat || x isa Integer
            return Float64(x)
        # For ForwardDiff.Dual and other AD types with .value field
        elseif hasfield(typeof(x), :value)
            return Float64(x.value)
        else
            # Fallback: try to convert directly
            return Float64(x)
        end
    end
    
    # Convert entire segment to Float64 types for efficient storage
    segment_f64 = [(Time(to_float64(t.jd1), to_float64(t.jd2), t.scale, t.format), 
                    map(to_float64, pv)) 
                   for (t, pv) in segment]
    push!(sc.history, segment_f64)
    return sc
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
    time_promoted = Time(Tnew(sc.time.jd1), Tnew(sc.time.jd2),
                        getfield(sc.time, :scale), getfield(sc.time, :format))
    
    # Promote the mass to new type
    mass_promoted = Tnew(sc.mass)
    
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
        sc.cad_model
    )
end

