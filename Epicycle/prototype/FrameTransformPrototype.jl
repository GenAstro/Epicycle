# Frame Transform Prototype
# Prototyping frame transformation design before integration into AstroFrames

using StaticArrays
using LinearAlgebra
using SatelliteToolboxTransformations
using SatelliteToolboxTransformations: OrbitStateVector, EopIau1980, EopIau2000A, 
                                        fetch_iers_eop, sv_ecef_to_eci, sv_eci_to_ecef
using AstroModels
using AstroModels: Spacecraft
using AstroStates
using AstroStates: CartesianState
using AstroUniverse
using AstroEpochs
using AstroEpochs: Time, TAI, UTC, ISOT
using EpicycleBase: AbstractPoint
using Base.Threads: ReentrantLock

# =============================================================================
# Core Type Hierarchy
# =============================================================================

abstract type AbstractReferenceSystem end
struct FK5 <: AbstractReferenceSystem end
struct IAU2006 <: AbstractReferenceSystem end
struct Inherited <: AbstractReferenceSystem end

abstract type AbstractOrientationType end
struct Inertial <: AbstractOrientationType end
struct Fixed <: AbstractOrientationType end
struct Derived <: AbstractOrientationType end

abstract type AbstractAxes end

# =============================================================================
# Axes Types
# =============================================================================

# Inertial frames
struct ICRF <: AbstractAxes end   # IAU2006 (bias-corrected from GCRF)
struct CIRS <: AbstractAxes end   # IAU2006
struct GCRF <: AbstractAxes end   # FK5 inertial
struct J2000 <: AbstractAxes end  # FK5 mean equator/equinox of J2000
struct MODEq <: AbstractAxes end   # FK5 mean-of-date (equatorial)
struct TODEq <: AbstractAxes end   # FK5 true-of-date (equatorial)

# Earth-fixed frames
struct ITRF <: AbstractAxes end
struct TIRS <: AbstractAxes end
struct PEF <: AbstractAxes end

# Orbit-relative frames
struct LVLH <: AbstractAxes 
    reference::Spacecraft  # Spacecraft whose orbit defines the LVLH frame
end

# =============================================================================
# CoordinateFrame and Coordinate Structs
# =============================================================================

"""
    CoordinateFrame

Defines a reference frame with an origin point and orientation axes.
"""
struct CoordinateFrame
    origin::AbstractPoint
    axes::AbstractAxes
end

"""
    Coordinate{T<:Real}

Represents a state in space-time with position, velocity, acceleration,
time, and coordinate frame. Parametric type T supports Float64, ForwardDiff.Dual, etc.
for automatic differentiation compatibility.
"""
struct Coordinate{T<:Real}
    pos::SVector{3,T}
    vel::SVector{3,T}
    acc::SVector{3,T}
    time::Time
    frame::CoordinateFrame
end

# Convenience constructors
function Coordinate(pos::SVector{3,T}, vel::SVector{3,T}, time::Time, frame::CoordinateFrame) where T<:Real
    Coordinate{T}(
        pos,
        vel,
        SVector{3,T}(zero(T), zero(T), zero(T)),
        time,
        frame
    )
end

# Allow construction from any vector-like input
function Coordinate(pos, vel, time::Time, frame::CoordinateFrame)
    T = promote_type(eltype(pos), eltype(vel))
    Coordinate(
        SVector{3,T}(pos),
        SVector{3,T}(vel),
        SVector{3,T}(zero(T), zero(T), zero(T)),
        time,
        frame
    )
end

# Display method for Coordinate (requested compact format)
function Base.show(io::IO, coord::Coordinate{T}) where T
    # Format time as: ISOT (SCALE), Julian date
    isot = string(coord.time.isot)
    scale_tag = uppercase(String(coord.time.scale))
    jd_str = string(coord.time.jd)
    time_str = "$isot ($scale_tag), $jd_str"
    
    # Format origin compactly (just the name)
    origin_str = if hasfield(typeof(coord.frame.origin), :name)
        coord.frame.origin.name
    else
        string(coord.frame.origin)
    end
    
    # Format frame and axes information
    frame_name = "$(nameof(typeof(coord.frame.axes)))Frame"
    axes_str = if coord.frame.axes isa LVLH
        "LVLH(ref: $(coord.frame.axes.reference.name))"
    else
        string(nameof(typeof(coord.frame.axes)))
    end
    
    println(io, "Coordinate{$T}:")
    println(io, "  Time: $time_str")
    println(io, "  State:")
    println(io, "      r : $(coord.pos)  km")
    println(io, "      v : $(coord.vel)  km/s")
    if !iszero(coord.acc)
        println(io, "      a : $(coord.acc)  km/s²")
    end
    println(io, "  Frame : $frame_name")
    println(io, "      Origin = $origin_str")
    print(io, "      Axes   = $axes_str")
end

# =============================================================================
# Edge Transform Types (how to execute transformations)
# =============================================================================

abstract type EdgeTransform end

"""
    STBTransform
    
Transformation using SatelliteToolboxTransformations.jl function.
Contains all metadata needed to call STB functions.
"""
struct STBTransform <: EdgeTransform
    stb_function::Function      # sv_ecef_to_ecef, sv_eci_to_ecef, etc.
    val_from::Val               # ITRF(), PEF(), etc.
    val_to::Val                 # PEF(), ITRF(), etc.
    eop_type::Symbol            # :IAU1980, :IAU2000A, or :none
end

"""
    CustomTransform
    
Simple custom transformation (rotation only).
"""
struct CustomTransform <: EdgeTransform
    transform_fn::Function      # Your implementation
end

"""
    OrbitRelativeTransform
    
Transformation to/from orbit-relative frames (LVLH, VNB, TNW).
Requires the orbit state to define the frame.
"""
struct OrbitRelativeTransform <: EdgeTransform
    transform_fn::Function
end

# =============================================================================
# Placeholder Transform Functions (to be implemented)
# =============================================================================

function cirs_to_icrf(coord::Coordinate, target_frame::CoordinateFrame)
    # CIRS to ICRF: Frame bias correction
    # The bias between CIRS and ICRF is very small (~mas level)
    # For now, treating as identity. Full implementation would apply
    # IAU 2006 frame bias matrix.
    # TODO: Implement proper bias matrix from IERS conventions
    return Coordinate(coord.pos, coord.vel, coord.time, target_frame)
end

function icrf_to_cirs(coord::Coordinate, target_frame::CoordinateFrame)
    # ICRF to CIRS: Inverse frame bias correction
    # For now, treating as identity (bias is ~mas level)
    # TODO: Implement inverse bias matrix
    return Coordinate(coord.pos, coord.vel, coord.time, target_frame)
end

"""
    compute_lvlh_rotation_matrix(r, v)

Compute LVLH (Local Vertical Local Horizontal) rotation matrix from orbit state.

LVLH frame definition:
- z-axis: radial direction (r̂ = r/|r|) pointing away from central body
- y-axis: -ĥ (negative orbit normal) = -(r × v)/|r × v|  
- x-axis: ŷ × ẑ (completes right-hand system, roughly velocity direction)

Returns 3×3 rotation matrix R where columns are [x̂ ŷ ẑ] expressed in inertial frame.
To rotate from inertial to LVLH: v_lvlh = R' * v_inertial
To rotate from LVLH to inertial: v_inertial = R * v_lvlh
"""
function compute_lvlh_rotation_matrix(r::SVector{3,T}, v::SVector{3,T}) where T
    # Radial direction (z-axis of LVLH)
    r_norm = norm(r)
    z_hat = r / r_norm
    
    # Orbit normal (angular momentum direction)
    h = cross(r, v)
    h_norm = norm(h)
    h_hat = h / h_norm
    
    # y-axis: negative orbit normal (points "south" for prograde orbits)
    y_hat = -h_hat
    
    # x-axis: completes right-hand system (roughly along-track)
    x_hat = cross(y_hat, z_hat)
    
    # Rotation matrix with LVLH axes as columns
    return hcat(x_hat, y_hat, z_hat)
end

"""
    compute_lvlh_angular_velocity(r, v)

Compute angular velocity of LVLH frame.
ω = (r × v) / |r|²
"""
function compute_lvlh_angular_velocity(r::SVector{3,T}, v::SVector{3,T}) where T
    r_norm_sq = dot(r, r)
    return cross(r, v) / r_norm_sq
end

function icrf_to_lvlh(coord::Coordinate, orbit_state::Coordinate, target_frame::CoordinateFrame)
    # Build LVLH rotation matrix from orbit state
    R = compute_lvlh_rotation_matrix(orbit_state.pos, orbit_state.vel)
    
    # Rotate position from ICRF to LVLH
    pos_lvlh = R' * coord.pos
    
    # Rotate velocity and correct for frame rotation
    # v_lvlh = R' * v_icrf - ω_lvlh × r_lvlh
    vel_icrf_rotated = R' * coord.vel
    omega_lvlh = compute_lvlh_angular_velocity(orbit_state.pos, orbit_state.vel)
    omega_lvlh_rotated = R' * omega_lvlh
    vel_lvlh = vel_icrf_rotated - cross(omega_lvlh_rotated, pos_lvlh)
    
    # Rotate acceleration (if non-zero)
    # For now, simple rotation without Coriolis/centrifugal terms
    acc_lvlh = R' * coord.acc
    
    return Coordinate(pos_lvlh, vel_lvlh, acc_lvlh, coord.time, target_frame)
end

function lvlh_to_icrf(coord::Coordinate, orbit_state::Coordinate, target_frame::CoordinateFrame)
    # Build LVLH rotation matrix from orbit state
    R = compute_lvlh_rotation_matrix(orbit_state.pos, orbit_state.vel)
    
    # Rotate position from LVLH to ICRF
    pos_icrf = R * coord.pos
    
    # Rotate velocity and correct for frame rotation
    # v_icrf = R * (v_lvlh + ω_lvlh × r_lvlh)
    omega_lvlh = compute_lvlh_angular_velocity(orbit_state.pos, orbit_state.vel)
    omega_lvlh_rotated = R' * omega_lvlh
    vel_with_rotation = coord.vel + cross(omega_lvlh_rotated, coord.pos)
    vel_icrf = R * vel_with_rotation
    
    # Rotate acceleration
    acc_icrf = R * coord.acc
    
    return Coordinate(pos_icrf, vel_icrf, acc_icrf, coord.time, target_frame)
end

# =============================================================================
# Graph Structure, Metadata, and Caches
# =============================================================================

"""
    EdgeTransformMetadata

Metadata used to evaluate whether an edge/path is valid in the current transform context.
Also used in cache keys and stored as a snapshot in cache entries.
"""
struct EdgeTransformMetadata
    has_eop_iau1980::Bool
    has_eop_iau2000a::Bool
    has_spacecraft_reference::Bool
end

EdgeTransformMetadata() = EdgeTransformMetadata(true, true, true)

"""
    PathCacheKey

Cache key for both L1 and L2 caches.
"""
struct PathCacheKey
    from_axes::DataType
    to_axes::DataType
    edge_transform_metadata::EdgeTransformMetadata
    graph_version::UInt64
end

"""
    PrecomputedPathEntry

L2 cache entry: stores precomputed path + metadata snapshot.
"""
struct PrecomputedPathEntry
    path::Vector{DataType}
    edge_transform_metadata::EdgeTransformMetadata
    graph_version::UInt64
end

"""
    SimPlanEntry

L1 cache entry: stores compiled execution plan + metadata snapshot.
"""
struct SimPlanEntry
    path::Vector{DataType}
    edges::Vector{EdgeTransform}
    edge_transform_metadata::EdgeTransformMetadata
    graph_version::UInt64
end

"""
    EDGE_TRANSFORMS

Metadata for executing each edge transformation.
Maps (FromAxes, ToAxes) → EdgeTransform with execution details.
"""
const EDGE_TRANSFORMS = Dict{Tuple{DataType, DataType}, EdgeTransform}(
    # ITRF ↔ TIRS (modern IAU2006)
    (ITRF, TIRS) => STBTransform(
        SatelliteToolboxTransformations.sv_ecef_to_ecef,
        Val(:ITRF), Val(:TIRS), :IAU2000A
    ),
    (TIRS, ITRF) => STBTransform(
        SatelliteToolboxTransformations.sv_ecef_to_ecef,
        Val(:TIRS), Val(:ITRF), :IAU2000A
    ),
    
    # TIRS ↔ CIRS (modern IAU2006)
    (TIRS, CIRS) => STBTransform(
        SatelliteToolboxTransformations.sv_ecef_to_eci,
        Val(:TIRS), Val(:CIRS), :none
    ),
    (CIRS, TIRS) => STBTransform(
        SatelliteToolboxTransformations.sv_eci_to_ecef,
        Val(:CIRS), Val(:TIRS), :none
    ),
    
    # CIRS ↔ ICRF
    (CIRS, ICRF) => CustomTransform(cirs_to_icrf),
    (ICRF, CIRS) => CustomTransform(icrf_to_cirs),
    
    # ICRF ↔ LVLH (orbit-relative)
    (ICRF, LVLH) => OrbitRelativeTransform(icrf_to_lvlh),
    (LVLH, ICRF) => OrbitRelativeTransform(lvlh_to_icrf),
    
    # ITRF ↔ PEF (legacy FK5)
    (ITRF, PEF) => STBTransform(
        SatelliteToolboxTransformations.sv_ecef_to_ecef,
        Val(:ITRF), Val(:PEF), :IAU1980
    ),
    (PEF, ITRF) => STBTransform(
        SatelliteToolboxTransformations.sv_ecef_to_ecef,
        Val(:PEF), Val(:ITRF), :IAU1980
    ),
    
    # ITRF ↔ GCRF (FK5 direct transform - FAST PATH)
    (ITRF, GCRF) => STBTransform(
        sv_ecef_to_eci,
        Val(:ITRF), Val(:GCRF), :IAU1980
    ),
    (GCRF, ITRF) => STBTransform(
        sv_eci_to_ecef,
        Val(:GCRF), Val(:ITRF), :IAU1980
    ),

    # ITRF ↔ J2000 (FK5 mean equator/equinox - FAST PATH)
    (ITRF, J2000) => STBTransform(
        sv_ecef_to_eci,
        Val(:ITRF), Val(:J2000), :IAU1980
    ),
    (J2000, ITRF) => STBTransform(
        sv_eci_to_ecef,
        Val(:J2000), Val(:ITRF), :IAU1980
    ),

    # ITRF ↔ MODEq (FK5 mean-of-date, equatorial)
    # Intentionally not added to FAST_PATH_IMPLEMENTATIONS so this exercises graph walking
    (ITRF, MODEq) => STBTransform(
        sv_ecef_to_eci,
        Val(:ITRF), Val(:MOD), :IAU1980
    ),
    (MODEq, ITRF) => STBTransform(
        sv_eci_to_ecef,
        Val(:MOD), Val(:ITRF), :IAU1980
    ),

    # ITRF ↔ TODEq (FK5 true-of-date, equatorial)
    (ITRF, TODEq) => STBTransform(
        sv_ecef_to_eci,
        Val(:ITRF), Val(:TOD), :IAU1980
    ),
    (TODEq, ITRF) => STBTransform(
        sv_eci_to_ecef,
        Val(:TOD), Val(:ITRF), :IAU1980
    ),
    
    # GCRF ↔ LVLH (orbit-relative, FK5 inertial)
    # Note: GCRF is inertial, so same rotation logic as ICRF
    (GCRF, LVLH) => OrbitRelativeTransform(icrf_to_lvlh),
    (LVLH, GCRF) => OrbitRelativeTransform(lvlh_to_icrf),

    # J2000 ↔ LVLH (orbit-relative, FK5 inertial)
    (J2000, LVLH) => OrbitRelativeTransform(icrf_to_lvlh),
    (LVLH, J2000) => OrbitRelativeTransform(lvlh_to_icrf),
)

"""
    FAST_PATH_IMPLEMENTATIONS

Direct fast-path edge mappings checked before cache lookup / graph walk.
"""
const FAST_PATH_IMPLEMENTATIONS = Dict{Tuple{DataType, DataType}, EdgeTransform}(
    (ITRF, GCRF) => EDGE_TRANSFORMS[(ITRF, GCRF)],
    (GCRF, ITRF) => EDGE_TRANSFORMS[(GCRF, ITRF)],
    (ITRF, J2000) => EDGE_TRANSFORMS[(ITRF, J2000)],
    (J2000, ITRF) => EDGE_TRANSFORMS[(J2000, ITRF)],
)

"""
    PRECOMPUTED_PATH_CACHE

L2 cache for precomputed paths.
"""
const PRECOMPUTED_PATH_CACHE = Dict{PathCacheKey, PrecomputedPathEntry}()

"""
    SIM_PLAN_CACHE

L1 cache for simulation/session compiled execution plans.
"""
const SIM_PLAN_CACHE = Dict{PathCacheKey, SimPlanEntry}()

"""
    GRAPH_VERSION

Version of the registered transform graph. Increment on each register call.
"""
const GRAPH_VERSION = Ref{UInt64}(0)

const GRAPH_LOCK = ReentrantLock()
const CACHE_LOCK = ReentrantLock()

function clear_transform_caches!()
    lock(CACHE_LOCK)
    try
        empty!(PRECOMPUTED_PATH_CACHE)
        empty!(SIM_PLAN_CACHE)
    finally
        unlock(CACHE_LOCK)
    end
    return nothing
end

function graph_version()
    lock(GRAPH_LOCK)
    try
        return GRAPH_VERSION[]
    finally
        unlock(GRAPH_LOCK)
    end
end

"""
    register_axes_transformation!(from_axes, to_axes, edge; fast_path=false)

Register a new edge transform and increment graph version.
Edge metadata is immutable after registration (no update API provided).
"""
function register_axes_transformation!(from_axes::Type, to_axes::Type,
                                       edge::EdgeTransform; fast_path::Bool=false)
    lock(GRAPH_LOCK)
    try
        EDGE_TRANSFORMS[(from_axes, to_axes)] = edge
        if fast_path
            FAST_PATH_IMPLEMENTATIONS[(from_axes, to_axes)] = edge
        end
        GRAPH_VERSION[] += 1
    finally
        unlock(GRAPH_LOCK)
    end
    return nothing
end

function edge_transform_metadata_from_context(coord::Coordinate,
                                              target_frame::CoordinateFrame,
                                              context)
    has_spacecraft_reference = (coord.frame.axes isa LVLH) || (target_frame.axes isa LVLH)
    return EdgeTransformMetadata(
        true,
        true,
        has_spacecraft_reference,
    )
end

function required_edge_transform_metadata(edge::EdgeTransform)
    if edge isa STBTransform
        return EdgeTransformMetadata(
            edge.eop_type == :IAU1980,
            edge.eop_type == :IAU2000A,
            false,
        )
    elseif edge isa OrbitRelativeTransform
        return EdgeTransformMetadata(false, false, true)
    else
        return EdgeTransformMetadata(false, false, false)
    end
end

function metadata_satisfies(available::EdgeTransformMetadata,
                            required::EdgeTransformMetadata)
    (!required.has_eop_iau1980 || available.has_eop_iau1980) &&
    (!required.has_eop_iau2000a || available.has_eop_iau2000a) &&
    (!required.has_spacecraft_reference || available.has_spacecraft_reference)
end

function combine_edge_transform_metadata(a::EdgeTransformMetadata,
                                         b::EdgeTransformMetadata)
    return EdgeTransformMetadata(
        a.has_eop_iau1980 || b.has_eop_iau1980,
        a.has_eop_iau2000a || b.has_eop_iau2000a,
        a.has_spacecraft_reference || b.has_spacecraft_reference,
    )
end

function neighbors_for_axes(from_axes::Type,
                            available_metadata::EdgeTransformMetadata)
    lock(GRAPH_LOCK)
    try
        neighbors = DataType[]
        for ((from, to), edge) in EDGE_TRANSFORMS
            if from == from_axes
                req = required_edge_transform_metadata(edge)
                if metadata_satisfies(available_metadata, req)
                    push!(neighbors, to)
                end
            end
        end
        return neighbors
    finally
        unlock(GRAPH_LOCK)
    end
end

function path_edge_transform_metadata(path::Vector{DataType})
    meta = EdgeTransformMetadata(false, false, false)
    lock(GRAPH_LOCK)
    try
        for i in 1:(length(path)-1)
            edge = EDGE_TRANSFORMS[(path[i], path[i+1])]
            meta = combine_edge_transform_metadata(meta, required_edge_transform_metadata(edge))
        end
    finally
        unlock(GRAPH_LOCK)
    end
    return meta
end

"""
    TransformContext

Holds runtime data needed by various transform types.
"""
mutable struct TransformContext
    eop_iau1980::Union{Nothing, EopIau1980}
    eop_iau2000a::Union{Nothing, EopIau2000A}
    orbit_state::Union{Nothing, Coordinate}  # For orbit-relative frames
    bodies::Dict{Symbol, Any}                 # For planetary frames (future)
    model_family::Symbol                       # :auto, :fk5, :iau2006
    warn_on_ambiguity::Bool
end

function TransformContext(; eop_iau1980::Union{Nothing,EopIau1980}=nothing,
                            eop_iau2000a::Union{Nothing,EopIau2000A}=nothing,
                            orbit_state::Union{Nothing,Coordinate}=nothing,
                            bodies::Dict{Symbol,Any}=Dict{Symbol,Any}(),
                            model_family::Symbol=:auto,
                            warn_on_ambiguity::Bool=true)
    if !(model_family in (:auto, :fk5, :iau2006))
        error("TransformContext model_family must be one of :auto, :fk5, :iau2006")
    end
    return TransformContext(eop_iau1980, eop_iau2000a, orbit_state, bodies, model_family, warn_on_ambiguity)
end

const SHARED_DEFAULT_CONTEXT = TransformContext()

const DEFAULT_EOP_LOCK = ReentrantLock()
const DEFAULT_EOP_IAU1980 = Ref{Union{Nothing, EopIau1980}}(nothing)
const DEFAULT_EOP_IAU2000A = Ref{Union{Nothing, EopIau2000A}}(nothing)

function get_default_eop_iau1980()
    lock(DEFAULT_EOP_LOCK)
    try
        if DEFAULT_EOP_IAU1980[] === nothing
            DEFAULT_EOP_IAU1980[] = fetch_iers_eop(Val(:IAU1980))
        end
        return DEFAULT_EOP_IAU1980[]
    finally
        unlock(DEFAULT_EOP_LOCK)
    end
end

function get_default_eop_iau2000a()
    lock(DEFAULT_EOP_LOCK)
    try
        if DEFAULT_EOP_IAU2000A[] === nothing
            DEFAULT_EOP_IAU2000A[] = fetch_iers_eop(Val(:IAU2000A))
        end
        return DEFAULT_EOP_IAU2000A[]
    finally
        unlock(DEFAULT_EOP_LOCK)
    end
end

function ensure_context_eop!(context::TransformContext, eop_type::Symbol)
    if eop_type == :IAU1980
        if context.eop_iau1980 === nothing
            context.eop_iau1980 = get_default_eop_iau1980()
        end
        return context.eop_iau1980
    elseif eop_type == :IAU2000A
        if context.eop_iau2000a === nothing
            context.eop_iau2000a = get_default_eop_iau2000a()
        end
        return context.eop_iau2000a
    else
        return nothing
    end
end

function resolve_transform_context(context::Union{Nothing, TransformContext})
    return context === nothing ? SHARED_DEFAULT_CONTEXT : context
end

const AMBIGUITY_WARN_LOCK = ReentrantLock()
const AMBIGUITY_WARNED_KEYS = Set{Tuple{DataType,DataType,Symbol}}()

@inline function warn_ambiguous_model_once(from_axes::Type, to_axes::Type,
                                           model_family::Symbol, message::String)
    key = (from_axes, to_axes, model_family)
    should_warn = false
    lock(AMBIGUITY_WARN_LOCK)
    try
        if !(key in AMBIGUITY_WARNED_KEYS)
            push!(AMBIGUITY_WARNED_KEYS, key)
            should_warn = true
        end
    finally
        unlock(AMBIGUITY_WARN_LOCK)
    end

    if should_warn
        @warn message
    end
end

# =============================================================================
# Spacecraft State Extraction
# =============================================================================

"""
    spacecraft_to_coordinate(sc::Spacecraft, frame::CoordinateFrame)

Extract position and velocity from Spacecraft and convert to Coordinate.
Assumes spacecraft state is in an inertial frame (ICRF, GCRF, J2000).
"""
function spacecraft_to_coordinate(sc::Spacecraft, frame::CoordinateFrame)
    # Get Cartesian state from spacecraft
    cart_state = get_state(sc, Cartesian())
    state_vec = to_vector(cart_state)  # Returns 6-element vector [x,y,z,vx,vy,vz]
    
    pos = SVector{3,Float64}(state_vec[1:3])
    vel = SVector{3,Float64}(state_vec[4:6])
    
    return Coordinate(pos, vel, sc.time, frame)
end

# =============================================================================
# Coordinate ↔ OrbitStateVector Conversion
# =============================================================================

"""
Convert Epicycle Coordinate to STB OrbitStateVector.
"""
function to_orbit_state_vector(coord::Coordinate{T}) where T
    # STB uses Julian date (Float64) for time
    time_utc = coord.time.utc
    jd_utc = time_utc.jd
    
    return OrbitStateVector(
        jd_utc,
        coord.pos,
        coord.vel
    )
end

"""
Convert STB OrbitStateVector back to Epicycle Coordinate.
Keep original time and frame from template.
"""
function from_orbit_state_vector(sv::OrbitStateVector, time::Time, frame::CoordinateFrame)
    T = eltype(sv.r)
    return Coordinate(
        sv.r,
        sv.v,
        SVector{3,T}(0, 0, 0),  # STB doesn't provide acceleration
        time,
        frame
    )
end

# =============================================================================
# Transform Execution (dispatch on EdgeTransform type)
# =============================================================================

"""
Execute STB transformation.
"""
function execute_transform(coord::Coordinate, edge::STBTransform, 
                          target_frame::CoordinateFrame, context::TransformContext)
    # Get EOP data if needed
    eop = if edge.eop_type == :IAU1980
        ensure_context_eop!(context, :IAU1980)
    elseif edge.eop_type == :IAU2000A
        ensure_context_eop!(context, :IAU2000A)
    else
        nothing
    end
    
    # Check if EOP data is required but missing
    if edge.eop_type != :none && eop === nothing
        error("Transformation from $(edge.val_from) to $(edge.val_to) requires EOP data ($(edge.eop_type)), but none provided in context")
    end
    
    # Convert to OrbitStateVector
    sv = to_orbit_state_vector(coord)
    
    # Get Julian date for STB (UT1 time scale)
    jd_ut1 = coord.time.utc.jd  # Simplified - should apply UT1-UTC correction
    
    # Call STB function
    if eop !== nothing
        sv_new = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1, eop)
    else
        sv_new = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1)
    end
    
    # Convert back to Coordinate
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_itrf_to_j2000_fast(coord::Coordinate,
                                            target_frame::CoordinateFrame,
                                            context::TransformContext)
    eop = ensure_context_eop!(context, :IAU1980)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_ecef_to_eci(sv, Val(:ITRF), Val(:J2000), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_j2000_to_itrf_fast(coord::Coordinate,
                                            target_frame::CoordinateFrame,
                                            context::TransformContext)
    eop = ensure_context_eop!(context, :IAU1980)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_eci_to_ecef(sv, Val(:J2000), Val(:ITRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_itrf_to_gcrf_fast(coord::Coordinate,
                                           target_frame::CoordinateFrame,
                                           context::TransformContext)
    eop = ensure_context_eop!(context, :IAU1980)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_ecef_to_eci(sv, Val(:ITRF), Val(:GCRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_itrf_to_gcrf_iau_fast(coord::Coordinate,
                                               target_frame::CoordinateFrame,
                                               context::TransformContext)
    eop = ensure_context_eop!(context, :IAU2000A)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_ecef_to_eci(sv, Val(:ITRF), Val(:GCRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_gcrf_to_itrf_fast(coord::Coordinate,
                                           target_frame::CoordinateFrame,
                                           context::TransformContext)
    eop = ensure_context_eop!(context, :IAU1980)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_eci_to_ecef(sv, Val(:GCRF), Val(:ITRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_gcrf_to_itrf_iau_fast(coord::Coordinate,
                                               target_frame::CoordinateFrame,
                                               context::TransformContext)
    eop = ensure_context_eop!(context, :IAU2000A)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_eci_to_ecef(sv, Val(:GCRF), Val(:ITRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_itrf_to_cirs_fast(coord::Coordinate,
                                           target_frame::CoordinateFrame,
                                           context::TransformContext)
    eop = ensure_context_eop!(context, :IAU2000A)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_ecef_to_eci(sv, Val(:ITRF), Val(:CIRS), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

@inline function execute_cirs_to_itrf_fast(coord::Coordinate,
                                           target_frame::CoordinateFrame,
                                           context::TransformContext)
    eop = ensure_context_eop!(context, :IAU2000A)
    jd_ut1 = coord.time.utc.jd
    sv = OrbitStateVector(jd_ut1, coord.pos, coord.vel)
    sv_new = sv_eci_to_ecef(sv, Val(:CIRS), Val(:ITRF), jd_ut1, eop)
    return from_orbit_state_vector(sv_new, coord.time, target_frame)
end

"""
Execute a sequence of STB edges with one Coordinate->OrbitStateVector conversion
and one OrbitStateVector->Coordinate conversion.
"""
function execute_stb_sequence(coord::Coordinate,
                              edges::Vector{EdgeTransform},
                              target_frame::CoordinateFrame,
                              context::TransformContext)
    sv = to_orbit_state_vector(coord)
    jd_ut1 = coord.time.utc.jd

    for edge_any in edges
        edge = edge_any::STBTransform
        eop = if edge.eop_type == :IAU1980
            ensure_context_eop!(context, :IAU1980)
        elseif edge.eop_type == :IAU2000A
            ensure_context_eop!(context, :IAU2000A)
        else
            nothing
        end

        if edge.eop_type != :none && eop === nothing
            error("Transformation from $(edge.val_from) to $(edge.val_to) requires EOP data ($(edge.eop_type)), but none provided in context")
        end

        if eop === nothing
            sv = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1)
        else
            sv = edge.stb_function(sv, edge.val_from, edge.val_to, jd_ut1, eop)
        end
    end

    return from_orbit_state_vector(sv, coord.time, target_frame)
end

"""
Execute custom transformation.
"""
function execute_transform(coord::Coordinate, edge::CustomTransform,
                          target_frame::CoordinateFrame, context::TransformContext)
    return edge.transform_fn(coord, target_frame)
end

"""
Execute orbit-relative transformation.
Extracts orbit state from LVLH spacecraft reference.
"""
function execute_transform(coord::Coordinate, edge::OrbitRelativeTransform,
                          target_frame::CoordinateFrame, context::TransformContext)
    # Extract spacecraft reference from LVLH frame
    spacecraft = if target_frame.axes isa LVLH
        target_frame.axes.reference
    elseif coord.frame.axes isa LVLH
        coord.frame.axes.reference
    else
        error("OrbitRelativeTransform requires LVLH frame with Spacecraft reference")
    end
    
    # Convert spacecraft state to Coordinate (assumes inertial frame)
    # TODO: Handle frame conversion if spacecraft.coord_sys != inertial
    orbit_state = spacecraft_to_coordinate(spacecraft, coord.frame)
    
    return edge.transform_fn(coord, orbit_state, target_frame)
end

# =============================================================================
# Path Finding (BFS)
# =============================================================================

"""
    find_axes_path(from_axes::Type, to_axes::Type) -> Vector{DataType}

Find shortest path through transformation graph using BFS.
Returns sequence of axes types from source to target.
"""
function find_axes_path(from_axes::Type, to_axes::Type)
    # Direct connection?
    if haskey(EDGE_TRANSFORMS, (from_axes, to_axes))
        return [from_axes, to_axes]
    end
    
    # BFS
    queue = [(from_axes, [from_axes])]
    visited = Set{DataType}([from_axes])
    
    while !isempty(queue)
        (current, path) = popfirst!(queue)
        
        # Get neighbors from registered edge registry
        neighbors = neighbors_for_axes(current, EdgeTransformMetadata())
        
        for neighbor in neighbors
            if neighbor == to_axes
                return [path; neighbor]
            end
            
            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, (neighbor, [path; neighbor]))
            end
        end
    end
    
    error("No transformation path found from $(from_axes) to $(to_axes)")
end

function find_axes_path(from_axes::Type, to_axes::Type,
                        available_metadata::EdgeTransformMetadata)
    if haskey(EDGE_TRANSFORMS, (from_axes, to_axes))
        direct_edge = EDGE_TRANSFORMS[(from_axes, to_axes)]
        req = required_edge_transform_metadata(direct_edge)
        if metadata_satisfies(available_metadata, req)
            return [from_axes, to_axes]
        end
    end

    queue = [(from_axes, [from_axes])]
    visited = Set{DataType}([from_axes])

    while !isempty(queue)
        (current, path) = popfirst!(queue)
        neighbors = neighbors_for_axes(current, available_metadata)

        for neighbor in neighbors
            if neighbor == to_axes
                return [path; neighbor]
            end

            if neighbor ∉ visited
                push!(visited, neighbor)
                push!(queue, (neighbor, [path; neighbor]))
            end
        end
    end

    error("No transformation path found from $(from_axes) to $(to_axes) for metadata $(available_metadata)")
end

"""
    get_cached_path(from_axes::Type, to_axes::Type) -> Vector{DataType}

Get transformation path, using cache if available.
"""
function get_cached_path(from_axes::Type, to_axes::Type)
    return find_axes_path(from_axes, to_axes, EdgeTransformMetadata())
end

function get_precomputed_path(from_axes::Type, to_axes::Type,
                              available_metadata::EdgeTransformMetadata)
    current_version = graph_version()

    # Fast-path precedence
    fast_key = (from_axes, to_axes)
    lock(GRAPH_LOCK)
    try
        if haskey(FAST_PATH_IMPLEMENTATIONS, fast_key)
            edge = FAST_PATH_IMPLEMENTATIONS[fast_key]
            req = required_edge_transform_metadata(edge)
            if metadata_satisfies(available_metadata, req)
                path = [from_axes, to_axes]
                return PrecomputedPathEntry(path, req, current_version)
            end
        end
    finally
        unlock(GRAPH_LOCK)
    end

    key = PathCacheKey(from_axes, to_axes, available_metadata, current_version)

    lock(CACHE_LOCK)
    try
        if haskey(PRECOMPUTED_PATH_CACHE, key)
            return PRECOMPUTED_PATH_CACHE[key]
        end
    finally
        unlock(CACHE_LOCK)
    end

    path = find_axes_path(from_axes, to_axes, available_metadata)
    metadata_snapshot = path_edge_transform_metadata(path)
    entry = PrecomputedPathEntry(path, metadata_snapshot, current_version)

    lock(CACHE_LOCK)
    try
        PRECOMPUTED_PATH_CACHE[key] = entry
    finally
        unlock(CACHE_LOCK)
    end

    return entry
end

function get_sim_plan(from_axes::Type, to_axes::Type,
                      available_metadata::EdgeTransformMetadata)
    current_version = graph_version()
    key = PathCacheKey(from_axes, to_axes, available_metadata, current_version)

    lock(CACHE_LOCK)
    try
        if haskey(SIM_PLAN_CACHE, key)
            return SIM_PLAN_CACHE[key]
        end
    finally
        unlock(CACHE_LOCK)
    end

    precomputed = get_precomputed_path(from_axes, to_axes, available_metadata)
    edges = EdgeTransform[]
    lock(GRAPH_LOCK)
    try
        for i in 1:(length(precomputed.path)-1)
            push!(edges, EDGE_TRANSFORMS[(precomputed.path[i], precomputed.path[i+1])])
        end
    finally
        unlock(GRAPH_LOCK)
    end

    plan = SimPlanEntry(
        precomputed.path,
        edges,
        precomputed.edge_transform_metadata,
        precomputed.graph_version,
    )

    lock(CACHE_LOCK)
    try
        SIM_PLAN_CACHE[key] = plan
    finally
        unlock(CACHE_LOCK)
    end

    return plan
end

# =============================================================================
# High-Level Transform Function
# =============================================================================

"""
    transform(coord::Coordinate, target_frame::CoordinateFrame, [context])

Transform coordinate to target frame using graph-based routing.
Runtime context is resolved lazily when omitted.
"""
function transform(coord::Coordinate, target_frame::CoordinateFrame)
    return transform(coord, target_frame, SHARED_DEFAULT_CONTEXT)
end

function transform(coord::Coordinate, target_frame::CoordinateFrame, ::Nothing)
    return transform(coord, target_frame, SHARED_DEFAULT_CONTEXT)
end

function transform(coord::Coordinate, target_frame::CoordinateFrame,
                   context::TransformContext)
    context_resolved = context

    # Extract axes types
    from_axes = typeof(coord.frame.axes)
    to_axes = typeof(target_frame.axes)
    
    # Same frame? Nothing to do
    if from_axes == to_axes && coord.frame.origin == target_frame.origin
        return coord
    end

    same_origin = (coord.frame.origin == target_frame.origin)

    model_family = context_resolved.model_family

    if (from_axes == ITRF && to_axes == GCRF) || (from_axes == GCRF && to_axes == ITRF)
        if model_family == :auto && context_resolved.warn_on_ambiguity
            warn_ambiguous_model_once(
                from_axes,
                to_axes,
                model_family,
                "Ambiguous model-family transform for $(from_axes) -> $(to_axes). " *
                "Defaulting to FK5 (:IAU1980). Pass TransformContext(model_family=:fk5|:iau2006) to disambiguate."
            )
        end
    end

    if (from_axes == J2000 || to_axes == J2000) && model_family == :iau2006 && context_resolved.warn_on_ambiguity
        warn_ambiguous_model_once(
            from_axes,
            to_axes,
            model_family,
            "J2000 is treated as FK5-specific in this prototype. Requested model_family=:iau2006 is ignored for this pair."
        )
    end

    if same_origin && from_axes == ITRF && to_axes == J2000
        return execute_itrf_to_j2000_fast(coord, target_frame, context_resolved)
    elseif same_origin && from_axes == J2000 && to_axes == ITRF
        return execute_j2000_to_itrf_fast(coord, target_frame, context_resolved)
    elseif same_origin && from_axes == ITRF && to_axes == GCRF
        if model_family == :iau2006
            return execute_itrf_to_gcrf_iau_fast(coord, target_frame, context_resolved)
        else
            return execute_itrf_to_gcrf_fast(coord, target_frame, context_resolved)
        end
    elseif same_origin && from_axes == GCRF && to_axes == ITRF
        if model_family == :iau2006
            return execute_gcrf_to_itrf_iau_fast(coord, target_frame, context_resolved)
        else
            return execute_gcrf_to_itrf_fast(coord, target_frame, context_resolved)
        end
    elseif same_origin && from_axes == ITRF && to_axes == CIRS
        return execute_itrf_to_cirs_fast(coord, target_frame, context_resolved)
    elseif same_origin && from_axes == CIRS && to_axes == ITRF
        return execute_cirs_to_itrf_fast(coord, target_frame, context_resolved)
    end

    # Direct STB fast path: bypass planner/cache/graph for common one-edge transforms
    # (ITRF↔GCRF, ITRF↔J2000, etc.)
    fast_key = (from_axes, to_axes)
    if haskey(FAST_PATH_IMPLEMENTATIONS, fast_key)
        fast_edge = FAST_PATH_IMPLEMENTATIONS[fast_key]
        if fast_edge isa STBTransform
            return execute_transform(coord, fast_edge, target_frame, context_resolved)
        end
    end
    
    # Build edge transform metadata from runtime context and get compiled plan
    available_metadata = edge_transform_metadata_from_context(coord, target_frame, context_resolved)
    plan = get_sim_plan(from_axes, to_axes, available_metadata)
    path = plan.path

    # Fast execution for pure STB plans (common Earth-fixed/inertial chains):
    # avoid per-edge Coordinate<->OrbitStateVector conversions.
    if !isempty(plan.edges) && all(edge -> edge isa STBTransform, plan.edges)
        return execute_stb_sequence(coord, plan.edges, target_frame, context_resolved)
    end
    
    # Execute transformations along path
    current_coord = coord
    for i in 1:(length(path)-1)
        from = path[i]
        to = path[i+1]
        
        # Get edge transformation from compiled plan
        edge = plan.edges[i]
        
        # Create intermediate frame
        # Special handling for LVLH which requires a reference
        intermediate_axes = if to == LVLH
            # If transforming TO LVLH, use the target frame's LVLH (which has the reference)
            target_frame.axes
        elseif from == LVLH
            # If transforming FROM LVLH, use the source frame's LVLH
            coord.frame.axes
        else
            # For other axes types, construct with zero arguments
            to()
        end
        
        intermediate_frame = CoordinateFrame(target_frame.origin, intermediate_axes)
        
        # Execute transformation
        current_coord = execute_transform(current_coord, edge, intermediate_frame, context_resolved)
    end
    
    return current_coord
end

# Test/demo functions moved to simple_transform_test.jl


