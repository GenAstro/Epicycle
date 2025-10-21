# Copyright 2025 Gen Astro LLC. All Rights Reserved.
#
# This software is licensed under the GNU AGPL v3.0,
# WITHOUT ANY WARRANTY, including implied warranties of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# This file may also be used under a commercial license,
# if one has been purchased from Gen Astro LLC.
#
# By modifying this software, you agree to the terms of the
# Gen Astro LLC Contributor License Agreement.

"""
    OrbitState(state::AbstractVector, statetype::AbstractOrbitStateType)

A generic wrapper for orbital state vectors and their associated metadata.

# Arguments
- `state`: The numerical state vector (e.g., position/velocity, orbital elements). Must be a vector of real numbers.
- `statetype`: Marker instance indicating the state representation (e.g., `Cartesian()`, `Keplerian()`, etc.).
-              - to see all state types, use `subtypes(AbstractOrbitStateType)`.
# Example
```julia
state_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]  # Cartesian position and velocity
os = OrbitState(state_vec, Cartesian())
```

# Notes
- Allows easy switching between different state representations without type instability.
- The struct is internally parameterized for performance, and type safety, and differentiation, but users should construct it using the outer constructor as shown above.
"""
struct OrbitState{S<:AbstractOrbitStateType, T<:Real, V<:AbstractVector{T}}
    state::V
    statetype::S
end

# TODO Add Delaunay, that one was missed. 

# Show method for OrbitState
function show(io::IO, os::OrbitState)
    println(io, "OrbitState:")
    println(io, "  statetype: $(typeof(os.statetype))")
    state_obj = _construct_state(os.state, os.statetype)
    show(io, state_obj)
end

# Ensure snapshots don't alias the inner vector
function Base.deepcopy_internal(os::OrbitState, dict::IdDict)
    st_copy = Base.deepcopy_internal(getfield(os, :state), dict)  # copy 6-vector
    tag     = getfield(os, :statetype)                            # tags are immutable
    return OrbitState(st_copy, tag)
end

# Optional: make shallow copy non-aliasing too
Base.copy(os::OrbitState) = OrbitState(copy(getfield(os, :state)), getfield(os, :statetype))

# Helper: construct concrete state from vector and marker
function _construct_state(vec, ::Cartesian)
    CartesianState(vec)
end
function _construct_state(vec, ::Keplerian)
    KeplerianState(vec...)
end
function _construct_state(vec, ::Equinoctial)
    EquinoctialState(vec...)
end
function _construct_state(vec, ::SphericalRADEC)
    SphericalRADECState(vec...)
end
function _construct_state(vec, ::SphericalAZIFPA)
    SphericalAZIFPAState(vec...)
end
function _construct_state(vec, ::ModifiedEquinoctial)
    ModifiedEquinoctialState(vec...)
end
function _construct_state(vec, ::OutGoingAsymptote)
    OutGoingAsymptoteState(vec...)
end
function _construct_state(vec, ::IncomingAsymptote)
    IncomingAsymptoteState(vec...)
end
function _construct_state(vec, ::ModifiedKeplerian)
    ModifiedKeplerianState(vec...)
end

function _construct_state(vec, ::AlternateEquinoctial)
    AlternateEquinoctialState(vec...)
end

# Helper: get state vector from concrete state
_state_to_vector(cs::CartesianState) = cs.posvel
_state_to_vector(ks::KeplerianState) = [ks.sma, ks.ecc, ks.inc, ks.raan, ks.aop, ks.ta]
_state_to_vector(s::SphericalRADECState) = [s.r, s.dec, s.ra, s.v, s.decv, s.rav]
_state_to_vector(s::SphericalAZIFPAState) = [s.r, s.ra, s.dec, s.v, s.vazi, s.fpa]
_state_to_vector(ms::ModifiedEquinoctialState) = [ms.p, ms.f, ms.g, ms.h, ms.k, ms.L]
_state_to_vector(os::OutGoingAsymptoteState) = [os.rp, os.c3, os.rla, os.dla, os.bpa, os.ta]
_state_to_vector(is::IncomingAsymptoteState) = [is.rp, is.c3, is.rla, is.dla, is.bpa, is.ta]
_state_to_vector(ms::ModifiedKeplerianState) = [ms.rp, ms.ra, ms.inc, ms.raan, ms.aop, ms.ta]
_state_to_vector(es::EquinoctialState) = [es.a, es.h, es.k, es.p, es.q, es.mlong]
_state_to_vector(ae::AlternateEquinoctialState) = [ae.a, ae.h, ae.k, ae.altp, ae.altq, ae.mlong]

"""
    OrbitState(state_obj::AbstractState)

Construct an `OrbitState` from a concrete state type (e.g., `CartesianState`, `KeplerianState`, etc.).

# Arguments
- `state_obj::AbstractState`: The concrete state instance to wrap (such as `CartesianState`, `KeplerianState`, etc.).

# Returns
- An `OrbitState` instance containing the state vector and a marker indicating the state type.

# Example
```julia
cs = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0])
os = OrbitState(cs)  # os is an OrbitState with statetype Cartesian

ks = KeplerianState(8000.0, 0.05, 0.1, 0.2, 0.3, 0.4)
os_kep = OrbitState(ks)  # os_kep is an OrbitState with statetype Keplerian
```

# Notes
- This constructor extracts the state vector and the marker type from the concrete state and stores them in the `OrbitState` wrapper.
- Use this to enable generic state management and easy switching between different state representations.
"""
#function OrbitState(state_obj::AbstractOrbitState)
#    marker = _marker_type(state_obj)
#    vec = _state_to_vector(state_obj)
#    OrbitState(vec, marker)
#end

function OrbitState(state_obj::AbstractOrbitState)
    marker = state_type_to_tag(state_obj)  
    vec = _state_to_vector(state_obj)
    OrbitState(vec, marker)
end
"""
    to_state(os::OrbitState) -> AbstractState

Convert an `OrbitState` to its corresponding concrete state type (subtype of `AbstractState`), using the stored state vector and marker type.

# Arguments
- `os::OrbitState`: The `OrbitState` instance to convert.

# Returns
- The concrete state type (e.g., `CartesianState`, `KeplerianState`, etc.) corresponding to the `statetype` marker and the stored state vector.

# Example
```julia
state_vec = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
os = OrbitState(state_vec, Cartesian())
cs = to_state(os)  # cs is a CartesianState

# For other state types:
os_kep = OrbitState([8000.0, 0.05, 0.1, 0.2, 0.3, 0.4], Keplerian())
ks = to_state(os_kep)  # ks is a KeplerianState
```

# Notes
- This function does not perform any physical conversions; it simply reconstructs the concrete state type from the stored vector and marker.
- For conversions between different state types (e.g., `CartesianState` to `KeplerianState`), use the appropriate constructor or conversion function, providing additional parameters (such as `μ`) if required.
"""
to_state(os::OrbitState) = _construct_state(os.state, os.statetype)
