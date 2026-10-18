# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Nonsingular Keplerian elements: a; e cos ϖ; e sin ϖ; sin(i/2) sin Ω; sin(i/2) cos Ω; λ, with
# ϖ = Ω + ω the longitude of periapsis and λ = Ω + ω + M the mean longitude. Every element stays
# defined at e = 0 and i = 0, where Ω, ω and M separately do not. Converts through KeplerianState.
# =============================================================================

"""
    NonsingularKeplerian()

The tag for the `NonsingularKeplerianState` representation. An `OrbitState` carries one to say what
its six numbers are, and [`state_tag_to_type`](@ref) maps it to `NonsingularKeplerianState`.

# Example
```julia
state_tag_to_type(NonsingularKeplerian())     # NonsingularKeplerianState
```
"""
struct NonsingularKeplerian <: AbstractOrbitStateType end

"""
    NonsingularKeplerianState(sma, e1, e2, e3, e4, e5)

Nonsingular Keplerian elements of a closed orbit, defined at zero eccentricity and zero
inclination.

# Fields
- `sma`: semi-major axis a, km.
- `e1`: e cos(Ω + ω).
- `e2`: e sin(Ω + ω).
- `e3`: sin(i/2) sin Ω.
- `e4`: sin(i/2) cos Ω.
- `e5`: mean longitude Ω + ω + M, rad, in [0, 2π).

# Notes
- Defined for closed orbits, 0 ≤ e < 1. A hyperbola converts to a state of `NaN`s with a warning,
  as other singular conversions do.
- (e1, e2) and (e3, e4) carry the eccentricity and inclination with the direction of periapsis and
  of the node; at e = 0 or i = 0 the pair is zero and the direction is not needed.
- Conversions to the other representations go through `KeplerianState`, so a circular or
  equatorial orbit reads back with the Keplerian conventions for the angles undefined there.

# Examples
```julia
ns = NonsingularKeplerianState(KeplerianState(7000.0, 0.01, 0.5, 1.0, 2.0, 0.3), 398600.4418)
CartesianState(ns, 398600.4418)
```
"""
struct NonsingularKeplerianState{T<:Real} <: AbstractOrbitState
    sma::T
    e1::T
    e2::T
    e3::T
    e4::T
    e5::T
end

function show(io::IO, s::NonsingularKeplerianState)
    println(io, "NonsingularKeplerianState:")
    println(io, @sprintf("  sma        = %14.6f", s.sma))
    println(io, @sprintf("  e1         = %14.8f", s.e1))
    println(io, @sprintf("  e2         = %14.8f", s.e2))
    println(io, @sprintf("  e3         = %14.8f", s.e3))
    println(io, @sprintf("  e4         = %14.8f", s.e4))
    println(io, @sprintf("  e5   (deg) = %14.6f", rad2deg(s.e5)))
end

"""
    kep_to_nonsingular(kep::AbstractVector{<:Real}) -> Vector

Keplerian elements `[a, e, i, Ω, ω, ν]` to nonsingular Keplerian elements
`[a, e1, e2, e3, e4, e5]`. See [`NonsingularKeplerianState`](@ref).

# Notes
Defined for a closed orbit; for e ≥ 1 or a ≤ 0 it warns and returns `NaN`s.

# Returns
A 6-element vector.

# Examples
```julia
kep_to_nonsingular([7000.0, 0.01, 0.5, 1.0, 2.0, 0.3])
```
"""
function kep_to_nonsingular(kep::AbstractVector{<:Real})
    a, e, i, Ω, ω, ν = kep
    T = promote_type(eltype(kep), Float64)
    if !(0 <= e < 1) || a <= 0
        @warn "Conversion Failed: nonsingular Keplerian elements are defined for a closed orbit; got e = $e, a = $a."
        return fill(T(NaN), 6)
    end
    ϖ = Ω + ω
    M = _true_to_mean_anom(ν, e)
    return T[a, e * cos(ϖ), e * sin(ϖ), sin(i / 2) * sin(Ω), sin(i / 2) * cos(Ω), _wrap_2pi(ϖ + M)]
end

"""
    nonsingular_to_kep(ns::AbstractVector{<:Real}) -> Vector

Nonsingular Keplerian elements `[a, e1, e2, e3, e4, e5]` to Keplerian elements
`[a, e, i, Ω, ω, ν]`, angles in [0, 2π). See [`NonsingularKeplerianState`](@ref).

# Notes
At i = 0 the node is undefined and Ω is returned as 0; at e = 0 the periapsis is undefined and
ω is returned as −Ω, wrapped, so the argument of latitude carries the position. Elements that
describe no closed orbit (e ≥ 1 or sin(i/2) > 1) warn and return `NaN`s.

# Returns
A 6-element vector.

# Examples
```julia
nonsingular_to_kep(kep_to_nonsingular([7000.0, 0.01, 0.5, 1.0, 2.0, 0.3]))
```
"""
function nonsingular_to_kep(ns::AbstractVector{<:Real})
    a, e1, e2, e3, e4, e5 = ns
    T = promote_type(eltype(ns), Float64)
    e = sqrt(e1^2 + e2^2)
    s = sqrt(e3^2 + e4^2)
    if !(e < 1) || s > 1 || a <= 0
        @warn "Conversion Failed: nonsingular elements describe no closed orbit; got e = $e, sin(i/2) = $s, a = $a."
        return fill(T(NaN), 6)
    end
    ϖ = e > 0 ? atan(e2, e1) : zero(T)
    Ω = s > 0 ? atan(e3, e4) : zero(T)
    i = 2 * asin(s)
    E = _eccentric_anomaly(e5 - ϖ, e)
    ν = atan(sqrt(1 - e^2) * sin(E), cos(E) - e)
    return T[a, e, i, _wrap_2pi(Ω), _wrap_2pi(ϖ - Ω), _wrap_2pi(ν)]
end

# --- construction and conversion ------------------------------------------------

NonsingularKeplerianState(s::AbstractVector{T}, μ::Real) where {T<:Real} =
    (length(s) == 6 || throw(ArgumentError("A state vector must have 6 elements; got $(length(s)).")); NonsingularKeplerianState{T}(s[1], s[2], s[3], s[4], s[5], s[6]))
NonsingularKeplerianState(s::AbstractVector{T}) where {T<:Real} = NonsingularKeplerianState(s, zero(T))
NonsingularKeplerianState(a::Real, b::Real, c::Real, d::Real, e::Real, f::Real) =
    NonsingularKeplerianState(promote(a, b, c, d, e, f)...)
NonsingularKeplerianState(s::NonsingularKeplerianState, μ::Real) = s
NonsingularKeplerianState(s::NonsingularKeplerianState) = s

NonsingularKeplerianState(kep::KeplerianState, μ::Real) =
    NonsingularKeplerianState(kep_to_nonsingular(to_vector(kep)))
KeplerianState(ns::NonsingularKeplerianState, μ::Real) =
    KeplerianState(nonsingular_to_kep(to_vector(ns)), μ)
NonsingularKeplerianState(s::AbstractOrbitState, μ::Real) =
    NonsingularKeplerianState(KeplerianState(s, μ), μ)

to_vector(s::NonsingularKeplerianState) = [s.sma, s.e1, s.e2, s.e3, s.e4, s.e5]
_state_to_vector(s::NonsingularKeplerianState) = to_vector(s)
_construct_state(vec, ::NonsingularKeplerian) = NonsingularKeplerianState(vec...)
@inline state_tag_to_type(::NonsingularKeplerian) = NonsingularKeplerianState
@inline state_type_to_tag(::Type{<:NonsingularKeplerianState}) = NonsingularKeplerian()
