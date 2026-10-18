# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Delaunay elements: the angles l = M, g = ω, h = Ω and their conjugate momenta L = √(μa),
# G = L√(1 − e²), H = G cos i, the canonical set of two-body perturbation theory. Converts through
# KeplerianState.
# =============================================================================

"""
    Delaunay()

The tag for the `DelaunayState` representation. An `OrbitState` carries one to say what
its six numbers are, and [`state_tag_to_type`](@ref) maps it to `DelaunayState`.

# Example
```julia
state_tag_to_type(Delaunay())     # DelaunayState
```
"""
struct Delaunay <: AbstractOrbitStateType end

"""
    DelaunayState(l, g, h, L, G, H)

Delaunay elements of a closed orbit: three angles and their conjugate momenta.

# Fields
- `l`: mean anomaly M, rad, in [0, 2π).
- `g`: argument of periapsis ω, rad.
- `h`: right ascension of the ascending node Ω, rad.
- `L`: √(μa), km²/s.
- `G`: L√(1 − e²), the specific angular momentum, km²/s.
- `H`: G cos i, the angular momentum about the coordinate system's Z axis, km²/s.

# Notes
- Defined for closed orbits, 0 ≤ e < 1. A hyperbola converts to a state of `NaN`s with a warning,
  as other singular conversions do.
- The momenta depend on μ, so a conversion to or from `DelaunayState` needs the gravitational
  parameter of the central body.
- Singular where ω and Ω are: at e = 0 and i = 0 the angles take the Keplerian conventions.

# Examples
```julia
d = DelaunayState(KeplerianState(7000.0, 0.01, 0.5, 1.0, 2.0, 0.3), 398600.4418)
KeplerianState(d, 398600.4418)
```
"""
struct DelaunayState{T<:Real} <: AbstractOrbitState
    l::T
    g::T
    h::T
    L::T
    G::T
    H::T
end

function show(io::IO, s::DelaunayState)
    println(io, "DelaunayState:")
    println(io, @sprintf("  l    (deg) = %14.6f", rad2deg(s.l)))
    println(io, @sprintf("  g    (deg) = %14.6f", rad2deg(s.g)))
    println(io, @sprintf("  h    (deg) = %14.6f", rad2deg(s.h)))
    println(io, @sprintf("  L          = %14.6f", s.L))
    println(io, @sprintf("  G          = %14.6f", s.G))
    println(io, @sprintf("  H          = %14.6f", s.H))
end

"""
    kep_to_delaunay(kep::AbstractVector{<:Real}, μ::Real) -> Vector

Keplerian elements `[a, e, i, Ω, ω, ν]` to Delaunay elements `[l, g, h, L, G, H]`. See
[`DelaunayState`](@ref).

# Notes
Defined for a closed orbit and positive μ; otherwise it warns and returns `NaN`s.

# Returns
A 6-element vector, angles in radians and momenta in km²/s.

# Examples
```julia
kep_to_delaunay([7000.0, 0.01, 0.5, 1.0, 2.0, 0.3], 398600.4418)
```
"""
function kep_to_delaunay(kep::AbstractVector{<:Real}, μ::Real)
    a, e, i, Ω, ω, ν = kep
    T = promote_type(eltype(kep), typeof(μ), Float64)
    if !(0 <= e < 1) || a <= 0 || μ <= 0
        @warn "Conversion Failed: Delaunay elements are defined for a closed orbit and positive μ; got e = $e, a = $a, μ = $μ."
        return fill(T(NaN), 6)
    end
    L = sqrt(μ * a)
    G = L * sqrt(1 - e^2)
    return T[_wrap_2pi(_true_to_mean_anom(ν, e)), ω, Ω, L, G, G * cos(i)]
end

"""
    delaunay_to_kep(d::AbstractVector{<:Real}, μ::Real) -> Vector

Delaunay elements `[l, g, h, L, G, H]` to Keplerian elements `[a, e, i, Ω, ω, ν]`. See
[`DelaunayState`](@ref).

# Notes
Elements that describe no orbit (L ≤ 0, |G| > L or |H| > G) or a non-positive μ warn and return
`NaN`s.

# Returns
A 6-element vector, true anomaly in [0, 2π).

# Examples
```julia
delaunay_to_kep(kep_to_delaunay([7000.0, 0.01, 0.5, 1.0, 2.0, 0.3], 398600.4418), 398600.4418)
```
"""
function delaunay_to_kep(d::AbstractVector{<:Real}, μ::Real)
    l, g, h, L, G, H = d
    T = promote_type(eltype(d), typeof(μ), Float64)
    if μ <= 0 || L <= 0 || !(0 < G <= L) || abs(H) > G
        @warn "Conversion Failed: Delaunay elements describe no orbit; got L = $L, G = $G, H = $H, μ = $μ."
        return fill(T(NaN), 6)
    end
    a = L^2 / μ
    e = sqrt(max(1 - (G / L)^2, zero(T)))
    i = acos(clamp(H / G, -one(T), one(T)))
    E = _eccentric_anomaly(l, e)
    ν = atan(sqrt(1 - e^2) * sin(E), cos(E) - e)
    return T[a, e, i, h, g, _wrap_2pi(ν)]
end

# --- construction and conversion ------------------------------------------------

DelaunayState(s::AbstractVector{T}, μ::Real) where {T<:Real} =
    (length(s) == 6 || throw(ArgumentError("A state vector must have 6 elements; got $(length(s)).")); DelaunayState{T}(s[1], s[2], s[3], s[4], s[5], s[6]))
DelaunayState(s::AbstractVector{T}) where {T<:Real} = DelaunayState(s, zero(T))
DelaunayState(a::Real, b::Real, c::Real, d::Real, e::Real, f::Real) =
    DelaunayState(promote(a, b, c, d, e, f)...)
DelaunayState(s::DelaunayState, μ::Real) = s
DelaunayState(s::DelaunayState) = s

DelaunayState(kep::KeplerianState, μ::Real) = DelaunayState(kep_to_delaunay(to_vector(kep), μ))
KeplerianState(d::DelaunayState, μ::Real) = KeplerianState(delaunay_to_kep(to_vector(d), μ), μ)
DelaunayState(s::AbstractOrbitState, μ::Real) = DelaunayState(KeplerianState(s, μ), μ)

to_vector(s::DelaunayState) = [s.l, s.g, s.h, s.L, s.G, s.H]
_state_to_vector(s::DelaunayState) = to_vector(s)
_construct_state(vec, ::Delaunay) = DelaunayState(vec...)
@inline state_tag_to_type(::Delaunay) = DelaunayState
@inline state_type_to_tag(::Type{<:DelaunayState}) = Delaunay()

# --- every other representation reaches these two through KeplerianState -------------

for S in (:CartesianState, :SphericalRADECState, :SphericalAZIFPAState, :ModifiedEquinoctialState,
          :OutGoingAsymptoteState, :IncomingAsymptoteState, :ModifiedKeplerianState,
          :EquinoctialState, :AlternateEquinoctialState, :BrouwerMeanLongState,
          :BrouwerMeanShortState),
    N in (:NonsingularKeplerianState, :DelaunayState)
    @eval $S(s::$N, μ::Real) = $S(KeplerianState(s, μ), μ)
end
NonsingularKeplerianState(d::DelaunayState, μ::Real) = NonsingularKeplerianState(KeplerianState(d, μ), μ)
DelaunayState(n::NonsingularKeplerianState, μ::Real) = DelaunayState(KeplerianState(n, μ), μ)
