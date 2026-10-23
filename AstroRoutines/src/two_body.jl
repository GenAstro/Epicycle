# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Two-body quantities of a single state: the position r and velocity v of a body relative to a
# central body of gravitational parameter μ, in any consistent units (km, km/s and km³/s² in
# Epicycle). Each is a closed form of (r, v, μ), differentiable with ForwardDiff, and throws an
# ArgumentError naming the orbit's energy or eccentricity where the quantity is undefined. The
# vector algebra is written out so the package keeps no dependencies.

function _check_state_vectors(r::AbstractVector, v::AbstractVector)
    length(r) == 3 || throw(ArgumentError("position r must have 3 components; got $(length(r))."))
    length(v) == 3 || throw(ArgumentError("velocity v must have 3 components; got $(length(v))."))
    return nothing
end

_dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_norm3(a) = sqrt(_dot3(a, a))
_cross3(a, b) = [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]

"""
    specific_energy(r::AbstractVector, v::AbstractVector, μ::Real)

Specific orbital energy of a two-body state, `ε = v²/2 − μ/r`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
The energy per unit mass [km²/s²]: negative on a closed orbit, zero on a parabola, positive on a
hyperbola.

# Example
```jldoctest
specific_energy([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

-28.209205204146848
```
"""
function specific_energy(r::AbstractVector, v::AbstractVector, μ::Real)
    _check_state_vectors(r, v)
    return _dot3(v, v) / 2 - μ / _norm3(r)
end

"""
    angular_momentum(r::AbstractVector, v::AbstractVector)

Specific angular momentum vector of a two-body state, `h = r × v`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].

# Returns
The 3-component vector `h` [km²/s], normal to the orbit plane.

# Example
```jldoctest
angular_momentum([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9])

# output

3-element Vector{Float64}:
 -8610.0
 -6780.0
 51920.0
```
"""
function angular_momentum(r::AbstractVector, v::AbstractVector)
    _check_state_vectors(r, v)
    return _cross3(r, v)
end

"""
    angular_momentum_magnitude(r::AbstractVector, v::AbstractVector)

Magnitude of the specific angular momentum of a two-body state, `‖r × v‖`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].

# Returns
`‖h‖` [km²/s].

# Example
```jldoctest
angular_momentum_magnitude([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9])

# output

53063.98873058828
```
"""
function angular_momentum_magnitude(r::AbstractVector, v::AbstractVector)
    _check_state_vectors(r, v)
    return _norm3(_cross3(r, v))
end

"""
    orbit_normal(r::AbstractVector, v::AbstractVector)

Unit normal to the orbit plane of a two-body state, along `r × v`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].

# Returns
The 3-component unit vector `ĥ`.

# Notes
Throws `ArgumentError` when `r × v = 0`, a radial trajectory, which has no orbit plane.

# Example
```jldoctest
orbit_normal([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9])

# output

3-element Vector{Float64}:
 -0.16225693179067105
 -0.12777026684561554
  0.9784413354903183
```
"""
function orbit_normal(r::AbstractVector, v::AbstractVector)
    _check_state_vectors(r, v)
    h = _cross3(r, v)
    magnitude = _norm3(h)
    magnitude > 0 || throw(ArgumentError(
        "the orbit normal is undefined for a radial trajectory (r × v = 0)."))
    return h ./ magnitude
end

"""
    orbital_angular_velocity(r::AbstractVector, v::AbstractVector)

Angular velocity vector of the position of a two-body state about the central body,
`ω = (r × v) / r²`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].

# Returns
The 3-component vector `ω` [rad/s], along the orbit normal.

# Example
```jldoctest
orbital_angular_velocity([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9])

# output

3-element Vector{Float64}:
 -0.00017039382545022759
 -0.00013417771620819315
  0.0010275084108450426
```
"""
function orbital_angular_velocity(r::AbstractVector, v::AbstractVector)
    _check_state_vectors(r, v)
    return _cross3(r, v) ./ _dot3(r, r)
end

"""
    eccentricity_vector(r::AbstractVector, v::AbstractVector, μ::Real)

Eccentricity vector of a two-body state, `e = ((v² − μ/r)·r − (r·v)·v) / μ`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
The 3-component vector `e`, pointing to periapsis, with magnitude the eccentricity.

# Example
```jldoctest
eccentricity_vector([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

3-element Vector{Float64}:
 -0.005542985399239766
 -0.009541538750105955
 -0.0021651913908546376
```
"""
function eccentricity_vector(r::AbstractVector, v::AbstractVector, μ::Real)
    _check_state_vectors(r, v)
    return ((_dot3(v, v) - μ / _norm3(r)) .* r .- _dot3(r, v) .* v) ./ μ
end

"""
    orbital_period(r::AbstractVector, v::AbstractVector, μ::Real)

Period of the closed two-body orbit through a state, `T = 2π·√(a³/μ)` with `a = −μ/(2ε)`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
The period [s].

# Notes
Throws `ArgumentError` unless the specific energy is negative.

# Example
```jldoctest
orbital_period([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

5909.984957722748
```
"""
function orbital_period(r::AbstractVector, v::AbstractVector, μ::Real)
    ε = specific_energy(r, v, μ)
    ε < 0 || throw(ArgumentError(
        "the orbital period is defined only for a closed orbit (negative specific energy); this " *
        "orbit has specific energy $(ε)."))
    a = -μ / (2ε)
    return 2π * sqrt(a^3 / μ)
end

"""
    mean_motion(r::AbstractVector, v::AbstractVector, μ::Real)

Mean motion of the two-body orbit through a state, `n = √(μ/|a|³)` with `a = −μ/(2ε)`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
The mean motion [rad/s]. On a hyperbola it is the rate of the hyperbolic mean anomaly.

# Notes
Throws `ArgumentError` on a parabola (zero specific energy).

# Example
```jldoctest
mean_motion([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

0.0010631474279759658
```
"""
function mean_motion(r::AbstractVector, v::AbstractVector, μ::Real)
    ε = specific_energy(r, v, μ)
    ε == 0 && throw(ArgumentError("a parabolic orbit (zero specific energy) has no mean motion."))
    return sqrt(μ / abs(μ / (2ε))^3)
end

"""
    hyperbolic_excess_speed(r::AbstractVector, v::AbstractVector, μ::Real)

Speed at infinity of the escape trajectory through a state, `v∞ = √(2ε)`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
`v∞` [km/s].

# Notes
Throws `ArgumentError` unless the specific energy is positive.

# Example
```jldoctest
hyperbolic_excess_speed([7000.0, 300.0, 0.0], [0.0, 12.5, 1.0], 398600.4418)

# output

6.5930724754364824
```
"""
function hyperbolic_excess_speed(r::AbstractVector, v::AbstractVector, μ::Real)
    ε = specific_energy(r, v, μ)
    ε > 0 || throw(ArgumentError(
        "the hyperbolic excess speed is defined only on an escape orbit (positive specific " *
        "energy); this orbit has specific energy $(ε)."))
    return sqrt(2ε)
end

"""
    apoapsis_speed(r::AbstractVector, v::AbstractVector, μ::Real)

Speed at apoapsis of the closed two-body orbit through a state, `vₐ = μ(1 − e)/h`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
`vₐ` [km/s].

# Notes
Throws `ArgumentError` unless the eccentricity is below 1.

# Example
```jldoctest
apoapsis_speed([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

7.42722369879468
```
"""
function apoapsis_speed(r::AbstractVector, v::AbstractVector, μ::Real)
    e = _norm3(eccentricity_vector(r, v, μ))
    e < 1 || throw(ArgumentError(
        "the apoapsis speed is defined only for a closed orbit (eccentricity below 1); this orbit " *
        "has eccentricity $(e)."))
    return μ * (1 - e) / _norm3(_cross3(r, v))
end

"""
    periapsis_speed(r::AbstractVector, v::AbstractVector, μ::Real)

Speed at periapsis of the two-body orbit through a state, `vₚ = μ(1 + e)/h`.

# Arguments
- `r`: Position relative to the central body, 3 components [km].
- `v`: Velocity relative to the central body, 3 components [km/s].
- `μ`: Gravitational parameter of the central body [km³/s²].

# Returns
`vₚ` [km/s], for any conic.

# Example
```jldoctest
periapsis_speed([7000.0, 300.0, 1200.0], [-0.4, 7.4, 0.9], 398600.4418)

# output

7.596164151814829
```
"""
function periapsis_speed(r::AbstractVector, v::AbstractVector, μ::Real)
    e = _norm3(eccentricity_vector(r, v, μ))
    return μ * (1 + e) / _norm3(_cross3(r, v))
end
