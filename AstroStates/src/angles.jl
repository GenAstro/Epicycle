# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Angles are recovered with atan(y, x) throughout. acos(clamp(c)) followed by a sign flip loses
# precision near 0 and π (acos resolves about 1e-8 rad there) and has no derivative at those
# points, where ForwardDiff returns NaN or a wrong finite value: at periapsis, at a zero node, at
# a zero argument of periapsis.

"""
    _wrap_2pi(x)

`x` from `atan`, in (-π, π], moved to [0, 2π). A tiny negative `x` rounds to 2π when 2π is
added, so that case returns 0; `-0.0` returns `0.0`.
"""
function _wrap_2pi(x)
    x = float(x)
    twopi = 2 * oftype(x, π)
    y = x < 0 ? x + twopi : x + zero(x)
    return y >= twopi ? y - twopi : y
end

"""
    _mod_2pi(x)

Any angle `x` reduced to [0, 2π). `mod(x, 2π)` alone returns 2π itself for a tiny negative `x`,
which rounds to 2π when 2π is added; that case returns 0.
"""
function _mod_2pi(x)
    x = float(x)    # an integer anomaly cannot hold π
    twopi = 2 * oftype(x, π)
    y = mod(x, twopi)
    return y >= twopi ? zero(y) : y
end

"""
    _plane_angle(v, p̂, q̂)

The angle of `v` in the plane spanned by the unit vectors `p̂` and `q̂`, measured from `p̂` toward
`q̂`, in [0, 2π).
"""
_plane_angle(v, p̂, q̂) = _wrap_2pi(atan(dot(q̂, v), dot(p̂, v)))

"""
    _orbit_orientation(h̄, ē, e, tol) -> (i, Ω, ω, p̂, q̂)

Inclination, right ascension of the ascending node and argument of periapsis of an orbit with
angular momentum `h̄` and eccentricity vector `ē` of magnitude `e`, and the in-plane axes `p̂`
(the direction the anomaly is measured from) and `q̂ = ĥ × p̂`.

The conventions where an angle is undefined are GMAT's:
- equatorial, `sin i ≤ tol`: Ω = 0. Prograde, ω is measured from +x toward +y; retrograde
  (i = π), from +x toward -y, the direction of motion, which is what `kep_to_cart` inverts with
  Ω = 0.
- circular, `e ≤ tol`: ω = 0, and the anomaly is measured from the ascending node (the argument of
  latitude), or from +x when the orbit is also equatorial (the true longitude).
"""
function _orbit_orientation(h̄v::AbstractVector, ēv::AbstractVector, e, tol)
    T  = promote_type(eltype(h̄v), eltype(ēv), typeof(e))
    h̄  = SVector{3,T}(h̄v[1], h̄v[2], h̄v[3])
    ē  = SVector{3,T}(ēv[1], ēv[2], ēv[3])
    h  = norm(h̄)
    ĥ  = h̄ / h
    hxy = hypot(h̄[1], h̄[2])
    # atan(0, ±h) is exact, but its derivative through hypot(0, 0) is 0/0; an exactly equatorial
    # orbit has no inclination derivative to give.
    i = iszero(hxy) ? (h̄[3] > 0 ? zero(T) : T(π)) : atan(hxy, h̄[3])

    if hxy <= tol * h                          # equatorial: no node
        Ω = zero(T)
        s = ĥ[3] > 0 ? 1 : -1                  # retrograde measures in the direction of motion
        if e <= tol
            ω  = zero(T)
            p̂  = SVector{3,T}(1, 0, 0)
        else
            ω  = _wrap_2pi(atan(s * ē[2], ē[1]))
            p̂  = ē / e
        end
    else
        n̂ = SVector{3,T}(-h̄[2], h̄[1], 0) / hxy       # ascending node, ẑ × ĥ normalised
        Ω = _wrap_2pi(atan(n̂[2], n̂[1]))
        if e <= tol
            ω  = zero(T)
            p̂  = n̂
        else
            ω  = _plane_angle(ē, n̂, cross(ĥ, n̂))
            p̂  = ē / e
        end
    end
    return i, Ω, ω, p̂, cross(ĥ, p̂)
end

"""
    _eccentric_longitude(λ, h, k; tol=1e-14, maxiter=60)

Solve the equinoctial Kepler equation `λ = F + h cos F - k sin F` for the eccentric longitude `F`,
for an ellipse, `h² + k² < 1`. Newton from `F = λ`, kept inside the bracket
[λ - (|h| + |k|), λ + (|h| + |k|)], where the root lies because `|h cos F - k sin F| ≤ |h| + |k|`;
a step that leaves it is replaced by bisection. Working in `F` directly, rather than through the
eccentricity and the longitude of periapsis, keeps the derivative defined at a circular orbit.
"""
function _eccentric_longitude(λ, h, k; tol=1e-14, maxiter=60)
    bound = abs(h) + abs(k)
    lo, hi = λ - bound, λ + bound
    F = λ + zero(h) + zero(k)
    for _ in 1:maxiter
        sinF, cosF = sincos(F)
        f = F + h * cosF - k * sinF - λ
        if f < 0
            lo = F
        elseif f > 0
            hi = F
        end
        Fnew = F - f / (1 - h * sinF - k * cosF)
        if !(lo < Fnew < hi)
            Fnew = (lo + hi) / 2
        end
        converged = abs(Fnew - F) < tol * (1 + abs(F))
        F = Fnew
        converged && break
    end
    return F
end

"""
    _has_nan(state, μ...)

Whether the state or μ holds a NaN. A NaN input comes from a conversion that has already failed
and warned, so the vector functions return NaNs for it without a second warning.
"""
_has_nan(state, μ...) = any(isnan, state) || any(isnan, μ)

"""
    _inclination_out_of_range(i)

Whether an inclination lies outside [0, π] by more than rounding. The bounds sit a few ulps
outside the range: `Float32(π)` is above π, and `cart_to_kep` returns it for a retrograde
equatorial orbit; and a ForwardDiff dual compared with an equal value breaks the tie on its
partials, so a dual of exactly π would read as above `Float64(π)`. A NaN is in range, so it passes
on without a second warning.
"""
function _inclination_out_of_range(i)
    slack = 4 * eps(float(typeof(i)))
    return i < -slack || i > π + slack
end
