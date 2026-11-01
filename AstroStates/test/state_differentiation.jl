# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Differentiability of the conversions. For each representation, the ForwardDiff Jacobian of
# Cartesian → representation times the Jacobian of representation → Cartesian must be the
# identity. That catches a derivative that is NaN, infinite, or finite and wrong, which is what
# acos-based angle recovery gave at periapsis, apoapsis, and a zero node or argument of periapsis.

using Test
using ForwardDiff
using LinearAlgebra
using AstroStates

const μD = 398600.4415

# Each representation as a pair of functions on 6-vectors: from Cartesian, to Cartesian.
const _AD_PAIRS = (
    Keplerian           = (c -> cart_to_kep(c, μD),               k -> kep_to_cart(k, μD)),
    ModifiedKeplerian   = (c -> kep_to_modkep(cart_to_kep(c, μD)), m -> kep_to_cart(modkep_to_kep(m), μD)),
    SphericalRADEC      = (cart_to_sphradec,                      sphradec_to_cart),
    SphericalAZIFPA     = (cart_to_sphazfpa,                      sphazfpa_to_cart),
    ModifiedEquinoctial = (c -> cart_to_mee(c, μD),               m -> mee_to_cart(m, μD)),
    Equinoctial         = (c -> cart_to_equinoctial(c, μD),       q -> equinoctial_to_cart(q, μD)),
    AltEquinoctial      = (c -> equinoctial_to_alt_equinoctial(cart_to_equinoctial(c, μD)),
                           q -> equinoctial_to_cart(alt_equinoctial_to_equinoctial(q), μD)),
    OutGoingAsymptote   = (c -> cart_to_outasymptote(c, μD),      a -> kep_to_cart(outasymptote_to_kep(a, μD), μD)),
    IncomingAsymptote   = (c -> cart_to_inasymptote(c, μD),       a -> kep_to_cart(inasymptote_to_kep(a, μD), μD)),
)

# Keplerian points, among them the angle values acos could not differentiate.
const _AD_POINTS = (
    general       = [8000.0, 0.1, 0.5, 0.3, 0.4, 1.0],
    periapsis     = [8000.0, 0.1, 0.5, 0.3, 0.4, 0.0],
    apoapsis      = [8000.0, 0.1, 0.5, 0.3, 0.4, π],
    zero_aop      = [8000.0, 0.1, 0.5, 0.3, 0.0, 1.0],
    zero_raan     = [8000.0, 0.1, 0.5, 0.0, 0.4, 1.0],
    retrograde    = [8000.0, 0.1, 2.8, 0.3, 0.4, 1.0],
    hyperbolic    = [-20000.0, 1.5, 0.5, 0.3, 0.4, 0.5],
)

@testset "Conversions are differentiable: J(from Cartesian) * J(to Cartesian) = I" begin
    for (pname, k) in pairs(_AD_POINTS), (rname, (from, to)) in pairs(_AD_PAIRS)
        k[1] < 0 && rname in (:Equinoctial, :AltEquinoctial) && continue   # elliptic only
        c = kep_to_cart(k, μD)
        x = from(c)
        J1 = ForwardDiff.jacobian(from, c)
        J2 = ForwardDiff.jacobian(to, x)
        @testset "$rname at $pname" begin
            @test all(isfinite, J1)
            @test norm(J1 * J2 - I) < 1e-8
        end
    end
end

@testset "Derivatives pass through the state structs" begin
    c = kep_to_cart(_AD_POINTS.periapsis, μD)
    J = ForwardDiff.jacobian(x -> to_vector(KeplerianState(CartesianState(x), μD)), c)
    @test J ≈ ForwardDiff.jacobian(x -> cart_to_kep(x, μD), c)
    g = ForwardDiff.derivative(m -> to_vector(KeplerianState(CartesianState(c), m))[1], μD)
    @test isfinite(g)
end

# The identity above shows the two Jacobians are inverses; this anchors one of them to an
# independent truth. Each to-Cartesian Jacobian is checked against central finite differences.
# The to-Cartesian direction is the one to difference: its output is continuous, where an angle
# coming out of the from-Cartesian direction wraps at 2π.
function _central_difference(f, x)
    J = zeros(6, 6)
    for k in 1:6
        h = 1e-6 * max(1.0, abs(x[k]))
        xp = copy(x); xp[k] += h
        xm = copy(x); xm[k] -= h
        J[:, k] = (f(xp) .- f(xm)) ./ (2h)
    end
    return J
end

# Each row against its own size: position partials are thousands of times the velocity ones, so a
# whole-matrix norm would hide an error in a velocity row. A row that is zero, such as z on an
# equatorial orbit, is held to a small share of the whole.
function _rows_agree(J, Jd; rtol = 1e-6)
    floor = 1e-10 * norm(J)
    return all(norm(J[i, :] - Jd[i, :]) <= rtol * max(norm(J[i, :]), floor) for i in axes(J, 1))
end

@testset "Conversions to Cartesian match finite differences [AS-9]" begin
    for (pname, k) in pairs(_AD_POINTS), (rname, (from, to)) in pairs(_AD_PAIRS)
        k[1] < 0 && rname in (:Equinoctial, :AltEquinoctial) && continue   # elliptic only
        x = from(kep_to_cart(k, μD))
        J  = ForwardDiff.jacobian(to, x)
        Jd = _central_difference(to, x)
        @testset "$rname at $pname" begin
            @test _rows_agree(J, Jd)
        end
    end
end

# The equinoctial sets define every element at a circular or equatorial orbit, so their
# derivatives exist there; a circular equatorial orbit, such as GEO, is where an optimizer often
# starts.
@testset "Equinoctial derivatives exist at circular and equatorial orbits [AS-9]" begin
    to_cart = (
        Equinoctial         = q -> equinoctial_to_cart(q, μD),
        AltEquinoctial      = q -> equinoctial_to_cart(alt_equinoctial_to_equinoctial(q), μD),
        ModifiedEquinoctial = m -> mee_to_cart(m, μD),
    )
    points = (
        geo             = [42164.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        circular        = [42164.0, 0.0, 0.0, 0.1, 0.2, 1.0],
        equatorial      = [42164.0, 0.01, 0.02, 0.0, 0.0, 1.0],
    )
    for (pname, x) in pairs(points), (rname, f) in pairs(to_cart)
        # the modified equinoctial set leads with p, not a; at these points p = a(1 - h² - k²)
        y = rname === :ModifiedEquinoctial ? [x[1] * (1 - x[2]^2 - x[3]^2); x[2:end]] : x
        @testset "$rname at $pname" begin
            J = ForwardDiff.jacobian(f, y)
            @test all(isfinite, J)
            @test _rows_agree(J, _central_difference(f, y))
        end
    end
    for f in (equinoctial_to_alt_equinoctial, alt_equinoctial_to_equinoctial)
        @test all(isfinite, ForwardDiff.jacobian(f, [42164.0, 0.0, 0.0, 0.0, 0.0, 1.0]))
    end
end

# A dual inclination of exactly 0 or π is in range. ForwardDiff breaks ties between equal values
# on the partials, so a comparison with Float64(π) would have refused it and returned a zero
# Jacobian. Finite differences would step outside [0, π], so the check is continuity: the Jacobian
# at the boundary matches the one just inside it.
@testset "Derivatives at i = 0 and i = π [AS-9]" begin
    for (i_edge, i_inside) in ((0.0, 1e-9), (Float64(π), π - 1e-9))
        k_edge   = [8000.0, 0.1, i_edge, 0.3, 0.4, 1.0]
        k_inside = [8000.0, 0.1, i_inside, 0.3, 0.4, 1.0]
        # Brouwer's osculating-to-mean iteration converges to 1e-8, so its Jacobians agree to 1e-5
        for (f, rtol) in ((k -> kep_to_cart(k, μD), 1e-8), (kep_to_modkep, 1e-8),
                          (k -> kep_to_brouwer_mean_long(k, μD), 1e-4))
            J = @test_logs match_mode=:any ForwardDiff.jacobian(f, k_edge)
            @test all(isfinite, J)
            @test norm(J) > 0
            @test norm(J - ForwardDiff.jacobian(f, k_inside)) <= rtol * norm(J)
        end
    end
end
