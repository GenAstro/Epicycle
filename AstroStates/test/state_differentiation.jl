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
