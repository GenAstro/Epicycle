# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

@testset "Parabolic Cartesian state cart_to_kep" begin
    # Standard Earth μ (ASCII name)
    mu = 398600.4415
    r  = 8000.0
    vmag = sqrt(2 * mu / r)          # Parabolic speed: escape speed at radius r

    # IMPORTANT: [a b c ...] creates a 1×N Matrix in Julia; use a Vector literal instead.
    cart = [r, 0.0, 0.0, 0.0, vmag, 0.0]

    @test length(cart) == 6

    kep = cart_to_kep(cart, mu)
    @test length(kep) == 6

    a, e, inc, Ω, ω, ν = kep

    # Eccentricity should be (very nearly) 1 for parabolic trajectory
    @test isapprox(e, 1.0; atol=1e-15, rtol=0)

    # Specific orbital energy ~ 0 ⇒ a -> ∞ (may be Inf or very large magnitude)
    @test (isinf(a))

    # The state types carry the same elements
    kstate = KeplerianState(CartesianState(cart), mu)
    @test kstate.ecc ≈ 1.0 atol=1e-15
    @test isinf(kstate.sma)
end
nothing
