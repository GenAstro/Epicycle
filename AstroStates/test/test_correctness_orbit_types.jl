# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Circular and parabolic orbits in the representations that define them [AS-5.2, AS-5.4],
# against closed forms. Elliptic and hyperbolic orbits are covered by the truth tables in
# runtests.jl, equatorial ones by runtests_kep_cart.jl.

using Test
using LinearAlgebra
using AstroStates

let μ = 398600.4415
    rel(a, b) = norm(a - b) / norm(b)

    @testset "Circular orbits [AS-5.2]" begin
        r, i, Ω, u = 7000.0, 0.5, 0.3, 1.0          # u is the argument of latitude
        cart = kep_to_cart([r, 0.0, i, Ω, 0.0, u], μ)

        mee = cart_to_mee(cart, μ)
        @test mee[1] ≈ r rtol = 1e-12                               # p = r
        @test hypot(mee[2], mee[3]) < 1e-12                         # f, g: no eccentricity
        @test hypot(mee[4], mee[5]) ≈ tan(i / 2) rtol = 1e-12       # h, k: the node
        @test mee[6] ≈ Ω + u rtol = 1e-12                           # true longitude
        @test rel(mee_to_cart(mee, μ), cart) < 1e-12

        eq = cart_to_equinoctial(cart, μ)
        @test eq[1] ≈ r rtol = 1e-12                                # a = r
        @test hypot(eq[2], eq[3]) < 1e-12
        @test hypot(eq[4], eq[5]) ≈ tan(i / 2) rtol = 1e-12
        @test eq[6] ≈ Ω + u rtol = 1e-12                            # mean longitude
        @test rel(equinoctial_to_cart(eq, μ), cart) < 1e-12

        modkep = kep_to_modkep(cart_to_kep(cart, μ))
        @test modkep[1] ≈ r rtol = 1e-12                            # rₚ = r
        @test modkep[2] ≈ r rtol = 1e-12                            # rₐ = r

        azfpa = cart_to_sphazfpa(cart)
        @test azfpa[1] ≈ r rtol = 1e-12
        @test azfpa[4] ≈ sqrt(μ / r) rtol = 1e-12                   # circular speed
        @test azfpa[6] ≈ π / 2 atol = 1e-12                         # horizontal flight
        @test rel(sphazfpa_to_cart(azfpa), cart) < 1e-12
        @test rel(sphradec_to_cart(cart_to_sphradec(cart)), cart) < 1e-12
    end

    @testset "Parabolic orbits [AS-5.4]" begin
        rp, i, Ω, ω = 7000.0, 0.5, 0.3, 0.4
        h̄ = sqrt(2μ * rp)                          # angular momentum of a parabola, p = 2rₚ
        for ν in (0.0, 1.0, -1.5)
            r = 2rp / (1 + cos(ν))
            # perifocal position and velocity, rotated by Ω, i, ω
            R = [cos(Ω) -sin(Ω) 0; sin(Ω) cos(Ω) 0; 0 0 1] *
                [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)] *
                [cos(ω) -sin(ω) 0; sin(ω) cos(ω) 0; 0 0 1]
            pos = R * [r * cos(ν), r * sin(ν), 0.0]
            vel = R * ((μ / h̄) .* [-sin(ν), 1 + cos(ν), 0.0])
            cart = vcat(pos, vel)

            mee = cart_to_mee(cart, μ)
            @test mee[1] ≈ 2rp rtol = 1e-12                         # p = 2rₚ
            @test hypot(mee[2], mee[3]) ≈ 1.0 rtol = 1e-12          # e = 1
            @test rel(mee_to_cart(mee, μ), cart) < 1e-12

            kep = cart_to_kep(cart, μ)
            @test isinf(kep[1])
            @test kep[2] ≈ 1.0 rtol = 1e-12
            @test rel(sphazfpa_to_cart(cart_to_sphazfpa(cart)), cart) < 1e-12
            @test rel(sphradec_to_cart(cart_to_sphradec(cart)), cart) < 1e-12
        end
    end
end
