# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The edge cases found in the 2026-09-30 production review, one testset each: states where a
# conversion used to return a wrong answer, a NaN without a warning, or an angle outside its
# documented range.

using Test
using AstroStates
using LinearAlgebra
using Logging

const μE = 398600.4415
_angdiff(a, b) = abs(mod(a - b + π, 2π) - π)
_quiet(f) = with_logger(f, NullLogger())

@testset "Conversion edge cases" begin

    @testset "SphericalRADECState names ra and dec as it stores them" begin
        # A point on +y, moving along +z: ra = 90°, dec = 0, rav anything, decv = 90°
        s = SphericalRADECState(CartesianState([0.0, 7000, 0, 0, 0, 7.5]), μE)
        @test s.ra ≈ π/2
        @test abs(s.dec) < 1e-15
        @test s.decv ≈ π/2
        @test fieldnames(SphericalRADECState) == (:r, :ra, :dec, :v, :rav, :decv)
        # Built field by field, it lands where its names say
        c = CartesianState(SphericalRADECState(7000.0, π/4, 0.0, 7.5, 3π/4, 0.0))
        @test c.position ≈ [7000cos(π/4), 7000sin(π/4), 0]
        @test c.velocity ≈ [7.5cos(3π/4), 7.5sin(3π/4), 0]
        @test to_vector(s) == [s.r, s.ra, s.dec, s.v, s.rav, s.decv]
    end

    @testset "retrograde equatorial orbits (i = π) convert, Ω pinned to 0" begin
        for cart in ([7000.0, 0, 0, 0, -7.5, 0],                    # elliptic
                     [7000.0, 0, 0, 0, -sqrt(μE / 7000), 0],         # circular
                     [5000.0, 5000, 0, 3, -5, 0])                   # periapsis off the x-axis
            k = cart_to_kep(cart, μE)
            @test all(isfinite, k)
            @test k[3] ≈ π
            @test k[4] == 0
            @test all(0 .<= k[4:6] .< 2π)
            @test isapprox(kep_to_cart(k, μE), cart; rtol = 1e-12)
        end
    end

    @testset "flight path angle is measured from the radial direction" begin
        # Horizontal circular motion is fpa = π/2; outward radial motion is 0
        @test cart_to_sphazfpa([7000.0, 0, 0, 0, 7.5, 0])[6] ≈ π/2
        @test cart_to_sphazfpa([7000.0, 0, 0, 7.5, 0, 0])[6] == 0
        @test cart_to_sphazfpa([7000.0, 0, 0, -7.5, 0, 0])[6] ≈ π
        # fpa = π/2 and azimuth π/2 (east) is horizontal eastward flight
        @test sphazfpa_to_cart([7000.0, 0, 0, 7.5, π/2, π/2]) ≈ [7000.0, 0, 0, 0, 7.5, 0] atol = 1e-12
    end

    @testset "angles have full precision near 0 and π" begin
        @test cart_to_kep(kep_to_cart([8000.0, 0.1, 0.5, 0.3, 0.4, 1e-9], μE), μE)[6] ≈ 1e-9 rtol = 1e-6
        @test cart_to_kep(kep_to_cart([8000.0, 0.1, 1e-7, 0.3, 0.4, 1.0], μE), μE)[3] ≈ 1e-7 rtol = 1e-8
    end

    @testset "MEE retrograde set (j = -1) round-trips" begin
        for k in ([8000.0, 0.2, 3.0, 1, 1, 1], [8000.0, 0.2, 0.3, 1, 1, 1])
            c = kep_to_cart(k, μE)
            for j in (1.0, -1.0)
                m = cart_to_mee(c, μE; j = j)
                @test 0 <= m[6] < 2π
                @test isapprox(mee_to_cart(m, μE; j = j), c; rtol = 1e-12)
            end
        end
    end

    @testset "equinoctial elements refuse orbits near i = π instead of losing precision" begin
        for di in (1e-5, 1e-6, 0.0)
            c = kep_to_cart([8000.0, 0.1, π - di, 0.3, 0.4, 1.0], μE)
            @test_logs (:warn, r"i ≈ π") @test all(isnan, cart_to_equinoctial(c, μE))
        end
        # Away from the singularity the round trip is exact to rounding
        c = kep_to_cart([8000.0, 0.1, π - 1e-3, 0.3, 0.4, 1.0], μE)
        @test isapprox(equinoctial_to_cart(cart_to_equinoctial(c, μE), μE), c; rtol = 1e-9)
    end

    @testset "a hyperbolic true anomaly beyond the asymptote is refused" begin
        @test_logs (:warn, r"beyond the asymptote") @test all(isnan, kep_to_cart([-10000, 2.0, 0.5, 0, 0, 2.5], μE))
        @test_logs (:warn, r"beyond the asymptote") @test all(isnan, mee_to_cart([20000.0, 2.0, 0, 0, 0, π], μE))
        # Just inside the asymptote converts
        @test all(isfinite, kep_to_cart([-10000, 2.0, 0.5, 0, 0, acos(-1/2) - 1e-3], μE))
    end

    @testset "modified Keplerian: a hyperbola needs |rₐ| > rₚ" begin
        for ra in (-5000.0, -7000.0)
            @test_logs (:warn, r"\|rₐ\| must exceed rₚ") @test all(isnan, modkep_to_kep([7000.0, ra, 0.1, 0, 0, 0]))
        end
        k = modkep_to_kep([7000.0, -9000.0, 0.1, 0, 0, 0])
        @test k[2] > 1 && k[1] < 0
    end

    @testset "asymptote elements validate the periapsis radius" begin
        @test_logs (:warn, r"must be positive") @test all(isnan, outasymptote_to_kep([-7000.0, 10, 0, 0.3, 0.5, 2.9], μE))
        @test_logs (:warn, r"exceeds the semi-major axis") @test all(isnan, inasymptote_to_kep([50000.0, -10, 0, 0.3, 0.5, 1.0], μE))
    end

    @testset "retrograde equatorial asymptote gives ω in [0, 2π)" begin
        for θ in (3π/2, 3π/2 + 1e-14, 3π/2 - 1e-14)
            k = _quiet(() -> outasymptote_to_kep([7000.0, 10, 0.4, 0.0, θ, 0.5], μE))
            @test all(0 .<= k[4:6] .< 2π)
        end
    end

    @testset "MEE degenerate states warn and return NaN" begin
        @test_logs (:warn, r"degenerate position or velocity") @test all(isnan, cart_to_mee([0.0, 0, 0, 1, 0, 0], μE))
        @test_logs (:warn, r"degenerate position or velocity") @test all(isnan, cart_to_mee([7000.0, 0, 0, 0, 0, 0], μE))
        @test_logs (:warn, r"degenerate angular momentum") @test all(isnan, cart_to_mee([7000.0, 0, 0, 1, 0, 0], μE))
        @test_logs (:warn, r"μ < tolerance") @test all(isnan, cart_to_mee([7000.0, 0, 0, 0, 7.5, 0], 0.0))
        @test_logs (:warn, r"Semi-latus rectum") @test all(isnan, mee_to_cart([0.0, 0, 0, 0, 0, 0], μE))
    end

    @testset "equinoctial Kepler solve converges up to e = 0.999" begin
        for e in (0.99, 0.999), M in range(0, 2π; length = 181)
            c = equinoctial_to_cart([20000.0, e * sin(1.0), e * cos(1.0), 0.1, 0.2, M], μE)
            @test all(isfinite, c)
        end
    end

    @testset "the Kepler solve keeps a root it lands on exactly" begin
        # GMAT's elliptic truth state: Newton reaches the root with f = 0 on its fourth step. The
        # solver used to discard it for a bisection midpoint and stop 7e-14 rad away, 0.3 µm here.
        h, k, λ = 0.1879385241571815, 0.0684040286651338, deg2rad(28.36066564454829)
        e, φ = hypot(h, k), atan(h, k)
        E = AstroStates._eccentric_anomaly(λ - φ, e)
        @test abs(E - e * sin(E) - mod(λ - φ, 2π)) < 4eps()
        c = equinoctial_to_cart([8000.0, h, k, -0.08626412365266437, 0.07238419434078193, λ], μE)
        gmat = [6759.747343616322723, 1115.043329211011041, 1344.722777534846955,
                -2.660243619064134, 7.541202154282467, 0.640887592324028]
        @test maximum(abs.(c[1:3] - gmat[1:3])) < 1e-11        # km
        @test maximum(abs.(c[4:6] - gmat[4:6])) < 1e-14        # km/s
    end

    @testset "equinoctial elements are elliptic only" begin
        @test_logs (:warn, r"elliptic orbits only") @test all(isnan, equinoctial_to_cart([-7000.0, 0.1, 0.1, 0.1, 0.1, 1], μE))
        @test_logs (:warn, r"parabolic or hyperbolic") @test all(isnan, cart_to_equinoctial(kep_to_cart([-7000.0, 1.5, 0.3, 0, 0, 0], μE), μE))
    end

    @testset "angles come back in their documented ranges" begin
        c = [7000.0, -100, -300, -1, -7.5, -0.2]
        r = cart_to_sphradec(c)
        @test 0 <= r[2] < 2π && 0 <= r[5] < 2π
        a = cart_to_sphazfpa(c)
        @test 0 <= a[2] < 2π && 0 <= a[5] < 2π && 0 <= a[6] <= π
        @test 0 <= cart_to_equinoctial(kep_to_cart([8000.0, 0.1, 0.5, 6.2, 6.2, 6.2], μE), μE)[6] < 2π
        for k in (cart_to_kep(c, μE), outasymptote_to_kep(cart_to_outasymptote(c, μE), μE))
            @test all(0 .<= k[4:6] .< 2π)
        end
    end

    @testset "element types are kept: Float32 in, Float32 out" begin
        c32 = Float32[7000, 100, 300, 0.1, 7.5, 0.2]
        for f in (x -> cart_to_kep(x, 398600f0), cart_to_sphradec, cart_to_sphazfpa,
                  x -> cart_to_mee(x, 398600f0), x -> cart_to_equinoctial(x, 398600f0),
                  x -> cart_to_outasymptote(x, 398600f0))
            @test eltype(f(c32)) == Float32
        end
        @test eltype(kep_to_cart(Float32[7000, 0.01, 0.1, 0, 0, 0], 398600f0)) == Float32
        @test KeplerianState(CartesianState(c32), 398600f0) isa KeplerianState{Float32}
    end

    @testset "any AbstractVector is accepted" begin
        v = [7000.0, 0, 100, 0, 7.5, 0.1, 99.0]
        @test cart_to_kep(view(v, 1:6), μE) == cart_to_kep(v[1:6], μE)
        @test cart_to_sphradec(view(v, 1:6)) == cart_to_sphradec(v[1:6])
    end

    @testset "constructors promote mixed element types" begin
        @test KeplerianState(7000.0, 0.01, 0, 0, 0, 0) isa KeplerianState{Float64}
        @test CartesianState([7000, 0, 0], [0.0, 7.5, 0.0]) isa CartesianState{Float64}
        @test KeplerianState(KeplerianState(7000.0, 0.01, 0.1, 0, 0, 0)) isa KeplerianState
    end

    @testset "CartesianState compares by value" begin
        a = CartesianState([7000.0, 0, 0, 0, 7.5, 0])
        @test a == CartesianState([7000.0, 0, 0, 0, 7.5, 0])
        @test hash(a) == hash(CartesianState([7000.0, 0, 0, 0, 7.5, 0]))
        @test a ≈ CartesianState([7000.0 + 1e-9, 0, 0, 0, 7.5, 0])
        @test !(a ≈ KeplerianState(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    end

    @testset "malformed inputs throw ArgumentError" begin
        @test_throws ArgumentError OrbitState([1.0, 2, 3], Keplerian())
        @test_throws ArgumentError KeplerianState([1.0, 2, 3], μE)
        @test_throws ArgumentError CartesianState([1.0, 2, 3])
    end

    @testset "a straight-line trajectory (r ∥ v) warns and returns NaN" begin
        line = CartesianState([7000.0, 7000, 7000, 7, 7, 7])
        @test_logs (:warn, r"degenerate angular momentum") @test all(isnan, to_vector(KeplerianState(line, μE)))
        @test_logs (:warn, r"zero angular momentum") @test all(isnan, to_vector(OutGoingAsymptoteState(line, μE)))
        # Spherical coordinates need no orbit, so they round-trip
        @test CartesianState(SphericalRADECState(line)) ≈ line
    end
end
