# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Delaunay elements (#11, section 5), and the alternate equinoctial set at e = 0 and i = 0.
#
# Truth, strongest first:
# - External reference: GMAT R2026a Delaunay elements for four orbits, from
#   AstroCallbacks/test/truth/orbit_parameters_truth.script.
# - Analytic: each element from its definition in Keplerian terms, round trips through Cartesian,
#   and the alternate equinoctial set staying finite and continuous at e = 0 and i = 0.
# - Differentiability: ForwardDiff Jacobians of Cartesian → elements → Cartesian are the identity.
# =============================================================================

using Test
using ForwardDiff
using LinearAlgebra
using Logging
using AstroStates

const _ES_MU = 398600.4418

const _ES_STATES = (
    leo        = [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9],
    molniya    = [4637.031328726552, 178.53697947901998, -5679.055240387162,
                  6.252424682730313, 6.928411997008257, 2.5730558589825416],
    polar      = [1731.6965430023222, -2999.386395771405, -7427.270441381914,
                  -3.0611710106575205, 5.302103721115798, -2.5783940176845728],
    retrograde = [-5967.240072871029, 10335.562987173626, 6.230821664531257e-13,
                  3.6791756283823123, 3.4619800616646623, -2.8389753744567603])

# GMAT Delaunay elements, angles in degrees: l, g, h, L, G, H.
const _ES_DELAUNAY_TRUTH = (
    leo        = [122.5442363250099, 291.2029829261602, 308.2188917629123,
                  53067.34411774273, 53063.98873058828, 51920.0],
    molniya    = [3.137014995390693, 270.0, 40.00000000000001,
                  102969.7613471061, 69258.16876405635, 31010.97447117744],
    polar      = [204.2134756141651, 45.00000000000003, 120.0,
                  54676.35058688537, 54402.28194078811, 3.331179022254901e-12],
    retrograde = [76.9807494122247, 80.00000000000003, 300.0,
                  69160.72080017674, 67763.39048140965, -58684.8176034654])

_es_angle_diff(a, b) = abs(rem(a - b, 2π, RoundNearest))
_true_to_mean(ν, e) = (E = atan(sqrt(1 - e^2) * sin(ν), e + cos(ν)); E - e * sin(E))

@testset "Delaunay elements agree with an external reference (GMAT)" begin
    for name in keys(_ES_STATES)
        @testset "$name" begin
            d = DelaunayState(CartesianState(_ES_STATES[name]), _ES_MU)
            t = _ES_DELAUNAY_TRUTH[name]
            for (k, x) in enumerate((d.l, d.g, d.h))
                @test _es_angle_diff(x, deg2rad(t[k])) < 1e-12
            end
            @test d.L ≈ t[4] rtol = 1e-13
            @test d.G ≈ t[5] rtol = 1e-13
            @test d.H ≈ t[6] rtol = 1e-12 atol = 1e-8
        end
    end
end

@testset "elements follow their definitions" begin
    for name in keys(_ES_STATES)
        @testset "$name" begin
            cart = CartesianState(_ES_STATES[name])
            a, e, i, Ω, ω, ν = to_vector(KeplerianState(cart, _ES_MU))
            E = atan(sqrt(1 - e^2) * sin(ν), e + cos(ν))
            M = mod(E - e * sin(E), 2π)

            alt = AlternateEquinoctialState(cart, _ES_MU)
            @test alt.a ≈ a
            @test alt.h ≈ e * sin(Ω + ω) atol = 1e-14
            @test alt.k ≈ e * cos(Ω + ω) atol = 1e-14
            @test alt.altp ≈ sin(i / 2) * sin(Ω) atol = 1e-14
            @test alt.altq ≈ sin(i / 2) * cos(Ω) atol = 1e-14
            @test _es_angle_diff(alt.mlong, Ω + ω + M) < 1e-12

            d = DelaunayState(cart, _ES_MU)
            @test d.L ≈ sqrt(_ES_MU * a)
            @test d.G ≈ norm(cross(_ES_STATES[name][1:3], _ES_STATES[name][4:6])) rtol = 1e-13
            @test d.H ≈ d.G * cos(i) atol = 1e-8
        end
    end
end

@testset "round trips through Cartesian and between representations" begin
    for name in keys(_ES_STATES)
        cart = CartesianState(_ES_STATES[name])
        for S in (AlternateEquinoctialState, DelaunayState)
            back = CartesianState(S(cart, _ES_MU), _ES_MU)
            @test to_vector(back) ≈ _ES_STATES[name] rtol = 1e-12 atol = 1e-9
        end
        # Delaunay through the other representations and back, by way of KeplerianState
        d = DelaunayState(cart, _ES_MU)
        for S in (EquinoctialState, AlternateEquinoctialState, ModifiedEquinoctialState, SphericalRADECState)
            again = DelaunayState(S(d, _ES_MU), _ES_MU)
            @test to_vector(again)[4:6] ≈ to_vector(d)[4:6] rtol = 1e-11 atol = 1e-8
            for k in 1:3
                @test _es_angle_diff(to_vector(again)[k], to_vector(d)[k]) < 1e-10
            end
        end
    end
end

@testset "the alternate equinoctial set is finite and continuous at e = 0 and i = 0" begin
    # Circular and equatorial: Ω, ω and ν are each undefined, the mean longitude is not.
    for (e, i) in ((0.0, 0.3), (0.01, 0.0), (0.0, 0.0))
        λ = 1.2
        kep = KeplerianState(7000.0, e, i, 0.0, 0.0, λ)
        alt = AlternateEquinoctialState(kep, _ES_MU)
        @test all(isfinite, to_vector(alt))
        @test _es_angle_diff(alt.mlong, _true_to_mean(λ, e)) < 1e-12
        back = CartesianState(alt, _ES_MU)
        @test to_vector(back) ≈ to_vector(CartesianState(kep, _ES_MU)) rtol = 1e-12 atol = 1e-9
    end
    # A small change in eccentricity gives a small change in the elements, through e = 0.
    near = [to_vector(AlternateEquinoctialState(KeplerianState(7000.0, e, 0.3, 0.4, 0.5, 0.6), _ES_MU))
            for e in (1e-8, 1e-12)]
    @test maximum(abs, near[1] .- near[2]) < 1e-7
end

@testset "a hyperbola warns and converts to NaN" begin
    hyp = KeplerianState(-20000.0, 1.4, 0.5, 1.0, 0.3, 0.4)
    for S in (AlternateEquinoctialState, DelaunayState)
        s = with_logger(NullLogger()) do
            S(hyp, _ES_MU)
        end
        @test all(isnan, to_vector(s))
        @test_logs (:warn, r"hyperbolic|closed orbit") S(hyp, _ES_MU)
    end
end

@testset "OrbitState tags, vector constructors and display" begin
    for (S, tag) in ((AlternateEquinoctialState, AlternateEquinoctial()), (DelaunayState, Delaunay()))
        @test state_tag_to_type(tag) === S
        @test state_type_to_tag(S) === tag
        s = S(CartesianState(_ES_STATES[:leo]), _ES_MU)
        @test S(to_vector(s)) == s
        @test OrbitState(s).statetype === tag
        @test occursin(string(nameof(S)), sprint(show, s))
        @test_throws ArgumentError S([1.0, 2.0, 3.0])
    end
    @test DelaunayState(1.0, 2.0, 3, 53000.0, 52000.0, 51000) isa DelaunayState{Float64}
end

_es_to_alt(c)   = equinoctial_to_alt_equinoctial(cart_to_equinoctial(c, _ES_MU))
_es_from_alt(n) = equinoctial_to_cart(alt_equinoctial_to_equinoctial(n), _ES_MU)

@testset "differentiable: Cartesian → elements → Cartesian has an identity Jacobian" begin
    to_d(c)    = kep_to_delaunay(cart_to_kep(c, _ES_MU), _ES_MU)
    from_d(d)  = kep_to_cart(delaunay_to_kep(d, _ES_MU), _ES_MU)
    for name in keys(_ES_STATES), (f, g) in ((_es_to_alt, _es_from_alt), (to_d, from_d))
        c = _ES_STATES[name]
        J = ForwardDiff.jacobian(g, f(c)) * ForwardDiff.jacobian(f, c)
        @test J ≈ I(6) atol = 1e-8
    end
end

@testset "the alternate equinoctial set is differentiable at e = 0 and i = 0" begin
    # The point of the set: a circular or equatorial orbit has a finite Jacobian, which a solver
    # varying these elements needs. Truth: the identity, and central differences.
    for (e, i) in ((0.0, 0.3), (0.01, 0.0), (0.0, 0.0))
        c = kep_to_cart([7000.0, e, i, 0.4, 0.5, 0.6], _ES_MU)
        n = _es_to_alt(c)
        Jfrom = ForwardDiff.jacobian(_es_from_alt, n)
        @test all(isfinite, Jfrom)
        @test Jfrom * ForwardDiff.jacobian(_es_to_alt, c) ≈ I(6) atol = 1e-8
        δ = 1e-6 .* max.(abs.(n), 1e-3)
        Jcd = reduce(hcat, [(_es_from_alt(n .+ δ[k] .* (1:6 .== k)) .-
                             _es_from_alt(n .- δ[k] .* (1:6 .== k))) ./ (2δ[k]) for k in 1:6])
        @test Jfrom ≈ Jcd rtol = 1e-6
    end
end

@testset "circular and equatorial orbits read back with the Keplerian conventions" begin
    # Truth: the conversion from CartesianState, which sets ω = 0 on a circular orbit and Ω = 0
    # on an equatorial one.
    for (e, i) in ((0.0, 0.3), (0.01, 0.0), (0.0, 0.0))
        cart = CartesianState(kep_to_cart([7000.0, e, i, 0.4, 0.5, 0.6], _ES_MU))
        from_alt  = to_vector(KeplerianState(AlternateEquinoctialState(cart, _ES_MU), _ES_MU))
        from_cart = to_vector(KeplerianState(cart, _ES_MU))
        @test from_alt[1:3] ≈ from_cart[1:3] rtol = 1e-12 atol = 1e-12
        for k in 4:6
            @test _es_angle_diff(from_alt[k], from_cart[k]) < 1e-10
        end
    end
end

@testset "the alternate equinoctial set is singular at i = π" begin
    singular = [7000.0, 0.01, 0.0, 0.0, 1.0, 0.6]          # sin(i/2) = 1
    @test all(isnan, with_logger(NullLogger()) do
        alt_equinoctial_to_equinoctial(singular)
    end)
    @test_logs (:warn, r"180°") alt_equinoctial_to_equinoctial(singular)
end

nothing
