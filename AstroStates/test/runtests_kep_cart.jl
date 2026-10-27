# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

@testset "cart_to_kep Input Validation" begin
    # Invalid length
    err = try
        cart_to_kep([1.0,2.0,3.0,4.0,5.0,6.0,7.0], μ)
        nothing
    catch e
        e
    end
    @test isa(err, ErrorException)
    @test occursin("Input vector must have exactly six elements", err.msg)

    # Degenerate position vector
    @testset "Degenerate r" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([0.0,0.0,0.0, 1.0,1.0,1.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Orbit is singular due to degenerate position", log_output)
    end

    # Degenerate velocity vector
    @testset "Degenerate v" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([1.0,1.0,1.0, 0.0,0.0,0.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion failed: Orbit is singular due to degenerate position or velocity vector", log_output)
    end

    # r ∥ v (zero angular momentum)
    @testset "Parallel r and v" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([1.0,1.0,1.0, 2.0,2.0,2.0], μ)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion failed: Orbit is singular due to degenerate angular momentum", log_output)
    end

    # Invalid μ
    @testset "Zero mu" begin
        io = IOBuffer()
        Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
            result = cart_to_kep([7000.0,0.0,0.0, 0.0,7.5,0.0], 0.0)
            @test all(isnan.(result))
        end
        log_output = String(take!(io))
        @test occursin("Conversion failed: μ < tolerance", log_output)
    end
end

@testset "kep_to_cart Input Validation" begin
    # Too many elements
    err = try
        kep_to_cart([7000.0, 0.0, 0.0, 0.0, pi, pi/2, 3.0], μ)
        nothing
    catch e
        e
    end
    @test isa(err, ErrorException)
    @test occursin("Input vector must have exactly six elements: a, e, i, Ω, ω, ν", err.msg)

    # Parabolic orbit
    io = IOBuffer()
    Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
        state = kep_to_cart([7000.0, 1.0, 0.0, 0.0, pi, pi/2], μ)
        @test all(isnan.(state))
    end
    output = String(take!(io))
    @test occursin("Conversion failed: Orbit is parabolic or singular.", output)

    # μ = 0 case
    io = IOBuffer()
    Logging.with_logger(ConsoleLogger(io, Logging.Warn)) do
        state = kep_to_cart([7000.0, 0.1, 0.0, 0.0, pi, pi/2], 0.0)
        @test all(isnan.(state))
    end
    output = String(take!(io))
    @test occursin("Conversion failed: μ < tolerance.", output)
end

@testset "cart_to_kep Quadrant Tests [AS-4, AS-6]" begin
    sma = 42000.0
    ecc = 0.05
    i_tests = [pi/4, pi/2, 3pi/4, pi]
    angle_tests = [pi/4.1, 3.1*pi/4.0, 5.0*pi/4.1, 7.0*pi/4.0]

    for i in i_tests
        kep = [sma, ecc, i, 0.1, 1.1, 0.5]
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        # At i = π the node is undefined: Ω is pinned to 0, and ω becomes the longitude of
        # periapsis measured in the direction of motion, ω - Ω.
        expected = isapprox(i, pi) ? [sma, ecc, i, 0.0, 1.1 - 0.1, 0.5] : kep
        @test isapproxvec_percent(expected, kep2; tol = 1e-13)
        @test isapproxvec_percent(cart, kep_to_cart(kep2, μ); tol = 1e-13)
    end

    for Ω in angle_tests, ω in angle_tests, ν in angle_tests
        kep = [sma, ecc, pi/3, Ω, ω, ν]
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        @test isapproxvec_percent(kep, kep2; tol = 1e-13)
        @test isapproxvec_percent(cart, kep_to_cart(kep2, μ); tol = 1e-13)
    end
end

# Where an angle is undefined, the conventions are those of `_orbit_orientation`: equatorial,
# Ω = 0 and ω becomes the longitude of periapsis; circular, ω = 0 and ν becomes the argument of
# latitude; circular equatorial, Ω = ω = 0 and ν becomes the true longitude.
@testset "cart to kep Special Case Orbits Tests [AS-5, AS-6]" begin
    wrap(x) = mod(x, 2pi)
    hyperbolic = [-6000.0, 2.16, pi/10, pi/20, 3pi/2, 3.5pi/2]
    elliptic   = [10000.0, 0.1, 1.0, 2.0, 3.0, 2pi/7]
    cases = [
        # label, elements in, elements expected back
        (:elliptic_inclined,   elliptic, elliptic),
        (:hyperbolic,          hyperbolic, hyperbolic),
        (:elliptic_equatorial, [7000.0, 0.08, 0.0, 0.5, 3pi/2, pi/7],
                               [7000.0, 0.08, 0.0, 0.0, wrap(0.5 + 3pi/2), pi/7]),
        (:circular_inclined,   [10000.0, 0.0, 1.0, 0.1, 0.2, 2pi/7],
                               [10000.0, 0.0, 1.0, 0.1, 0.0, 0.2 + 2pi/7]),
        (:circular_equatorial, [10000.0, 0.0, 0.0, 2.0, 4.0, pi/7],
                               [10000.0, 0.0, 0.0, 0.0, 0.0, wrap(2.0 + 4.0 + pi/7)]),
    ]
    @testset "$label" for (label, kep, expected) in cases
        cart = kep_to_cart(kep, μ)
        kep2 = cart_to_kep(cart, μ)
        @test isapproxvec_percent(expected, kep2; tol = 1e-12)
        @test isapproxvec_percent(cart, kep_to_cart(kep2, μ); tol = 1e-12)
    end
end
