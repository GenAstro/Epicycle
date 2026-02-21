using Test
using AstroStates
using LinearAlgebra

@testset "cart_to_equinoctial Edge Cases and Warnings" begin
    μ = 398600.4415
    tol = 1e-12
    
    @testset "Invalid input length" begin
        # Line 39: Wrong number of elements
        cart_short = [1.0, 2.0, 3.0]
        @test_throws ErrorException AstroStates.cart_to_equinoctial(cart_short, μ; tol=tol)
        
        cart_long = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        @test_throws ErrorException AstroStates.cart_to_equinoctial(cart_long, μ; tol=tol)
    end
    
    @testset "Degenerate position vector (r < tol)" begin
        # Lines 49-50: Position magnitude less than tolerance
        cart = [1e-15, 0.0, 0.0, 0.0, 7.5, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Invalid gravitational parameter (μ < tol)" begin
        # Lines 53-54: μ less than tolerance
        cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, 1e-15; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Radial/degenerate orbit (h_mag < tol)" begin
        # Lines 62-63: Angular momentum near zero (radial orbit)
        # Position and velocity vectors parallel (radial trajectory)
        r = 10000.0
        v_radial = 5.0  # Velocity in same direction as position
        cart = [r, 0.0, 0.0, v_radial, 0.0, 0.0]  # r and v both along x-axis
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Hyperbolic orbit (e > 1 - tol)" begin
        # Lines 61-62: Parabolic or hyperbolic orbit
        # Create hyperbolic orbit
        r = 10000.0
        v_hyp = sqrt(2 * μ / r) * 1.5  # Hyperbolic velocity
        cart = [r, 0.0, 0.0, 0.0, v_hyp, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Parabolic orbit (e ≈ 1)" begin
        # Lines 61-62: Parabolic orbit case
        r = 10000.0
        v_para = sqrt(2 * μ / r)  # Escape velocity (parabolic)
        cart = [r, 0.0, 0.0, 0.0, v_para, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Near-zero periapsis (singular conic)" begin
        # Lines 70-71: Radius of periapsis less than tolerance
        # Need high eccentricity but still e < 1 - tol to avoid the earlier check
        # Using a state with very small periapsis but definitely elliptic
        # Construct using known elliptic elements with tiny rp
        
        # Small periapsis, large apoapsis (highly eccentric ellipse)
        rp_target = 1e-14  # Very small but > 0
        ra = 20000.0
        a = (rp_target + ra) / 2
        e = (ra - rp_target) / (ra + rp_target)
        
        # At periapsis: r = rp, velocity perpendicular
        r_vec = [rp_target, 0.0, 0.0]
        v_mag = sqrt(μ * (2 / rp_target - 1 / a))
        v_vec = [0.0, v_mag, 0.0]
        
        cart = [r_vec..., v_vec...]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Retrograde equatorial orbit (i ≈ π)" begin
        # Lines 80-81: Inclination close to π
        # Create retrograde equatorial orbit
        r = 7000.0
        v_circ = sqrt(μ / r)
        # Velocity in -y direction creates i ≈ π orbit
        cart = [r, 0.0, 0.0, 0.0, -v_circ, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Valid elliptic orbit" begin
        # Valid test case
        cart = [6778.0, 0.0, 0.0, 0.0, 7.66, 0.0]
        result = AstroStates.cart_to_equinoctial(cart, μ; tol=tol)
        
        @test length(result) == 6
        @test all(isfinite.(result))
        
        a, h, k, p, q, λ = result
        @test a > 0  # Semi-major axis positive
        @test sqrt(h^2 + k^2) < 1  # Eccentricity magnitude < 1
    end
end
