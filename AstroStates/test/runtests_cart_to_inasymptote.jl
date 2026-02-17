using Test
using AstroStates
using LinearAlgebra

@testset "cart_to_inasymptote Edge Cases and Warnings" begin
    μ = 398600.4415
    tol = 1e-12
    
    @testset "Invalid input length" begin
        # Test with wrong number of elements
        cart_short = [1.0, 2.0, 3.0]
        @test_throws ErrorException AstroStates.cart_to_inasymptote(cart_short, μ; tol=tol)
        
        cart_long = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
        @test_throws ErrorException AstroStates.cart_to_inasymptote(cart_long, μ; tol=tol)
    end
    
    @testset "Degenerate position vector (r < tol)" begin
        # Position vector with very small magnitude
        cart = [1e-15, 0.0, 0.0, 0.0, 12.0, 0.0]
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Degenerate velocity vector (v < tol)" begin
        # Velocity vector with very small magnitude
        cart = [10000.0, 0.0, 0.0, 1e-15, 0.0, 0.0]
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Zero angular momentum (h < tol)" begin
        # Radial trajectory - position and velocity are collinear
        cart = [10000.0, 0.0, 0.0, 12.0, 0.0, 0.0]
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Circular orbit (e < tol)" begin
        # Circular orbit with e ≈ 0
        r = 7000.0
        v_circ = sqrt(μ / r)
        cart = [r, 0.0, 0.0, 0.0, v_circ, 0.0]
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Parabolic orbit (C₃ ≈ 0)" begin
        # Parabolic trajectory with C₃ ≈ 0
        r = 10000.0
        v_para = sqrt(2 * μ / r)  # Escape velocity
        cart = [r, 0.0, 0.0, 0.0, v_para, 0.0]
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Valid hyperbolic incoming asymptote" begin
        # Valid test case from existing tests
        truth_cart = [4486.62555067877, -2278.090328380258, 8463.948027345208, 
                      -8.831137842558398, 2.38254720654504, -0.1278435615976432]
        result = AstroStates.cart_to_inasymptote(truth_cart, μ; tol=tol)
        
        # Should return valid 6-element vector
        @test length(result) == 6
        @test all(isfinite.(result))
        
        # Check that values are reasonable
        rₚ, C₃, λₐ, δₐ, θᵦ, ν = result
        @test rₚ > 0  # Periapsis radius should be positive
        @test C₃ > 0  # Hyperbolic orbit has positive C₃
        @test 0 ≤ λₐ ≤ 2π  # Right ascension in [0, 2π]
        @test -π/2 ≤ δₐ ≤ π/2  # Declination in [-π/2, π/2]
        @test 0 ≤ θᵦ ≤ 2π  # B-plane angle in [0, 2π]
        @test 0 ≤ ν ≤ 2π  # True anomaly in [0, 2π]
    end
    
    @testset "True anomaly quadrant handling (ν adjustment at end)" begin
        # Test case where dot(v̄, r̄) < 0 to trigger ν = 2π - ν
        # This happens when velocity and position are in opposite hemispheres
        # (velocity pointing backwards relative to position)
        
        # Create a hyperbolic state after periapsis passage
        # where velocity has negative radial component
        r_vec = [10000.0, 5000.0, 2000.0]
        # Velocity pointing somewhat opposite to position direction
        v_vec = [-8.0, 3.0, 1.0]
        cart = [r_vec..., v_vec...]
        
        result = AstroStates.cart_to_inasymptote(cart, μ; tol=tol)
        
        if all(isfinite.(result))
            ν = result[6]
            # Should be in valid range
            @test 0 ≤ ν ≤ 2π
        end
    end
    
    @testset "Negative angle wraparound (θᵦ < 0 and λₐ < 0)" begin
        # To hit line 96 (λₐ < 0 → λₐ + 2π): 
        # Create state with rla in range (π, 2π) which maps to negative in atan
        
        rp = 7000.0
        c3 = 10.0
        rla = 3π/2  # This is in (π, 2π), should trigger λₐ < 0 before wraparound
        dla = deg2rad(20.0)
        bpa = π/4
        ta = deg2rad(45.0)
        
        inasym1 = IncomingAsymptoteState(rp, c3, rla, dla, bpa, ta)
        cart1 = AstroStates.to_vector(CartesianState(inasym1, μ))
        result1 = AstroStates.cart_to_inasymptote(cart1, μ; tol=tol)
        
        if all(isfinite.(result1))
            λₐ = result1[3]
            @test 0 ≤ λₐ ≤ 2π  # Should be wrapped to [0, 2π]
        end
        
        # Another test with rla in 3rd quadrant
        rla2 = 5π/4
        inasym2 = IncomingAsymptoteState(rp, c3, rla2, dla, bpa, ta)
        cart2 = AstroStates.to_vector(CartesianState(inasym2, μ))
        result2 = AstroStates.cart_to_inasymptote(cart2, μ; tol=tol)
        
        if all(isfinite.(result2))
            λₐ = result2[3]
            @test 0 ≤ λₐ ≤ 2π
        end
        
        # To hit line 95 (θᵦ < 0 → θᵦ + 2π):
        # Try various B-plane angles that might result in negative atan result
        
        bpa3 = 3π/2  # B-plane angle in range that might give negative θᵦ
        inasym3 = IncomingAsymptoteState(rp, c3, π/4, dla, bpa3, ta)
        cart3 = AstroStates.to_vector(CartesianState(inasym3, μ))
        result3 = AstroStates.cart_to_inasymptote(cart3, μ; tol=tol)
        
        if all(isfinite.(result3))
            θᵦ = result3[5]
            @test 0 ≤ θᵦ ≤ 2π  # Should be wrapped to [0, 2π]
        end
        
        # Another B-plane angle test
        bpa4 = 7π/4
        inasym4 = IncomingAsymptoteState(rp, c3, π/3, dla, bpa4, ta)
        cart4 = AstroStates.to_vector(CartesianState(inasym4, μ))
        result4 = AstroStates.cart_to_inasymptote(cart4, μ; tol=tol)
        
        if all(isfinite.(result4))
            θᵦ = result4[5]
            @test 0 ≤ θᵦ ≤ 2π
        end
    end
end
