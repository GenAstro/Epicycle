using Test
using AstroStates
using LinearAlgebra

@testset "inasymptote_to_kep Edge Cases and Warnings" begin
    μ = 398600.4415
    tol = 1e-12
    
    @testset "Invalid input length" begin
        # Line 36: Wrong number of elements
        inasym_short = [7000.0, 10.0, 0.0]
        @test_throws ErrorException AstroStates.inasymptote_to_kep(inasym_short, μ; tol=tol)
        
        inasym_long = [7000.0, 10.0, 0.0, π/4, π/2, π/2, 0.0]
        @test_throws ErrorException AstroStates.inasymptote_to_kep(inasym_long, μ; tol=tol)
    end
    
    @testset "Nearly parabolic orbit (abs(C₃) < tol)" begin
        # Lines 50-51: C₃ ≈ 0
        rₚ = 7000.0
        c₃ = 1e-15  # Nearly zero
        λₐ = 0.0
        δₐ = π/4
        θᵦ = π/2
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Nearly circular orbit (e < tol)" begin
        # Lines 54-55: e ≈ 0
        # e = 1 - rₚ/a, where a = -μ/c₃
        # For e ≈ 0, need rₚ ≈ a, so rₚ ≈ -μ/c₃
        c₃ = 10.0
        a = -μ / c₃
        rₚ = a * (1 - 1e-15)  # Makes e very small
        λₐ = 0.0
        δₐ = π/4
        θᵦ = π/2
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Asymptote aligned with z-axis" begin
        # Lines 67-68: δₐ = ±π/2 makes asymptote point along z
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = 0.0
        δₐ = π/2  # Asymptote pointing straight up
        θᵦ = π/4
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        @test all(isnan.(result))
    end
    
    @testset "Elliptical orbit branch (c₃ <= -tol)" begin
        # Line 113: c₃ < 0 for elliptical orbit
        rₚ = 7000.0
        c₃ = -5.0  # Negative C₃ for elliptical
        λₐ = π/4
        δₐ = π/6
        θᵦ = π/3
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        # Should return valid result for elliptical orbit
        @test length(result) == 6
        @test all(isfinite.(result))
        
        # Calculate expected values
        a_expected = -μ / c₃
        e_expected = 1 - rₚ / a_expected
        
        a, e = result[1], result[2]
        @test a ≈ a_expected rtol=1e-10
        @test e ≈ e_expected rtol=1e-10
    end
    
    @testset "Hyperbolic orbit branch (c₃ > tol)" begin
        # Lines 115-117: c₃ > 0 for hyperbolic orbit
        rₚ = 7000.0
        c₃ = 10.0  # Positive C₃ for hyperbolic
        λₐ = π/4
        δₐ = π/6
        θᵦ = π/3
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        # Should return valid result for hyperbolic orbit
        @test length(result) == 6
        @test all(isfinite.(result))
        
        # Calculate expected values
        a_expected = -μ / c₃
        e_expected = 1 - rₚ / a_expected
        
        a, e = result[1], result[2]
        @test a ≈ a_expected rtol=1e-10
        @test e ≈ e_expected rtol=1e-10
    end
    
    
    @testset "Equatorial prograde orbit (i < tol)" begin
        # Lines 119, 121-123: e >= tol and i < tol
        # For i = 0, need ĥ = ẑ
        # With δₐ = 0, λₐ = 0: ŝ = [1, 0, 0]
        # Then Ê = [0, 1, 0] and N̂ = [0, 0, 1] = ẑ
        # With θᵦ = π/2, ami = 0, so ĥ = cos(0)*N̂ = ẑ
        
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = 0.0
        δₐ = 0.0  # Asymptote in equatorial plane
        θᵦ = π/2  # Makes ĥ = ẑ
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        @test all(isfinite.(result))
        i = result[3]
        Ω = result[4]
        ω = result[5]
        @test i < tol  # Should be nearly zero
        @test Ω ≈ 0.0 atol=1e-10  # Lines 121
        @test isfinite(ω)  # Lines 122-123 compute ω
    end
    
    @testset "Equatorial retrograde orbit (i >= π - tol)" begin
        # Lines 113, 115-117: e >= tol and i >= π - tol
        # For i = π, need ĥ = -ẑ
        # With δₐ = 0, λₐ = 0: ŝ = [1, 0, 0]
        # Then Ê = [0, 1, 0] and N̂ = [0, 0, 1] = ẑ
        # With θᵦ = -π/2, ami = π, so ĥ = cos(π)*N̂ = -ẑ
        
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = 0.0
        δₐ = 0.0  # Asymptote in equatorial plane
        θᵦ = -π/2  # Makes ĥ = -ẑ (retrograde)
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        @test all(isfinite.(result))
        i = result[3]
        Ω = result[4]
        ω = result[5]
        @test i >= π - tol  # Should be nearly π (line 113)
        @test Ω ≈ 0.0 atol=1e-10  # Line 115
        @test isfinite(ω)  # Lines 116-117 compute ω with -acos
    end
    
    @testset "General inclined orbit" begin
        # Line 102: General inclined case (e >= tol && i >= tol && i < π - tol)
        # Use typical inclined hyperbolic orbit
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = π/4
        δₐ = π/6  # Gives inclined orbit
        θᵦ = π/4  # Non-zero, non-extreme value
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        @test all(isfinite.(result))
        i = result[3]
        Ω = result[4]
        ω = result[5]
        # Should be inclined (not equatorial)
        @test i >= tol && i < π - tol
        @test isfinite(Ω) && isfinite(ω)
    end
    
    @testset "General inclined with ω quadrant adjustment" begin
        # Line 113: ω = ê[3] < 0 ? 2π - ω : ω
        # Need eccentricity vector with ê[3] < 0
        # This happens when ω is in 3rd or 4th quadrant (π < ω < 2π)
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = π
        δₐ = -π/6  # Negative declination
        θᵦ = π/6
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        @test all(isfinite.(result))
        # Just verify it executes successfully
        @test length(result) == 6
    end
    
    @testset "Equatorial prograde with ω quadrant adjustment" begin
        # Line 119: ω = ê[2] < 0 ? 2π - ω : ω (for equatorial prograde)
        # Need ê[2] < 0, which means eccentricity vector pointing in -y direction
        # For equatorial prograde (i ≈ 0), use λₐ = π to flip asymptote direction
        rₚ = 7000.0
        c₃ = 10.0
        λₐ = π  # Asymptote in -x direction
        δₐ = 0.0  # Equatorial
        θᵦ = π/2  # Makes ĥ = ẑ (prograde)
        ν = π/2
        
        inasym = [rₚ, c₃, λₐ, δₐ, θᵦ, ν]
        result = AstroStates.inasymptote_to_kep(inasym, μ; tol=tol)
        
        @test all(isfinite.(result))
        i = result[3]
        @test i < tol  # Equatorial prograde
        # Just verify it executes the quadrant adjustment
        @test length(result) == 6
    end
end
