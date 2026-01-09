using Test
using LinearAlgebra

# Test round-trip conversions and verify against built-in conversions
@testset "to_posvel/set_posvel! round-trip" begin
    state_types = [
        (KeplerianState(7000.0, 0.001, 0.1, 0.2, 0.3, 0.4), Keplerian()),
        (ModifiedEquinoctialState(7000.0, 0.001, 0.002, 0.003, 0.004, 0.5), ModifiedEquinoctial()),
        (CartesianState([7000.0, 100.0, 200.0, 0.5, 7.5, 0.1]), Cartesian()),
        (EquinoctialState(7000.0, 0.001, 0.002, 0.003, 0.004, 0.5), Equinoctial()),
        (SphericalRADECState(7000.0, 0.1, 0.2, 7.5, 0.01, 0.02), SphericalRADEC()),
    ]
    
    for (initial_state, expected_type) in state_types
        sc = Spacecraft(
            state = initial_state,
            coord_sys = CoordinateSystem(earth, ICRFAxes())
        )
        
        # Get original posvel
        original_posvel = to_posvel(sc)
        
        # For Keplerian, verify exact match with built-in conversion
        if initial_state isa KeplerianState
            expected_cart = CartesianState(initial_state, earth.mu)
            @test original_posvel ≈ expected_cart.posvel atol=1e-10
        end
        
        # Set it back
        set_posvel!(sc, original_posvel)
        
        # Should get the same thing back
        final_posvel = to_posvel(sc)
        @test final_posvel ≈ original_posvel rtol=1e-12
        
        # State type should be preserved
        @test sc.state.statetype == expected_type
    end
    
    # Test Cartesian fast path with exact equality
    sc_cart = Spacecraft(state = CartesianState([7000.0, 100.0, 200.0, 0.5, 7.5, 0.1]))
    @test to_posvel(sc_cart) == [7000.0, 100.0, 200.0, 0.5, 7.5, 0.1]
    set_posvel!(sc_cart, [7100.0, 50.0, 100.0, 0.3, 7.6, 0.2])
    @test sc_cart.state.state == [7100.0, 50.0, 100.0, 0.3, 7.6, 0.2]
end

# Test multiple conversions
@testset "Multiple state type conversions" begin
    # Start with Keplerian
    sc = Spacecraft(
        state = KeplerianState(7000.0, 0.01, 0.2, 0.3, 0.4, 0.5),
        coord_sys = CoordinateSystem(earth, ICRFAxes())
    )
    
    # Convert through Cartesian several times
    for i in 1:3
        pv = to_posvel(sc)
        pv[1] += 10.0  # Small perturbation
        set_posvel!(sc, pv)
        
        # State type should remain Keplerian
        @test sc.state.statetype == Keplerian()
    end
    
    # Final state should be different from initial
    final_kep = get_state(sc, Keplerian())
    @test final_kep.sma > 7000.0
end

# Test error handling
@testset "to_posvel/set_posvel! error cases" begin
    # Test set_posvel! with wrong length vector
    sc = Spacecraft(state = KeplerianState(7000.0, 0.001, 0.1, 0.2, 0.3, 0.4))
    @test_throws ArgumentError set_posvel!(sc, [1.0, 2.0, 3.0])  # Too short
    @test_throws ArgumentError set_posvel!(sc, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])  # Too long
    
    # Test error when μ is required but coord system origin is not a celestial body
    # Use a spacecraft as the origin (has no mu field)
    origin_sc = Spacecraft(state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]))
    coords_no_mu = CoordinateSystem(origin_sc, ICRFAxes())
    
    # to_posvel: Keplerian -> Cartesian requires μ
    sc_kep = Spacecraft(
        state = KeplerianState(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4),
        coord_sys = coords_no_mu
    )
    @test_throws ErrorException to_posvel(sc_kep)
    
    # set_posvel!: Cartesian -> Keplerian requires μ
    sc_kep2 = Spacecraft(
        state = KeplerianState(7000.0, 0.01, 0.1, 0.2, 0.3, 0.4),
        coord_sys = coords_no_mu
    )
    @test_throws ErrorException set_posvel!(sc_kep2, [7100.0, 50.0, 100.0, 0.3, 7.6, 0.2])
end

# Test numeric type preservation
@testset "Numeric type preservation in conversions" begin
    # Float64 spacecraft
    sc = Spacecraft(
        state = KeplerianState(7000.0, 0.001, 0.1, 0.2, 0.3, 0.4),
        mass = 1500.0,
        coord_sys = CoordinateSystem(earth, ICRFAxes())
    )
    
    new_pv = [7100.0, 50.0, 100.0, 0.3, 7.6, 0.2]
    set_posvel!(sc, new_pv)
    
    # State should still be Float64
    @test eltype(sc.state.state) === Float64
    @test typeof(sc.mass) === Float64
end
