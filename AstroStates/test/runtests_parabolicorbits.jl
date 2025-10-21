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

    # Also exercise state constructor if available
    try
        cstate = CartesianState(cart, mu)
        kstate = KeplerianState(cstate, mu)
        @test kstate.ecc ≈ 1.0 atol=1e-10
    catch e
        @info "State type constructors not applicable for parabolic case: $(e)"
    end
end
nothing