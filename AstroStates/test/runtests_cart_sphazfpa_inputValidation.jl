using Test
using AstroStates

mu = 398600.4415
bad7 = collect(1.0:7.0)  # wrong length

# Helper to assert an error message contains substring
function test_throws_msg(f, args...; substr)
    err = try
        f(args...)
        error("Expected throw did not occur")
    catch e
        @test occursin(substr, sprint(showerror, e))
        e
    end
    return err
end

@testset "cart_to_sphazfpa / sphazfpa_to_cart input validation" begin
    # cart_to_sphazfpa: wrong length
    test_throws_msg(cart_to_sphazfpa, bad7;
        substr = "Input vector must have six elements: [x, y, z, vx, vy, vz]")

    # sphazfpa_to_cart: wrong length
    test_throws_msg(sphazfpa_to_cart, bad7;
        substr = "Input vector must have six elements: [r, λ, δ, v, αₚ, ψ]")

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]
    good_sphazfpa  = cart_to_sphazfpa(good_cart)
    @test length(good_sphazfpa) == 6
    back_cart = sphazfpa_to_cart(good_sphazfpa)
    @test length(back_cart) == 6
end

# r < tol should warn and return NaNs
@testset "cart_to_sphazfpa: r below tol returns NaNs with warning" begin
    cart_r_small = [0.0, 0.0, 0.0, 0.1, 0.2, 0.3]  # r = 0, v > 0
    @test_logs (:warn, r"Conversion failed: Position magnitude r = .* below tolerance") begin
        out = cart_to_sphazfpa(cart_r_small; tol=1e-9)
        @test length(out) == 6
        @test all(isnan, out)
    end
end

# v < tol should warn and return NaNs
@testset "cart_to_sphazfpa: v below tol returns NaNs with warning" begin
    cart_v_small = [7000.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # r > 0, v = 0
    @test_logs (:warn, r"Conversion failed: Velocity magnitude v = .* below tolerance") begin
        out = cart_to_sphazfpa(cart_v_small; tol=1e-9)
        @test length(out) == 6
        @test all(isnan, out)
    end
end