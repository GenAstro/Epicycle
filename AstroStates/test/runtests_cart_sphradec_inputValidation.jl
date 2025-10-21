using Test
using AstroStates

mu = 398600.4415
ad7 = collect(1.0:7.0)  # wrong length

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

@testset "cart_to_sphradec / sphradec_to_cart input validation" begin
    # cart_to_sphradec: wrong length
    test_throws_msg(cart_to_sphradec, bad7;
        substr = "Input vector must have six elements: [x, y, z, vx, vy, vz]")

    # sphradec_to_cart: wrong length
    test_throws_msg(sphradec_to_cart, bad7;
        substr = "Input vector must have exactly six elements: [r, λᵣ, δᵣ, v, λᵥ, δᵥ]")

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]
    good_sphradec  = cart_to_sphradec(good_cart)
    @test length(good_sphradec) == 6
    back_cart = sphradec_to_cart(good_sphradec)
    @test length(back_cart) == 6
end

# Targeted: r < tol should warn and return NaNs
@testset "cart_to_sphradec: r below tol returns NaNs with warning" begin
    cart_r_zero = [0.0, 0.0, 0.0, 0.1, 0.2, 0.3]  # r = 0, v > 0
    @test_logs (:warn, r"Conversion failed: Radius is zero\.") begin
        out = cart_to_sphradec(cart_r_zero; tol=1e-9)
        @test length(out) == 6
        @test all(isnan, out)
    end
end

# Targeted: v < tol should warn and return NaNs
@testset "cart_to_sphradec: v below tol returns NaNs with warning" begin
    cart_v_zero = [7000.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # r > 0, v = 0
    @test_logs (:warn, r"Conversion failed: Velocity is zero\.") begin
        out = cart_to_sphradec(cart_v_zero; tol=1e-9)
        @test length(out) == 6
        @test all(isnan, out)
    end
end