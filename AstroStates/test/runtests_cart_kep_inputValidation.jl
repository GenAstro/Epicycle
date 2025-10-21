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

@testset "cart_to_kep / kep_to_cart input validation" begin
    # cart_to_kep: wrong length
    test_throws_msg(cart_to_kep, bad7, mu;
        substr = "Input vector must have exactly six elements: [x, y, z, vx, vy, vz].")

    # kep_to_cart: wrong length
    test_throws_msg(kep_to_cart, bad7, mu;
        substr = "Input vector must have exactly six elements: a, e, i, Ω, ω, ν.")

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]
    good_kep  = cart_to_kep(good_cart, mu)
    @test length(good_kep) == 6
    back_cart = kep_to_cart(good_kep, mu)
    @test length(back_cart) == 6
end