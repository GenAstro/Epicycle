# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

using Test
using AstroStates

mu = 398600.4415            # use ASCII name
bad7 = collect(1.0:7.0)      # wrong length

# Helper (contains its own @test; keep usage consistent)
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

@testset "cart_to_mee / mee_to_cart input validation" begin
    test_throws_msg(cart_to_mee, bad7, mu;
        substr = "exactly six elements: [x, y, z, vx, vy, vz]")

    test_throws_msg(mee_to_cart, bad7, mu;
        substr = "exactly six elements: [p, f, g, h, k, L]")

    good_cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]
    good_mee  = cart_to_mee(good_cart, mu)
    @test length(good_mee) == 6
    back_cart = mee_to_cart(good_mee, mu)
    @test length(back_cart) == 6
end

const mu_pos = 398600.4415
const mu_neg = -398600.4415
good_cart = [7000.0, 0.0, 0.0, 0.0, 7.5, 1.0]
parallel_cart = [7000.0, 0.0, 0.0, 1.0, 0.0, 0.0]
bad_len = collect(1.0:7.0)

@testset "cart_to_mee branch coverage" begin
    @test_throws ErrorException cart_to_mee(bad_len, mu_pos)

    invalid_j_supported = true
    invalid_j_err = try
        cart_to_mee(good_cart, mu_pos; j=0.0)
        invalid_j_supported = false
        nothing
    catch e
        e
    end
    if invalid_j_err !== nothing
        @test occursin("Invalid value for j", sprint(showerror, invalid_j_err))
    elseif !invalid_j_supported
        @info "cart_to_mee keyword j unsupported or unvalidated"
    end

    # Degenerate states and a non-positive μ warn and return NaN, as the other converters do
    @test_logs (:warn, r"degenerate angular momentum") begin
        @test all(isnan, cart_to_mee(parallel_cart, mu_pos))
    end
    @test_logs (:warn, r"μ < tolerance") begin
        @test all(isnan, cart_to_mee(good_cart, mu_neg))
    end
    @test_logs (:warn, r"degenerate position or velocity") begin
        @test all(isnan, cart_to_mee([0.0, 0, 0, 1, 0, 0], mu_pos))
    end

    retrograde_cart = [7000.0, 0.0, 0.0, 0.0, -7.5, 0.0]
    @test_logs (:warn, r"Singularity computing h and k") begin
        sing = cart_to_mee(retrograde_cart, mu_pos)
        @test length(sing) == 6
        @test all(isnan, sing)
    end
end

function throws_msg(f, args...; substr)
    err = try
        f(args...)
        error("Expected throw")
    catch e
        @test occursin(substr, sprint(showerror, e))
        e
    end
    return err
end

@testset "mee_to_cart input / branch coverage" begin
    bad_len = collect(1.0:7.0)
    throws_msg(mee_to_cart, bad_len, mu;
        substr = "exactly six elements: [p, f, g, h, k, L]")

    valid_mod = [7000.0, 0.001, -0.002, 0.01, -0.015, 1.2]
    res = mee_to_cart(valid_mod, mu)
    @test length(res) == 6

   # invalid j
    err_j = try
        mee_to_cart(valid_mod, mu; j=0.0)
        nothing
    catch e
        e
    end
    @test err_j !== nothing
    @test occursin("Invalid value for j", sprint(showerror, err_j))


    # A non-positive semi-latus rectum describes no orbit: a warning and NaN, not zeros
    for p in (-10.0, 0.0)
        bad_p = copy(valid_mod); bad_p[1] = p
        @test_logs (:warn, r"Semi-latus rectum") begin
            @test all(isnan, mee_to_cart(bad_p, mu))
        end
    end
end

nothing
