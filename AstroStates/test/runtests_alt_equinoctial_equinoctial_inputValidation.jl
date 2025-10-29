using Test
using AstroStates

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

@testset "alt_equinoctial_to_equinoctial / equinoctial_to_alt_equinoctial input validation" begin
    # alt_equinoctial_to_equinoctial: wrong length
    test_throws_msg(alt_equinoctial_to_equinoctial, bad7;
        substr = "Input must be a 6-element alternate equinoctial vector: [a, h, k, altp, altq, λ]")

    # equinoctial_to_alt_equinoctial: wrong length
    test_throws_msg(equinoctial_to_alt_equinoctial, bad7;
        substr = "Input must be a 6-element equinoctial vector: [a, h, k, p, q, λ]")

    # Test singularity: inclination ≈ π (180°) in alt_equinoctial_to_equinoctial
    # When altp^2 + altq^2 ≈ 1, then i = 2*asin(1) = π
    singular_alt_eq = [7000.0, 0.01, 0.0, 0.99999999999999999, 0.0, π/4]  # altp ≈ 1, so i ≈ π
    result = alt_equinoctial_to_equinoctial(singular_alt_eq)
    @test all(isnan.(result))

    # Test singularity: inclination ≈ π (180°) in equinoctial_to_alt_equinoctial
    # When p^2 + q^2 is very large, then i = 2*atan(large) ≈ π
    # Need value that gets within the default tolerance of 1e-12
    large_pq = 1e13  # Makes i = 2*atan(large) even closer to π
    singular_eq = [7000.0, 0.01, 0.0, large_pq, 0.0, π/4]
    result = equinoctial_to_alt_equinoctial(singular_eq)
    @test all(isnan.(result))

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_alt_eq = [7000.0, 0.01, 0.0, 0.05, 0.02, π/4]
    good_eq = alt_equinoctial_to_equinoctial(good_alt_eq)
    @test length(good_eq) == 6
    back_alt_eq = equinoctial_to_alt_equinoctial(good_eq)
    @test length(back_alt_eq) == 6
end
nothing