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

@testset "kep_to_modkep / modkep_to_kep input validation" begin
    # kep_to_modkep: wrong length
    test_throws_msg(kep_to_modkep, bad7;
        substr = "Input vector must contain six elements: [a, e, i, Ω, ω, ν]")

    # modkep_to_kep: wrong length
    test_throws_msg(modkep_to_kep, bad7;
        substr = "Input vector must contain six elements: [rₚ, rₐ, i, Ω, ω, ν]")

    # Test parabolic orbit (e ≈ 1)
    parabolic_kep = [7000.0, 1.0, π/4, 0.0, 0.0, π/3]
    result = kep_to_modkep(parabolic_kep)
    @test all(isnan.(result))

    # Test negative eccentricity
    neg_ecc_kep = [7000.0, -0.1, π/4, 0.0, 0.0, π/3]
    result = kep_to_modkep(neg_ecc_kep)
    @test all(isnan.(result))

    # Test incompatible a and e: positive a with hyperbolic e
    bad_ae_kep = [7000.0, 1.5, π/4, 0.0, 0.0, π/3]
    result = kep_to_modkep(bad_ae_kep)
    @test all(isnan.(result))

    # Test incompatible a and e: negative a with elliptic e
    bad_ae_kep2 = [-7000.0, 0.5, π/4, 0.0, 0.0, π/3]
    result = kep_to_modkep(bad_ae_kep2)
    @test all(isnan.(result))

    # Test singular conic (rp ≈ 0)
    singular_kep = [1e-15, 1.0 - 1e-15, π/4, 0.0, 0.0, π/3]  # rp = a*(1-e) ≈ 0
    result = kep_to_modkep(singular_kep)
    @test all(isnan.(result))

    # modkep_to_kep specific tests to cover untested lines

    # Test zero radius of apoapsis (abs(rₐ) < tol)
    zero_ra_modkep = [7000.0, 0.0, π/4, 0.0, 0.0, π/3]
    result = modkep_to_kep(zero_ra_modkep)
    @test all(isnan.(result))

    # Test inconsistent Modified Keplerian state: rₐ < rₚ but rₐ > 0
    inconsistent_modkep = [8000.0, 7000.0, π/4, 0.0, 0.0, π/3]  # rₐ < rₚ but both positive
    result = modkep_to_kep(inconsistent_modkep)
    @test all(isnan.(result))

    # Test singular conic: rₚ ≤ tol
    singular_rp_modkep = [1e-15, 7000.0, π/4, 0.0, 0.0, π/3]
    result = modkep_to_kep(singular_rp_modkep)
    @test all(isnan.(result))

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_kep = [7000.0, 0.01, π/4, 0.0, 0.0, π/3]
    good_modkep = kep_to_modkep(good_kep)
    @test length(good_modkep) == 6
    back_kep = modkep_to_kep(good_modkep)
    @test length(back_kep) == 6
end
nothing