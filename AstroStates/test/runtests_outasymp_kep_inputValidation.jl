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

@testset "out_to_kep / out_to_cart input validation" begin
    # out_to_kep: wrong length
    test_throws_msg(outasymptote_to_kep, bad7, mu;
        substr = "Input must be a 6-element vector: [a, e, i, Ω, ω, ν]")

    # Positive controls (shape only; numeric content not asserted deeply here)
    good_out = [6400.000000000003, -49.8250551875, deg2rad(250.633213963147),
                deg2rad(-11.10476195323574), deg2rad(96.50681841871349), deg2rad(300.0)]
    back_kep = outasymptote_to_kep(good_out, mu)
    @test length(back_kep) == 6
end

# Nearly parabolic c3 ~ 0: force branch with very large tol
@testset "outasymptote_to_kep: near-parabolic warns and returns NaNs" begin
    good_out = [6400.0, 0.0, deg2rad(10.0), deg2rad(10.0), deg2rad(20.0), deg2rad(30.0)]
    @test_logs (:warn, r"nearly parabolic") begin
        v = outasymptote_to_kep(good_out, mu; tol=Inf)
        @test length(v) == 6
        @test all(isnan, v)
    end
end

# Asymptote vector aligned with z-axis → warn + NaNs
@testset "outasymptote_to_kep: asymptote aligned with z-axis warns and returns NaNs" begin
    aligned_out = [6400.0, -50.0, 0.0, π/2, 0.0, 0.0]  # set DEC = ±π/2 to align ŝ with ẑ
    @test_logs (:warn, r"Asymptote vector is aligned with z-axis") begin
        v = outasymptote_to_kep(aligned_out, mu)
        @test length(v) == 6
        @test all(isnan, v)
    end
end

# Equatorial prograde branch (i ≈ 0): expect Ω set to 0.0
@testset "outasymptote_to_kep: equatorial prograde (i≈0) sets Ω=0" begin
    # Build a hyperbolic, equatorial, prograde Cartesian state at periapsis, then derive out-asymptote
    rp = 7000.0
    e = 1.3
    vperi = sqrt(mu * (1 + e) / rp)
    cart_pro = [rp, 0.0, 0.0, 0.0, vperi, 0.0]            # h along +z ⇒ i≈0 (prograde)
    out_pro = cart_to_outasymptote(cart_pro, mu)
    kep = outasymptote_to_kep(out_pro, mu)
    @test length(kep) == 6
    a, e_out, i_out, Ω_out, ω_out, ν_out = kep
    @test isapprox(i_out, 0.0; atol=1e-8)
    @test Ω_out == 0.0
    @test 0.0 ≤ ω_out ≤ 2π
end

# Nearly circular: pass in a circular orbit (e ≈ 0) → warn and return NaNs
#=
@testset "outasymptote_to_kep: nearly circular warns and returns NaNs" begin
    r = 7000.0
    v_circ = sqrt(mu / r)
    cart_circ = [r, 0.0, 0.0, 0.0, 0.0, v_circ]  # circular, equatorial, prograde
    out_circ = cart_to_outasymptote(cart_circ, mu)
    @test_logs (:warn, r"Conversion failed: Orbit is nearly circular\.") begin
        v = outasymptote_to_kep(out_circ, mu)  # default tol; e ≈ 0 should trigger
        @test length(v) == 6
        @test all(isnan, v)
    end
end


# Line-of-nodes undefined: straight-line trajectory (r ∥ v) → h = 0 ⇒ n = 0
@testset "outasymptote_to_kep: line-of-nodes undefined warns and returns NaNs" begin
    cart_line = [7000.0, 0.0, 0.0, 1.0, 0.0, 0.0]  # r × v = 0 ⇒ undefined plane
    out_line = cart_to_outasymptote(cart_line, mu)
    @test_logs (:warn, r"Conversion failed: Line-of-nodes is undefined\.") begin
        v = outasymptote_to_kep(out_line, mu)
        @test length(v) == 6
        @test all(isnan, v)
    end
end

=#

# Equatorial retrograde branch (i ≈ π): expect Ω set to 0.0
#=
@testset "outasymptote_to_kep: equatorial retrograde (i≈π) sets Ω=0" begin
    rp = 7000.0
    e = 1.3
    vperi = sqrt(mu * (1 + e) / rp)
    cart_ret = [rp, 0.0, 0.0, 0.0, -vperi, 0.0]           # h along -z ⇒ i≈π (retrograde)
    out_ret = cart_to_outasymptote(cart_ret, mu)
    kep = outasymptote_to_kep(out_ret, mu)
    @test length(kep) == 6
    a, e_out, i_out, Ω_out, ω_out, ν_out = kep
    @test isapprox(i_out, π; atol=1e-8)
    @test Ω_out == 0.0
    @test 0.0 ≤ ω_out ≤ 2π
end
=#