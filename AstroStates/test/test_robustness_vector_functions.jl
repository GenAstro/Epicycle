# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Every unphysical branch of the vector functions warns and returns six NaNs of the input's
# element type [AS-7], and every malformed input throws naming what it received [AS-8].

using Test
using AstroStates
using ForwardDiff

let μ = 398600.4415
    good     = [7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]
    zero_r   = [0.0, 0.0, 0.0, 0.0, 7.5, 0.0]
    zero_v   = [7000.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    radial   = [7000.0, 0.0, 0.0, 7.5, 0.0, 0.0]
    retro_eq = [7000.0, 0.0, 0.0, 0.0, -7.5, 0.0]                  # i = π
    circular = [7000.0, 0.0, 0.0, 0.0, sqrt(μ / 7000.0), 0.0]
    parabola = [7000.0, 0.0, 0.0, 0.0, sqrt(2μ / 7000.0), 0.0]
    hyperbola = [7000.0, 0.0, 0.0, 0.0, 12.0, 0.0]
    # A hyperbola in the x-z plane whose asymptotes lie along ±z: e = 2 puts the asymptotes
    # 120° from periapsis, so periapsis at ω = -30° and ω = 210° aims the outgoing and the
    # incoming asymptote along +z.
    zaxis_out = kep_to_cart([-7000.0, 2.0, π/2, 0.0, 11π/6, 0.0], μ)
    zaxis_in  = kep_to_cart([-7000.0, 2.0, π/2, 0.0, 7π/6, 0.0], μ)

    # (label, function of the state vector, state vector, expected warning)
    cases = [
        ("cart_to_kep μ",          c -> cart_to_kep(c, 0.0), good, r"μ < tolerance"),
        ("cart_to_kep r",          c -> cart_to_kep(c, μ), zero_r, r"degenerate position or velocity"),
        ("cart_to_kep h",          c -> cart_to_kep(c, μ), radial, r"degenerate angular momentum"),
        ("kep_to_cart μ",          k -> kep_to_cart(k, 0.0), [7000.0, 0.01, 0.1, 0.0, 0.0, 0.0], r"μ < tolerance"),
        ("kep_to_cart parabolic",  k -> kep_to_cart(k, μ), [7000.0, 1.0, 0.1, 0.0, 0.0, 0.0], r"parabolic or singular"),
        ("kep_to_cart asymptote",  k -> kep_to_cart(k, μ), [-7000.0, 2.0, 0.1, 0.0, 0.0, 2.5], r"asymptote"),
        ("kep_to_cart e < 0",      k -> kep_to_cart(k, μ), [7000.0, -0.1, 0.1, 0.0, 0.0, 0.0], r"cannot be negative"),
        ("kep_to_cart a > 0, e > 1", k -> kep_to_cart(k, μ), [7000.0, 1.5, 0.1, 0.0, 0.0, 0.0], r"cannot be positive"),
        ("kep_to_cart a < 0, e < 1", k -> kep_to_cart(k, μ), [-7000.0, 0.5, 0.1, 0.0, 0.0, 0.0], r"cannot be negative if"),
        ("kep_to_cart i < 0",      k -> kep_to_cart(k, μ), [7000.0, 0.1, -0.1, 0.0, 0.0, 0.0], r"Inclination .* outside"),
        ("kep_to_cart i > π",      k -> kep_to_cart(k, μ), [7000.0, 0.1, 3.2, 0.0, 0.0, 0.0], r"Inclination .* outside"),
        ("kep_to_modkep i < 0",    k -> kep_to_modkep(k), [7000.0, 0.1, -0.1, 0.0, 0.0, 0.0], r"Inclination .* outside"),
        ("modkep_to_kep i < 0",    m -> modkep_to_kep(m), [7000.0, 8000.0, -0.1, 0.0, 0.0, 0.0], r"Inclination .* outside"),
        # e = 1.5 puts the asymptote at 131.8°; ν = 2.5 rad is 143.2°, beyond it
        ("kep_to_modkep beyond asymptote", k -> kep_to_modkep(k), [-20000.0, 1.5, 0.5, 0.3, 0.4, 2.5], r"beyond the asymptote"),
        ("modkep_to_kep beyond asymptote", m -> modkep_to_kep(m), [10000.0, -50000.0, 0.5, 0.3, 0.4, 2.5], r"beyond the asymptote"),
        ("cart_to_mee μ",          c -> cart_to_mee(c, 0.0), good, r"μ < tolerance"),
        ("cart_to_mee r",          c -> cart_to_mee(c, μ), zero_r, r"degenerate position or velocity"),
        ("cart_to_mee h",          c -> cart_to_mee(c, μ), radial, r"degenerate angular momentum"),
        ("cart_to_mee i = π",      c -> cart_to_mee(c, μ), retro_eq, r"modified equinoctial elements are singular"),
        ("mee_to_cart μ",          m -> mee_to_cart(m, 0.0), [7000.0, 0.01, 0.0, 0.1, 0.0, 0.0], r"μ < tolerance"),
        ("mee_to_cart p",          m -> mee_to_cart(m, μ), [-1.0, 0.0, 0.0, 0.0, 0.0, 0.0], r"Semi-latus rectum"),
        ("mee_to_cart asymptote",  m -> mee_to_cart(m, μ), [7000.0, 2.0, 0.0, 0.0, 0.0, π], r"asymptote"),
        ("cart_to_equinoctial r",  c -> cart_to_equinoctial(c, μ), zero_r, r"Position magnitude"),
        ("cart_to_equinoctial μ",  c -> cart_to_equinoctial(c, 0.0), good, r"Gravitational parameter"),
        ("cart_to_equinoctial h",  c -> cart_to_equinoctial(c, μ), radial, r"Angular momentum near zero"),
        ("cart_to_equinoctial hyperbolic", c -> cart_to_equinoctial(c, μ), hyperbola, r"parabolic or hyperbolic"),
        ("cart_to_equinoctial i = π", c -> cart_to_equinoctial(c, μ), retro_eq, r"not defined for i ≈ π"),
        ("equinoctial_to_cart μ",  e -> equinoctial_to_cart(e, 0.0), [7000.0, 0.01, 0.0, 0.0, 0.0, 0.0], r"Gravitational parameter"),
        ("equinoctial_to_cart a",  e -> equinoctial_to_cart(e, μ), [-7000.0, 0.01, 0.0, 0.0, 0.0, 0.0], r"Semi-major axis"),
        ("equinoctial_to_cart e",  e -> equinoctial_to_cart(e, μ), [7000.0, 1.5, 0.0, 0.0, 0.0, 0.0], r"Eccentricity"),
        ("equinoctial_to_alt i = π", e -> equinoctial_to_alt_equinoctial(e), [7000.0, 0.0, 0.0, 1e13, 0.0, 0.0], r"i ≈ π"),
        ("alt_to_equinoctial i = π", a -> alt_equinoctial_to_equinoctial(a), [7000.0, 0.0, 0.0, 1.0, 0.0, 0.0], r"i ≈ π"),
        ("alt_to_equinoctial altp > 1", a -> alt_equinoctial_to_equinoctial(a), [7000.0, 0.0, 0.0, 1.2, 0.0, 0.0], r"exceeds 1"),
        ("kep_to_modkep parabolic", k -> kep_to_modkep(k), [7000.0, 1.0, 0.1, 0.0, 0.0, 0.0], r"parabolic"),
        ("kep_to_modkep e < 0",     k -> kep_to_modkep(k), [7000.0, -0.1, 0.1, 0.0, 0.0, 0.0], r"cannot be negative"),
        ("kep_to_modkep a > 0, e > 1", k -> kep_to_modkep(k), [7000.0, 1.5, 0.1, 0.0, 0.0, 0.0], r"cannot be positive"),
        ("kep_to_modkep a < 0, e < 1", k -> kep_to_modkep(k), [-7000.0, 0.5, 0.1, 0.0, 0.0, 0.0], r"cannot be negative if"),
        ("modkep_to_kep rₐ = 0",    m -> modkep_to_kep(m), [7000.0, 0.0, 0.1, 0.0, 0.0, 0.0], r"must not be zero"),
        ("modkep_to_kep rₐ < rₚ",   m -> modkep_to_kep(m), [8000.0, 7000.0, 0.1, 0.0, 0.0, 0.0], r"If rₐ < rₚ"),
        ("modkep_to_kep rₚ = 0",    m -> modkep_to_kep(m), [0.0, 8000.0, 0.1, 0.0, 0.0, 0.0], r"Singular conic"),
        ("modkep_to_kep hyperbola", m -> modkep_to_kep(m), [8000.0, -7000.0, 0.1, 0.0, 0.0, 0.0], r"For a hyperbola"),
        ("cart_to_sphazfpa r",      c -> cart_to_sphazfpa(c), zero_r, r"Position magnitude"),
        ("cart_to_sphazfpa v",      c -> cart_to_sphazfpa(c), zero_v, r"Velocity magnitude"),
        ("cart_to_sphradec r",      c -> cart_to_sphradec(c), zero_r, r"Radius is zero"),
        ("cart_to_sphradec v",      c -> cart_to_sphradec(c), zero_v, r"Velocity is zero"),
        # p = a(1 - e²) below tol on an ellipse: an orbit collapsed to a point
        ("kep_to_cart collapsed",   k -> kep_to_cart(k, μ), [1e-13, 0.5, 0.1, 0.0, 0.0, 0.0], r"parabolic or singular"),
        ("kep_to_modkep a = 0",     k -> kep_to_modkep(k), [0.0, 0.2, 0.3, 0.4, 0.5, 0.6], r"collapsed orbit"),
        ("modkep_to_kep rₐ = Inf",  m -> modkep_to_kep(m), [7000.0, Inf, 0.3, 0.4, 0.5, 0.6], r"infinite"),
        ("sphradec_to_cart r < 0",  s -> sphradec_to_cart(s), [-7000.0, 0.0, 0.0, 7.5, 0.0, 0.0], r"cannot be negative"),
        ("sphradec_to_cart v < 0",  s -> sphradec_to_cart(s), [7000.0, 0.0, 0.0, -7.5, 0.0, 0.0], r"cannot be negative"),
        ("sphazfpa_to_cart r < 0",  s -> sphazfpa_to_cart(s), [-7000.0, 0.0, 0.0, 7.5, 0.0, 0.0], r"cannot be negative"),
        ("sphazfpa_to_cart v < 0",  s -> sphazfpa_to_cart(s), [7000.0, 0.0, 0.0, -7.5, 0.0, 0.0], r"cannot be negative"),
    ]
    for (name, to_asym) in (("outasymptote", cart_to_outasymptote), ("inasymptote", cart_to_inasymptote))
        append!(cases, [
            ("cart_to_$name μ < 0",     c -> to_asym(c, -1.0), [1.0, 0.0, 0.0, 0.0, 1.0, 0.0], r"μ < tolerance"),
            ("cart_to_$name r",         c -> to_asym(c, μ), zero_r, r"degenerate position or velocity"),
            ("cart_to_$name h",         c -> to_asym(c, μ), radial, r"zero angular momentum"),
            ("cart_to_$name circular",  c -> to_asym(c, μ), circular, r"circular"),
            ("cart_to_$name parabolic", c -> to_asym(c, μ), parabola, r"parabolic"),
            ("cart_to_$name z-axis",    c -> to_asym(c, μ), name == "outasymptote" ? zaxis_out : zaxis_in,
             r"aligned with the z-axis"),
        ])
    end
    for (name, to_kep) in (("outasymptote", outasymptote_to_kep), ("inasymptote", inasymptote_to_kep))
        append!(cases, [
            ("$(name)_to_kep μ < 0",    a -> to_kep(a, -μ), [7000.0, 10.0, 0.3, 0.2, 0.5, 0.0], r"μ < tolerance"),
            ("$(name)_to_kep C₃ = 0",   a -> to_kep(a, μ), [7000.0, 0.0, 0.0, 0.0, 0.0, 0.0], r"nearly parabolic"),
            ("$(name)_to_kep rₚ ≤ 0",   a -> to_kep(a, μ), [-1.0, 10.0, 0.0, 0.0, 0.0, 0.0], r"must be positive"),
            ("$(name)_to_kep rₚ > a",   a -> to_kep(a, μ), [50000.0, -10.0, 0.0, 0.0, 0.0, 0.0], r"exceeds the semi-major axis"),
            ("$(name)_to_kep circular", a -> to_kep(a, μ), [7000.0, -μ / 7000.0, 0.0, 0.0, 0.0, 0.0], r"nearly circular"),
            ("$(name)_to_kep z-axis",   a -> to_kep(a, μ), [7000.0, 10.0, 0.0, π/2, 0.0, 0.0], r"aligned with z-axis"),
            # e = 1.18 puts the asymptote at 148.3°; ν = 3 rad is 171.9°, beyond it
            ("$(name)_to_kep beyond asymptote", a -> to_kep(a, μ), [7000.0, 10.0, 0.3, 0.2, 0.5, 3.0], r"beyond the asymptote"),
        ])
    end

    @testset "Unphysical conversions warn and return NaNs [AS-7]" begin
        @testset "$label" for (label, f, x, warning) in cases
            out = @test_logs (:warn, Regex("Conversion failed.*" * warning.pattern)) f(x)
            @test length(out) == 6 && all(isnan, out)
            # the NaNs keep the input's element type, so a dual number stays one [AS-9]
            duals = ForwardDiff.Dual.(x, 1.0)
            out_dual = @test_logs (:warn,) f(duals)
            @test eltype(out_dual) <: ForwardDiff.Dual
        end
    end

    vector_functions = [
        ("cart_to_kep", v -> cart_to_kep(v, μ)),         ("kep_to_cart", v -> kep_to_cart(v, μ)),
        ("cart_to_sphradec", cart_to_sphradec),           ("sphradec_to_cart", sphradec_to_cart),
        ("cart_to_sphazfpa", cart_to_sphazfpa),           ("sphazfpa_to_cart", sphazfpa_to_cart),
        ("cart_to_mee", v -> cart_to_mee(v, μ)),          ("mee_to_cart", v -> mee_to_cart(v, μ)),
        ("cart_to_equinoctial", v -> cart_to_equinoctial(v, μ)),
        ("equinoctial_to_cart", v -> equinoctial_to_cart(v, μ)),
        ("equinoctial_to_alt_equinoctial", equinoctial_to_alt_equinoctial),
        ("alt_equinoctial_to_equinoctial", alt_equinoctial_to_equinoctial),
        ("kep_to_modkep", kep_to_modkep),                 ("modkep_to_kep", modkep_to_kep),
        ("cart_to_outasymptote", v -> cart_to_outasymptote(v, μ)),
        ("outasymptote_to_kep", v -> outasymptote_to_kep(v, μ)),
        ("cart_to_inasymptote", v -> cart_to_inasymptote(v, μ)),
        ("inasymptote_to_kep", v -> inasymptote_to_kep(v, μ)),
        ("kep_to_brouwer_mean_short", v -> kep_to_brouwer_mean_short(v, μ)),
        ("brouwer_mean_short_to_kep", v -> brouwer_mean_short_to_kep(v, μ)),
        ("kep_to_brouwer_mean_long", v -> kep_to_brouwer_mean_long(v, μ)),
        ("brouwer_mean_long_to_kep", v -> brouwer_mean_long_to_kep(v, μ)),
    ]

    @testset "Modified equinoctial warns near its singularity and still converts [AS-4]" begin
        near_retro = kep_to_cart([8000.0, 0.01, π - 1e-5, 0.3, 0.4, 1.0], μ)
        mee = @test_logs (:warn, r"near their singularity") cart_to_mee(near_retro, μ)
        @test all(isfinite, mee)
        @test_logs cart_to_mee(near_retro, μ; j = -1)          # the retrograde set is fine there
        @test_logs cart_to_mee(kep_to_cart([8000.0, 0.01, π - 1e-3, 0.3, 0.4, 1.0], μ), μ)
    end

    @testset "Each failure warns once, with its own cause [AS-7]" begin
        # NaNs passed on are not a second failure, in any vector function, nor is a NaN μ
        @testset "$name" for (name, f) in vector_functions
            out = @test_logs f(fill(NaN, 6))
            @test length(out) == 6 && all(isnan, out)
        end
        for f in (cart_to_kep, cart_to_outasymptote, brouwer_mean_long_to_kep, brouwer_mean_short_to_kep)
            out = @test_logs f([9000.0, 0.05, 0.9, 0.3, 0.4, 1.0], NaN)
            @test all(isnan, out)
        end
        # a Keplerian state that fails reports its own cause, not a second one in the next step
        @test_logs (:warn, r"cannot be negative") SphericalRADECState(
            KeplerianState(7000.0, -0.1, 0.3, 0.4, 0.5, 0.6), μ)
        # a chain reports the first cause only
        beyond = OutGoingAsymptoteState(7000.0, 10.0, 0.3, 0.2, 0.5, 3.0)
        @test_logs (:warn, r"beyond the asymptote") CartesianState(beyond, μ)
        # cart_to_kep's own parabola is reported as one by kep_to_cart
        kep = cart_to_kep(parabola, μ)
        @test_logs (:warn, r"parabolic") kep_to_cart(kep, μ)
        # the near-singularity advice names the other set
        near_pro = kep_to_cart([8000.0, 0.01, 1e-5, 0.3, 0.4, 1.0], μ)
        @test_logs (:warn, r"j = 1 suits") cart_to_mee(near_pro, μ; j = -1)
    end

    @testset "A Float32 retrograde equatorial orbit round-trips [AS-4, AS-5.3]" begin
        μ32 = Float32(μ)
        c32 = Float32[7000, 0, 0, 0, -7.5, 0]
        kep = cart_to_kep(c32, μ32)
        @test kep[3] == Float32(π)                       # above π, in Float32
        back = @test_logs kep_to_cart(kep, μ32)
        @test eltype(back) == Float32
        @test isapprox(back, c32; rtol = 1e-5)
        @test all(isfinite, to_vector(ModifiedKeplerianState(CartesianState(c32), μ32)))
    end

    @testset "Alternate equinoctial refuses near i = π where the equinoctial set does [AS-4, AS-7]" begin
        eq(Δ) = [8000.0, 0.01, 0.02, tan((π - Δ) / 2) * cos(0.3), tan((π - Δ) / 2) * sin(0.3), 1.0]
        for Δ in (1e-3, 2e-4)                            # clear of the threshold: a round trip
            e = eq(Δ)
            back = @test_logs alt_equinoctial_to_equinoctial(equinoctial_to_alt_equinoctial(e))
            @test isapprox(back, e; rtol = 1e-8)
        end
        for Δ in (1e-5, 1e-7, 1e-9, 1e-11)               # inside it: refused, never a silent error
            alt = [8000.0, 0.01, 0.02, cos(Δ / 2) * cos(0.3), cos(Δ / 2) * sin(0.3), 1.0]
            out = @test_logs (:warn, r"i ≈ π") alt_equinoctial_to_equinoctial(alt)
            @test all(isnan, out)
        end
        # Float32 rounding just over sin²(i/2) = 1 is an orbit at i ≈ π, not an impossible one
        alt32 = Float32[8000, 0.01, 0.02, cos(5f-5) * cos(0.3f0), cos(5f-5) * sin(0.3f0), 1]
        @test_logs (:warn, r"i ≈ π") alt_equinoctial_to_equinoctial(alt32)
    end

    @testset "Integer input converts as its floating-point value [AS-2, AS-13]" begin
        int_inputs = [
            (v -> cart_to_kep(v, μ), [7000, 0, 0, 0, 8, 1]),         (v -> kep_to_cart(v, μ), [7000, 0, 1, 0, 0, 1]),
            (cart_to_sphradec, [7000, 0, 0, 0, 8, 1]),                (sphradec_to_cart, [7000, 0, 0, 8, 1, 0]),
            (cart_to_sphazfpa, [7000, 0, 0, 0, 8, 1]),                (sphazfpa_to_cart, [7000, 0, 0, 8, 1, 1]),
            (v -> cart_to_mee(v, μ), [7000, 0, 0, 0, 8, 1]),          (v -> mee_to_cart(v, μ), [7000, 0, 0, 0, 0, 1]),
            (v -> cart_to_equinoctial(v, μ), [7000, 0, 0, 0, 8, 1]),  (v -> equinoctial_to_cart(v, μ), [7000, 0, 0, 0, 0, 1]),
            (equinoctial_to_alt_equinoctial, [7000, 0, 0, 0, 0, 1]),  (alt_equinoctial_to_equinoctial, [7000, 0, 0, 0, 0, 1]),
            (kep_to_modkep, [7000, 0, 1, 0, 0, 1]),                   (modkep_to_kep, [7000, 8000, 1, 0, 0, 1]),
            (v -> cart_to_outasymptote(v, μ), [7000, 0, 0, 0, 12, 1]), (v -> outasymptote_to_kep(v, μ), [7000, 10, 0, 0, 0, 0]),
            (v -> cart_to_inasymptote(v, μ), [7000, 0, 0, 0, 12, 1]),  (v -> inasymptote_to_kep(v, μ), [7000, 10, 0, 0, 0, 0]),
            (v -> kep_to_brouwer_mean_short(v, μ), [9000, 0, 1, 0, 0, 1]),
            (v -> brouwer_mean_short_to_kep(v, μ), [9000, 0, 1, 0, 0, 1]),
            (v -> kep_to_brouwer_mean_long(v, μ), [9000, 0, 1, 0, 0, 1]),
            (v -> brouwer_mean_long_to_kep(v, μ), [9000, 0, 1, 0, 0, 1]),
        ]
        for (f, x) in int_inputs
            out = f(x)
            @test eltype(out) <: AbstractFloat
            @test isequal(out, f(float.(x)))
        end
    end

    @testset "Angles come back in [0, 2π) [AS-11]" begin
        c = kep_to_cart([8000.0, 0.1, 0.5, 0.0, 0.0, -1e-17], μ)
        @test 0 <= cart_to_equinoctial(c, μ)[6] < 2π
        @test 0 <= cart_to_mee(c, μ)[6] < 2π
        @test 0 <= cart_to_kep(c, μ)[6] < 2π
    end

    @testset "Malformed input names what it received [AS-8]" begin
        @testset "$name" for (name, f) in vector_functions, n in (5, 7)
            err = try f(collect(1.0:n)); nothing catch e; e end
            @test err isa Exception
            @test occursin("got $n elements", sprint(showerror, err))
        end
    end
end
