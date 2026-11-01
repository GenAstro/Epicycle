# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# BROU — Brouwer-Lyddane mean element conversions (FR-BROU-1..5, PR-BROU-1, SR-AD-1)
# Reimplemented from GMAT StateConversionUtil.cpp. Classic @testset style.
using Test
using AstroStates
using ForwardDiff

const _BROU_MU = 398600.4415
_wrappi(x) = mod(x + pi, 2pi) - pi
_angd(a, b) = abs(_wrappi(a - b))
_rt_err(back, osc) = max(abs(back[1]-osc[1])/osc[1], abs(back[2]-osc[2]),
                         _angd(back[3],osc[3]), _angd(back[4],osc[4]),
                         _angd(back[5],osc[5]), _angd(back[6],osc[6]))

# osculating Keplerian test orbits [a(km), e, i, Ω, ω, ν] (radians), inside validated domain
_brou_cases = [
    ("LEO",       [7000.0,  0.01, deg2rad(28.5), deg2rad(80.0),  deg2rad(40.0),  deg2rad(120.0)]),
    ("MEO",       [12000.0, 0.05, deg2rad(55.0), deg2rad(200.0), deg2rad(10.0),  deg2rad(300.0)]),
    ("GTOish",    [24000.0, 0.30, deg2rad(20.0), deg2rad(10.0),  deg2rad(90.0),  deg2rad(45.0)]),
    ("nearpolar", [8000.0,  0.02, deg2rad(98.0), deg2rad(120.0), deg2rad(200.0), deg2rad(260.0)]),
]

@testset "BROU round-trip (osc→mean→osc) [FR-BROU-1..4, PR-BROU-1]" begin
    for (nm, osc) in _brou_cases
        for (o2m, m2o) in ((kep_to_brouwer_mean_long, brouwer_mean_long_to_kep),
                           (kep_to_brouwer_mean_short, brouwer_mean_short_to_kep))
            back = m2o(o2m(osc, _BROU_MU), _BROU_MU)
            @test _rt_err(back, osc) < 1e-6
        end
    end
end

# Extreme regimes: high eccentricity, retrograde pseudostate (i>175°), critical
# inclination (~63°/117°), small inclination. Round-trip is self-consistent, so it
# closes even where the theory is degraded; tolerances reflect that.
@testset "BROU round-trip — extreme regimes" begin
    extreme = [
        ("ecc0.6",   [10000.0, 0.60, deg2rad(45.0),  deg2rad(30.0),  deg2rad(50.0),  deg2rad(200.0)], 1e-6),
        ("ecc0.7",   [15000.0, 0.70, deg2rad(30.0),  deg2rad(10.0),  deg2rad(80.0),  deg2rad(150.0)], 1e-6),
        ("retro177", [8000.0,  0.01, deg2rad(177.0), deg2rad(120.0), deg2rad(200.0), deg2rad(260.0)], 1e-5),
        ("crit117",  [9000.0,  0.05, deg2rad(117.0), deg2rad(15.0),  deg2rad(70.0),  deg2rad(100.0)], 1e-5),
        ("lowinc",   [7500.0,  0.02, deg2rad(0.5),   deg2rad(60.0),  deg2rad(30.0),  deg2rad(220.0)], 1e-5),
    ]
    for (nm, osc, tol) in extreme
        for (o2m, m2o) in ((kep_to_brouwer_mean_long, brouwer_mean_long_to_kep),
                           (kep_to_brouwer_mean_short, brouwer_mean_short_to_kep))
            back = m2o(o2m(osc, _BROU_MU), _BROU_MU)
            @test _rt_err(back, osc) < tol
        end
    end
end

@testset "BROU graceful degradation (near-critical + high ecc)" begin
    bad = [15000.0, 0.70, deg2rad(63.4), deg2rad(10.0), deg2rad(80.0), deg2rad(150.0)]
    ml = kep_to_brouwer_mean_long(bad, _BROU_MU)   # warns + returns best estimate, must not throw
    @test all(isfinite, ml) && length(ml) == 6
end

# Near the critical inclinations (~63.4° and ~116.6°) |bisubc| ≥ 0.001 zeroes the long-period
# terms, so the LONG theory is degraded there and warns. SHORT has no such branch and stays
# accurate. bisubc carries the sign of θ⁵, so a test on its signed value would miss 116.6°.
@testset "BROU critical-inclination degradation (i≈63.4°)" begin
    osc = [9000.0, 0.05, deg2rad(63.4), deg2rad(15.0), deg2rad(70.0), deg2rad(100.0)]
    @test _rt_err(brouwer_mean_short_to_kep(kep_to_brouwer_mean_short(osc, _BROU_MU), _BROU_MU), osc) < 1e-5
    ml = kep_to_brouwer_mean_long(osc, _BROU_MU)
    @test all(isfinite, ml) && length(ml) == 6
end

# Where the long-period corrections are dropped, the orbit must still keep its orientation. With
# them dropped there is nothing in the theory that depends on Ω, so turning the node turns the
# orbit about z: the returned node moves with it and every other element stays put.
@testset "BROU critical inclination keeps the orientation [AS-10]" begin
    mean(Ω) = [9000.0, 0.05, deg2rad(63.4), Ω, 0.4, 1.0]
    o1 = @test_logs (:warn, r"near critical") brouwer_mean_long_to_kep(mean(0.3), _BROU_MU)
    o2 = @test_logs (:warn, r"near critical") brouwer_mean_long_to_kep(mean(0.8), _BROU_MU)
    @test _angd(o2[4] - o1[4], 0.5) < 1e-12
    @test isapprox(o1[[1, 2, 3, 5, 6]], o2[[1, 2, 3, 5, 6]]; rtol = 1e-12)
    @test _angd(o1[5], 0.4) < 1e-3                   # ω stays near its mean value, not near -Ω
end

@testset "BROU osculating→mean iteration warns on every exit but convergence [AS-7, AS-10]" begin
    # The residual stops falling: near-circular at the critical inclination
    stall = [7000.0, 0.0, deg2rad(63.4), 0.17, 0.0, 0.0]
    m = @test_logs (:warn, r"did not converge \(relative position residual") match_mode=:any kep_to_brouwer_mean_long(stall, _BROU_MU)
    @test all(isfinite, m)
    # An estimate leaves the theory's domain
    leave = [10000.0, 0.7, deg2rad(62.5), 0.17, 0.0, 0.0]
    m = @test_logs (:warn, r"did not converge \(near-singular") match_mode=:any kep_to_brouwer_mean_long(leave, _BROU_MU)
    @test all(isfinite, m)
    # A converged case is silent
    @test_logs kep_to_brouwer_mean_long([8000.0, 0.01, deg2rad(30.0), 0.17, 0.4, 1.0], _BROU_MU)
end

@testset "BROU retrograde critical inclination warns and converts [AS-7, AS-10]" begin
    # Unsuppressed, the long-period terms here put sin(i/2) above 1 and asin threw a DomainError
    for i in (116.4, 116.5, 116.6, 116.7), e in (0.1, 0.3)
        mean = [7000.0, e, deg2rad(i), 0.17, 0.0, 0.0]
        osc = @test_logs (:warn, r"near critical") match_mode=:any brouwer_mean_long_to_kep(mean, _BROU_MU)
        @test all(isfinite, osc)
        @test abs(osc[3] - mean[3]) < 1e-3
    end
end

@testset "BROU domain boundaries are excluded [AS-7, AS-10]" begin
    # e = 0.99 and a periapsis radius of exactly 3000 km are outside, in both directions
    for f in (brouwer_mean_long_to_kep, brouwer_mean_short_to_kep, kep_to_brouwer_mean_long, kep_to_brouwer_mean_short)
        out = @test_logs (:warn, r"ECC") match_mode=:any f([400000.0, 0.99, 0.3, 0.1, 0.2, 0.3], _BROU_MU)
        @test all(isnan, out)
        out = @test_logs (:warn, r"RadPer must be > 3000") match_mode=:any f([6000.0, 0.5, 0.3, 0.1, 0.2, 0.3], _BROU_MU)
        @test all(isnan, out)
    end
end

@testset "BROU negative mean eccentricity is checked after it is normalized [AS-7]" begin
    # e = -0.8 is e = 0.8 with periapsis turned round, a periapsis radius of 1400 km
    out = @test_logs (:warn, r"RadPer must be > 3000") brouwer_mean_short_to_kep([7000.0, -0.8, 0.3, 0.1, 0.2, 0.3], _BROU_MU)
    @test all(isnan, out)
    # the long-period theory takes a negative e the same way, and |e| > 1 is refused, not thrown
    out = @test_logs (:warn, r"ECC must be < 0.99 \(got -2.0\)") brouwer_mean_long_to_kep([7000.0, -2.0, 0.3, 0.1, 0.2, 0.3], _BROU_MU)
    @test all(isnan, out)
    mean    = [9000.0, -0.01, deg2rad(30.0), 0.2, 0.4, 1.0]
    flipped = [9000.0,  0.01, deg2rad(30.0), 0.2, 0.4 + π, 1.0 - π]
    @test _rt_err(brouwer_mean_long_to_kep(mean, _BROU_MU), brouwer_mean_long_to_kep(flipped, _BROU_MU)) < 1e-12
end

@testset "BROU struct + OrbitState wiring [FR-BROU-1..4]" begin
    kep = KeplerianState(7000.0, 0.01, deg2rad(28.5), deg2rad(80.0), deg2rad(40.0), deg2rad(120.0))
    for BT in (BrouwerMeanLongState, BrouwerMeanShortState)
        b = BT(kep, _BROU_MU)
        @test b isa BT
        @test BT(b, _BROU_MU) === b                       # identity
        @test length(to_vector(b)) == 6
        kb = KeplerianState(b, _BROU_MU)
        @test maximum(abs.(to_vector(kb) .- to_vector(kep))) < 1e-6
        os = OrbitState(b)
        @test to_state(os) isa BT
        @test to_vector(os) == to_vector(b)
        io = IOBuffer(); show(io, b); @test occursin("Brouwer", String(take!(io)))
        cart = CartesianState(kep, _BROU_MU)
        b2 = BT(cart, _BROU_MU)
        @test maximum(abs.(to_vector(CartesianState(b2, _BROU_MU)) .- to_vector(cart))) < 1e-6
        v = to_vector(b); @test BT(v) isa BT              # no-μ vector ctor
    end
    # cross conversion Long <-> Short via Keplerian pivot
    bl = BrouwerMeanLongState(kep, _BROU_MU)
    bs = BrouwerMeanShortState(bl, _BROU_MU)
    @test bs isa BrouwerMeanShortState
end

@testset "BROU full permutation coverage (all state types ↔ Brouwer)" begin
    kep = KeplerianState(7000.0, 0.02, deg2rad(45.0), deg2rad(80.0), deg2rad(40.0), deg2rad(120.0))
    ell_sources = [CartesianState(kep, _BROU_MU), kep,
                   SphericalRADECState(kep, _BROU_MU), SphericalAZIFPAState(kep, _BROU_MU),
                   ModifiedEquinoctialState(kep, _BROU_MU), ModifiedKeplerianState(kep, _BROU_MU),
                   EquinoctialState(kep, _BROU_MU), AlternateEquinoctialState(kep, _BROU_MU)]
    targets = [CartesianState, KeplerianState, SphericalRADECState, SphericalAZIFPAState,
               ModifiedEquinoctialState, ModifiedKeplerianState, EquinoctialState,
               AlternateEquinoctialState, OutGoingAsymptoteState, IncomingAsymptoteState,
               BrouwerMeanLongState, BrouwerMeanShortState]
    for BT in (BrouwerMeanLongState, BrouwerMeanShortState)
        b = BT(kep, _BROU_MU)
        for T in targets                       # FROM Brouwer TO every representation
            out = T(b, _BROU_MU)
            @test out isa T
            @test all(isfinite, to_vector(out))
        end
        for s in ell_sources                   # TO Brouwer FROM every elliptic representation
            out = BT(s, _BROU_MU)
            @test out isa BT
            @test isapprox(to_vector(out), to_vector(b); rtol = 1e-8)
        end
    end
    # the originally-missing case, explicitly
    @test ModifiedEquinoctialState(BrouwerMeanLongState(kep, _BROU_MU), _BROU_MU) isa ModifiedEquinoctialState
end

@testset "BROU outside the domain warns and returns NaNs [AS-7, AS-10]" begin
    bad_mean_and_osc = [
        ([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], _BROU_MU, "(ECC|RadPer)"),         # e ≥ 0.99; the short core checks periapsis first
        ([3200.0, 0.5, 0.3, 0.0, 0.0, 0.0], _BROU_MU, "RadPer"),                 # periapsis 1600 km
        ([7000.0, 0.01, deg2rad(200.0), 0.0, 0.0, 0.0], _BROU_MU, "INC"),        # i > 180°
        ([7000.0, 0.01, 0.3, 0.0, 0.0, 0.0], 42000.0, "only to the Earth"),      # non-Earth μ
    ]
    for (elements, μ, reason) in bad_mean_and_osc
        for f in (kep_to_brouwer_mean_long, kep_to_brouwer_mean_short,
                  brouwer_mean_long_to_kep, brouwer_mean_short_to_kep)
            out = @test_logs (:warn, Regex("Conversion failed.*$reason")) f(elements, μ)
            @test length(out) == 6 && all(isnan, out)
        end
    end
    # a negative osculating eccentricity is checked on the osculating side only
    for f in (kep_to_brouwer_mean_long, kep_to_brouwer_mean_short)
        out = @test_logs (:warn, r"Conversion failed.*ECC") f([7000.0, -0.01, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
        @test all(isnan, out)
    end
    # the state types carry the NaNs through
    kep = KeplerianState(7000.0, 0.995, 0.3, 0.0, 0.0, 0.0)
    for BT in (BrouwerMeanLongState, BrouwerMeanShortState)
        b = @test_logs (:warn, r"Conversion failed") BT(kep, _BROU_MU)
        @test all(isnan, to_vector(b))
    end
    # AD: the NaNs keep the input's element type, so a dual number stays one
    duals = ForwardDiff.Dual.([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], 1.0)
    out = @test_logs (:warn,) kep_to_brouwer_mean_long(duals, _BROU_MU)
    @test eltype(out) <: ForwardDiff.Dual && all(isnan ∘ ForwardDiff.value, out)
end

@testset "BROU malformed input throws [AS-8]" begin
    for f in (kep_to_brouwer_mean_long, kep_to_brouwer_mean_short,
              brouwer_mean_long_to_kep, brouwer_mean_short_to_kep)
        err = try f([7000.0, 0.01, 0.3, 0.0, 0.0], _BROU_MU); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("got 5 elements", err.msg)
    end
end

@testset "BROU warnings (inside-Earth radper, critical inc)" begin
    # radper in (3000, 6378) -> inside-Earth warning path; still returns finite
    @test all(isfinite, kep_to_brouwer_mean_short([6500.0, 0.02, deg2rad(45.0), 0.1, 0.2, 0.3], _BROU_MU))
    @test all(isfinite, kep_to_brouwer_mean_long([6500.0, 0.02, deg2rad(45.0), 0.1, 0.2, 0.3], _BROU_MU))
    @test all(isfinite, brouwer_mean_short_to_kep([6500.0, 0.02, deg2rad(45.0), 0.1, 0.2, 0.3], _BROU_MU))
    @test all(isfinite, brouwer_mean_long_to_kep([6500.0, 0.02, deg2rad(45.0), 0.1, 0.2, 0.3], _BROU_MU))
    # critical-inclination warning path (public long + core long bisubc branch)
    @test all(isfinite, kep_to_brouwer_mean_long([9000.0, 0.05, deg2rad(63.4), 0.1, 0.2, 0.3], _BROU_MU))
    @test all(isfinite, brouwer_mean_long_to_kep([9000.0, 0.05, deg2rad(63.4), 0.1, 0.2, 0.3], _BROU_MU))
end

@testset "BROU degenerate branches (circular / equatorial)" begin
    # near-equatorial mean -> osc exercises inc<=1e-7 / raan=0 branches (long core)
    @test all(isfinite, brouwer_mean_long_to_kep([7000.0, 0.01, 0.0,        deg2rad(0.0), deg2rad(30.0), deg2rad(60.0)], _BROU_MU))
    @test all(isfinite, brouwer_mean_short_to_kep([7000.0, 0.01, 0.0,       deg2rad(0.0), deg2rad(30.0), deg2rad(60.0)], _BROU_MU))
    # near-circular
    @test all(isfinite, brouwer_mean_long_to_kep([7000.0, 1e-9, deg2rad(45.0), 0.1, 0.2, deg2rad(60.0)], _BROU_MU))
    @test all(isfinite, brouwer_mean_short_to_kep([7000.0, 1e-9, deg2rad(45.0), 0.1, 0.2, deg2rad(60.0)], _BROU_MU))
    # tiny-ecc AND zero-inc -> osculating output is near-circular-equatorial (ecc>1e-11 & inc<=1e-7 branch)
    @test all(isfinite, brouwer_mean_long_to_kep([7000.0, 1e-9, 0.0, deg2rad(30.0), deg2rad(40.0), deg2rad(60.0)], _BROU_MU))
end

@testset "BROU differentiability [SR-AD-1]" begin
    fd_jac(f, x; h=1e-6) = hcat([(f(x .+ h .* ((1:length(x)).==k)) .- f(x .- h .* ((1:length(x)).==k)))./(2h) for k in 1:length(x)]...)
    meanv = [7100.0, 0.02, deg2rad(50.0), deg2rad(80.0), deg2rad(40.0), deg2rad(120.0)]
    oscv  = [7100.0, 0.02, deg2rad(50.0), deg2rad(80.0), deg2rad(40.0), deg2rad(120.0)]
    # direct mean->osc (closed form)
    for f in (m -> brouwer_mean_long_to_kep(m, _BROU_MU), m -> brouwer_mean_short_to_kep(m, _BROU_MU))
        J = ForwardDiff.jacobian(f, meanv)
        @test all(isfinite, J)
        @test maximum(abs.(J .- fd_jac(f, meanv))) < 1e-4
    end
    # iterative osc->mean: AD must flow through the fixed point (O-5)
    for f in (o -> kep_to_brouwer_mean_long(o, _BROU_MU), o -> kep_to_brouwer_mean_short(o, _BROU_MU))
        J = ForwardDiff.jacobian(f, oscv)
        @test all(isfinite, J)
        @test maximum(abs.(J .- fd_jac(f, oscv))) < 1e-3
    end
end

# ---------------------------------------------------------------------------
# TIER-2 INDEPENDENT TRUTH (GMAT). osc/long/short in DEGREES from GMAT DefaultSC.
# ---------------------------------------------------------------------------
_gmat_truth = NamedTuple[]
push!(_gmat_truth, (name="retrograde",
    osc  =[7191.938817629009, 0.002454974900598474, 90.85008005658096, 306.6148021947984, 314.1905515359845, 99.88774933205636],
    long =[7194.814398118431, 0.003503464010814326, 90.8499133356497,  306.6153120402385, 318.1042609106366, 95.64492369924062],
    short=[7194.807565972313, 0.002876142756243898, 90.84991202281407, 306.6153090500403, 335.0214533896338, 78.7278355299037]))
push!(_gmat_truth, (name="prograde",
    osc  =[7700.000000000001, 0.009999999999999952, 28.5, 45.0, 120.0000000000013, 44.99999999999866],
    long =[7698.210899541718, 0.008949857021775129, 28.48871302547095, 45.01587422472831, 117.3804245507798, 46.80430093128462],
    short=[7698.20743211596,  0.009333025431126935, 28.4883419997415,  45.01545072752478, 116.163027340027,  48.02187357619155]))

@testset "BROU vs GMAT truth [FR-BROU-1..4] (tier-2)" begin
    d2r(v) = [v[1], v[2], deg2rad(v[3]), deg2rad(v[4]), deg2rad(v[5]), deg2rad(v[6])]
    resid(got, ref) = [abs(got[1]-ref[1]), abs(got[2]-ref[2]),
                       rad2deg(_angd(got[3],ref[3])), rad2deg(_angd(got[4],ref[4])),
                       rad2deg(_angd(got[5],ref[5])), rad2deg(_angd(got[6],ref[6]))]
    for row in _gmat_truth
        osc = d2r(row.osc)
        ml = kep_to_brouwer_mean_long(osc, _BROU_MU);  rl = resid(ml, d2r(row.long))
        ms = kep_to_brouwer_mean_short(osc, _BROU_MU); rs = resid(ms, d2r(row.short))
        println("  [$(row.name)] LONG  dSMA(km)=$(rl[1])  dECC=$(rl[2])  dANG(deg)=$(maximum(rl[3:6]))")
        println("  [$(row.name)] SHORT dSMA(km)=$(rs[1])  dECC=$(rs[2])  dANG(deg)=$(maximum(rs[3:6]))")
        @test rl[1] < 1e-3 && rl[2] < 1e-8 && maximum(rl[3:6]) < 1e-5
        @test rs[1] < 1e-3 && rs[2] < 1e-8 && maximum(rs[3:6]) < 1e-5
    end
end

@testset "BROU state types against GMAT truth [AS-1.11, AS-1.12, AS-10]" begin
    d2r(v) = [v[1], v[2], deg2rad(v[3]), deg2rad(v[4]), deg2rad(v[5]), deg2rad(v[6])]
    for row in _gmat_truth
        kep = KeplerianState(d2r(row.osc))
        for (BT, ref) in ((BrouwerMeanLongState, row.long), (BrouwerMeanShortState, row.short))
            b = to_vector(BT(kep, _BROU_MU))
            ref = d2r(ref)
            @test abs(b[1] - ref[1]) < 1e-3 && abs(b[2] - ref[2]) < 1e-8
            @test maximum(_angd.(b[3:6], ref[3:6])) < deg2rad(1e-5)
        end
    end
end

# GMAT's short-period theory takes a negative mean eccentricity as the same orbit with
# periapsis turned half way round: e → -e, ω → ω + π, M → M - π.
@testset "BROU negative mean eccentricity" begin
    mean    = [7000.0, -0.01, deg2rad(30.0), 0.2, 0.4, 1.0]
    flipped = [7000.0,  0.01, deg2rad(30.0), 0.2, 0.4 + π, 1.0 - π]
    @test _rt_err(brouwer_mean_short_to_kep(mean, _BROU_MU),
                  brouwer_mean_short_to_kep(flipped, _BROU_MU)) < 1e-12
end
