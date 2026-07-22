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

# At the prograde critical inclination (~63.4°) GMAT's bisubc≥0.001 trigger (bisubc ∝ signed
# θ⁵) zeroes the long-period terms → the LONG theory is degraded there (GMAT warns). SHORT has
# no such branch and stays accurate. i≈117° does NOT trigger (θ⁵<0), so it converges. We
# faithfully reproduce this GMAT behavior.
@testset "BROU critical-inclination degradation (i≈63.4°)" begin
    osc = [9000.0, 0.05, deg2rad(63.4), deg2rad(15.0), deg2rad(70.0), deg2rad(100.0)]
    @test _rt_err(brouwer_mean_short_to_kep(kep_to_brouwer_mean_short(osc, _BROU_MU), _BROU_MU), osc) < 1e-5
    ml = kep_to_brouwer_mean_long(osc, _BROU_MU)
    @test all(isfinite, ml) && length(ml) == 6
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
            @test T(b, _BROU_MU) isa T
        end
        for s in ell_sources                   # TO Brouwer FROM every elliptic representation
            @test BT(s, _BROU_MU) isa BT
        end
    end
    # the originally-missing case, explicitly
    @test ModifiedEquinoctialState(BrouwerMeanLongState(kep, _BROU_MU), _BROU_MU) isa ModifiedEquinoctialState
end

@testset "BROU loud failure / domain [FR-BROU-5]" begin
    # public osc->mean entry checks
    @test_throws ErrorException kep_to_brouwer_mean_long([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_short([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_long([7000.0, -0.01, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_short([3200.0, 0.5, 0.3, 0.0, 0.0, 0.0], _BROU_MU)  # radper 1600
    @test_throws ErrorException kep_to_brouwer_mean_long([3200.0, 0.5, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_short([7000.0, 0.01, deg2rad(200.0), 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_long([7000.0, 0.01, deg2rad(200.0), 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException kep_to_brouwer_mean_long([7000.0, 0.01, 0.3, 0.0, 0.0, 0.0], 42000.0)   # non-Earth μ
    @test_throws ErrorException kep_to_brouwer_mean_short([7000.0, 0.01, 0.3, 0.0, 0.0, 0.0], 42000.0)
    # direct mean->osc core checks
    @test_throws ErrorException brouwer_mean_long_to_kep([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_short_to_kep([7000.0, 0.995, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_short_to_kep([3200.0, 0.5, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_long_to_kep([3200.0, 0.5, 0.3, 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_short_to_kep([7000.0, 0.01, deg2rad(200.0), 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_long_to_kep([7000.0, 0.01, deg2rad(200.0), 0.0, 0.0, 0.0], _BROU_MU)
    @test_throws ErrorException brouwer_mean_long_to_kep([7000.0, 0.01, 0.3, 0.0, 0.0, 0.0], 42000.0)
    @test_throws ErrorException brouwer_mean_short_to_kep([7000.0, 0.01, 0.3, 0.0, 0.0, 0.0], 42000.0)
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
