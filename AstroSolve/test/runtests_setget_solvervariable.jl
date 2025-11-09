using Test


using AstroSolve
using AstroFun
using AstroStates
using AstroEpochs
using AstroManeuvers

@testset "SolverVariable set/get (ManeuverCalc DeltaVVector)" begin
    sc = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
    )

    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
    calc = ManeuverCalc(man, sc, DeltaVVector())
    sv   = SolverVariable(calc = calc, name = "dv-full",
                          lower_bound = [-10.0, 0.0, 0.0], upper_bound = [10.0, 0.0, 0.0])

    # Initial get via AstroSolve API
    @test get_sol_var(sv) == [0.1, 0.2, 0.3]

    # Set via AstroSolve API (delegates to set_calc!)
    set_sol_var(sv, [0.2, 0.3, 0.4])
    @test get_sol_var(sv) == [0.2, 0.3, 0.4]
    @test man.element1 == 0.2 && man.element2 == 0.3 && man.element3 == 0.4

    # Length mismatch should error
    @test_throws ArgumentError set_sol_var(sv, [1.0, 2.0])
end

@testset "SolverVariable set/get (OrbitCalc Ecc)" begin
    sc = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
    )

    calc = OrbitCalc(sc, Ecc())
    sv   = SolverVariable(calc = calc, name = "ecc", lower_bound = -Inf, upper_bound = Inf)

    # Initial get matches current eccentricity
    sc_kep = get_state(sc, Keplerian())
    @test isapprox(get_sol_var(sv)[1], sc_kep.ecc; atol=1e-8)

    # Set via AstroSolve API; verify spacecraft and get() reflect change
    set_sol_var(sv, [0.02])
    sc_kep = get_state(sc, Keplerian())
    @test isapprox(sc_kep.ecc, 0.02; atol=1e-8)
    @test isapprox(get_sol_var(sv)[1], 0.02; atol=1e-8)
end

@testset "SolverVariable set/get (BodyCalc GravParam)" begin
    # Local mutable body with ASCII mu field to avoid global side effects
    b = CelestialBody(
        name = "TestBody",
        mu   = 3.986004418e14,
        equatorial_radius = 6378.137e3,
        flattening = 1/298.257223563,
        naifid = 399,
    )

    calc = BodyCalc(b, GravParam())
    sv   = SolverVariable(calc = calc, name = "mu", lower_bound = -Inf, upper_bound = Inf)

    # Initial get
    @test isapprox(get_sol_var(sv)[1], 3.986004418e14; atol=1e-8)

    # Set via AstroSolve API; verify body and get() reflect change
    set_sol_var(sv, [4.0e14])
    @test b.mu == 4.0e14
    @test get_sol_var(sv)[1] == 4.0e14
end

@testset "SolverVariable set/get (OrbitCalc IncomingAsymptoteFull)" begin
    ia0 = IncomingAsymptoteState(7000.0, 0.01, 1.0, 0.5, 0.2, 0.0)  # [rp,c3,rla,dla,bpa,ta]
    sc  = Spacecraft(state = ia0, time = Time("2020-09-21T12:23:12", TAI(), ISOT()))

    calc = OrbitCalc(sc, IncomingAsymptote())
    sv   = SolverVariable(calc = calc, name = "ia-full",
                          lower_bound = fill(-Inf, 6), upper_bound = fill(Inf, 6))

    # Initial get matches state
    vals0 = get_sol_var(sv)
    @test length(vals0) == 6
    @test all(isapprox.(vals0, to_vector(ia0); atol=1e-10))

    # Set via AstroSolve API; verify spacecraft and get() reflect change
    ia_new = [7100.0, 0.02, 1.1, 0.6, 0.25, 0.05]
    set_sol_var(sv, ia_new)

    vals = get_sol_var(sv)
    @test length(vals) == 6
    @test all(isapprox.(vals, ia_new; atol=1e-10))

    sc_ia = get_state(sc, IncomingAsymptote())
    @test all(isapprox.(to_vector(sc_ia), ia_new; atol=1e-10))
end

nothing