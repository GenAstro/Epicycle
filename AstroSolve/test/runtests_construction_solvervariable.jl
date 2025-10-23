using Test


using AstroSolve
using AstroFun
using AstroStates
using AstroEpochs
using AstroMan

@testset "SolverVariable calc-constructor" begin
    # Spacecraft context
    sat1 = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
    )

    # Maneuver models
    toi = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
    moi = ImpulsiveManeuver(axes = VNB(), element1 = 0.4, element2 = 0.5, element3 = 0.6)

    # Build calcs
    toi_calc = ManeuverCalc(toi, sat1, DeltaVVector())
    moi_calc = ManeuverCalc(moi, sat1, DeltaVVector())

    # Construct SolverVariables with new interface
    var_toi = SolverVariable(
        calc = toi_calc,
        name = "toi",
        lower_bound = [-10.0, 0.0, 0.0],
        upper_bound = [ 10.0, 0.0, 0.0],
    )

    var_moi = SolverVariable(
        calc = moi_calc,
        name = "moi",
        lower_bound = [-10.0, 0.0, 0.0],
        upper_bound = [ 10.0, 0.0, 0.0],
    )

    # Basic assertions (construction only)
    @test var_toi.numvars == 3
    @test var_toi.name == "toi"
    @test var_toi.calc isa ManeuverCalc
    @test var_toi.lower_bound == [-10.0, 0.0, 0.0]
    @test var_toi.upper_bound == [ 10.0, 0.0, 0.0]

    @test var_moi.numvars == 3
    @test var_moi.name == "moi"
    @test var_moi.calc isa ManeuverCalc
    @test var_moi.lower_bound == [-10.0, 0.0, 0.0]
    @test var_moi.upper_bound == [ 10.0, 0.0, 0.0]
end

@testset "SolverVariable: OrbitCalc IncomingAsymptoteFull" begin
    # Build a spacecraft in IncomingAsymptoteState
    ia0 = IncomingAsymptoteState(7000.0, 0.01, 1.0, 0.5, 0.2, 0.0)  # [rp, c3, rla, dla, bpa, ta]
    sc  = Spacecraft(state = ia0, time = Time("2020-09-21T12:23:12", TAI(), ISOT()))

    calc = OrbitCalc(sc, IncomingAsymptote())
    sv   = SolverVariable(calc = calc, name = "ia-full",
                          lower_bound = fill(-Inf, 6), upper_bound = fill(Inf, 6))

    @test sv.numvars == 6
    @test get_calc(sv.calc) == to_vector(ia0)

    # Set new asymptote values and verify spacecraft state updated
    ia_new = [7100.0, 0.02, 1.1, 0.6, 0.25, 0.05]
    set_calc!(sv.calc, ia_new)

    vals = get_calc(sv.calc)
    @test length(vals) == 6
    @test all(isapprox.(vals, ia_new; atol=1e-10))

    sc_ia = get_state(sc, IncomingAsymptote())
    @test all(isapprox.(to_vector(sc_ia), ia_new; atol=1e-10))
end

@testset "SolverVariable: ManeuverCalc DeltaVVector" begin
    # Spacecraft context (not mutated by set_calc!)
    sc = Spacecraft(state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
                    time  = Time("2020-09-21T12:23:12", TAI(), ISOT()))

    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
    calc = ManeuverCalc(man, sc, DeltaVVector())
    sv   = SolverVariable(calc = calc, name = "dv-full",
                          lower_bound = [-10.0, 0.0, 0.0], upper_bound = [10.0, 0.0, 0.0])

    @test sv.numvars == 3
    @test get_calc(sv.calc) ≈ [0.1, 0.2, 0.3]

    set_calc!(sv.calc, [0.2, 0.3, 0.4])
    @test get_calc(sv.calc) ≈ [0.2, 0.3, 0.4]
    @test man.element1 ≈ 0.2 && man.element2 ≈ 0.3 && man.element3 ≈ 0.4
end

@testset "SolverVariable: BodyCalc GravParam" begin
    # Local mutable body with ASCII mu field

    b = CelestialBody(name = "Moon", mu = 4000.0)

    calc = BodyCalc(b, GravParam())
    sv   = SolverVariable(calc = calc, name = "mu", lower_bound = -Inf, upper_bound = Inf)

    @test sv.numvars == 1
    @test get_calc(sv.calc) ≈ b.mu
    set_calc!(sv.calc, 4.0e14)
    @test b.mu == 4.0e14
    @test get_calc(sv.calc) == 4.0e14
end

nothing