using Test
using LinearAlgebra
using AstroFun
using AstroStates
using AstroEpochs

# Helper to normalize get_calc return to scalar if it’s a length-1 vector
scalarize(x) = x
scalarize(x::AbstractVector) = (length(x) == 1 ? x[1] : error("expected length-1 vector"))

# Fresh spacecraft for isolation
function make_sc()
    Spacecraft(
        state = SphericalRADECState(CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])),
        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
    )
end

# Rationale: convert_orbitcalc_state should fast-path when already the requested concrete type.
@testset "convert_orbitcalc_state: fast-path no conversion" begin
    cs = CoordinateSystem(earth, Inertial())
    cart = CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0])
    os = OrbitState(to_vector(cart), Cartesian())
    st = AstroFun.convert_orbitcalc_state(os, cs, Cartesian())
    @test st isa CartesianState
    @test to_vector(st) == to_vector(cart)
end

@testset "Orbit Calc Tests - Per Calc" begin
    # RAAN
    @testset "RAAN" begin
        sc = make_sc()
        calc = OrbitCalc(sc, RAAN())
        val = scalarize(get_calc(calc))
        # Test calc is set correctly
        sc_kep = get_state(sc,Keplerian())
        @test isapprox(scalarize(get_calc(calc)), sc_kep.raan; atol=1e-12)
        set_calc!(calc, pi/4)
        sc_kep = get_state(sc,Keplerian())
        # Test set_calc! updates spacecraft state
        @test isapprox(scalarize(get_calc(calc)), pi/4; atol=1e-12)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
        @test calc_input_statetag(calc.var) == Keplerian()
    end

    # PosMag
    @testset "PosMag" begin
        sc = make_sc()
        calc = OrbitCalc(sc, PosMag())
        # Test calc is set correctly (matches Cartesian |r|)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(scalarize(get_calc(calc)), norm(sc_cart.posvel[1:3]); atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, 10000.0)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(norm(sc_cart.posvel[1:3]), 10000.0; atol=1e-8)
        @test isapprox(scalarize(get_calc(calc)), 10000.0; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 1   
        @test calc_is_settable(calc.var) == true
    end

    # VelMag
    @testset "VelMag" begin
        sc = make_sc()
        calc = OrbitCalc(sc, VelMag())
        # Test calc is set correctly (matches Cartesian |v|)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(scalarize(get_calc(calc)), norm(sc_cart.posvel[4:6]); atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, 8.0)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(norm(sc_cart.posvel[4:6]), 8.0; atol=1e-8)
        @test isapprox(scalarize(get_calc(calc)), 8.0; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
        @test Base.invokelatest(calc_is_settable, VelMag()) === true
        @test Base.invokelatest(calc_numvars, VelMag()) === 1
    end

    # SMA
    @testset "SMA" begin
        sc = make_sc()
        calc = OrbitCalc(sc, SMA())
        # Test calc is set correctly (matches Keplerian sma)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(scalarize(get_calc(calc)), sc_kep.sma; atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, 10000.0)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(sc_kep.sma, 10000.0; atol=1e-8)
        @test isapprox(scalarize(get_calc(calc)), 10000.0; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
    end

    @testset "Ecc" begin
        sc = make_sc()
        calc = OrbitCalc(sc, Ecc())
        # Test calc is set correctly (matches Keplerian ecc)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(scalarize(get_calc(calc)), sc_kep.ecc; atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, 0.02)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(sc_kep.ecc, 0.02; atol=1e-8)
        @test isapprox(scalarize(get_calc(calc)), 0.02; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
    end

    # OutGoingRLA
    @testset "OutGoingRLA" begin
        sc = make_sc()
        calc = OrbitCalc(sc, OutGoingRLA())
        # Test calc is set correctly (matches OutGoingAsymptote rla)
        sc_out = get_state(sc, OutGoingAsymptote())
        @test isapprox(scalarize(get_calc(calc)), sc_out.rla; atol=1e-12)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, pi/3)
        sc_out = get_state(sc, OutGoingAsymptote())
        @test isapprox(sc_out.rla, pi/3; atol=1e-12)
        @test isapprox(scalarize(get_calc(calc)), pi/3; atol=1e-12)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
    end

    # PositionVector
    @testset "PositionVector" begin
        sc = make_sc()
        calc = OrbitCalc(sc, PositionVector())
        # Test calc is set correctly (matches Cartesian position vector)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(get_calc(calc), sc_cart.posvel[1:3]; atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        target = [7000.0, 300.0, 0.0]
        set_calc!(calc, target)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(sc_cart.posvel[1:3], target; atol=1e-8)
        @test isapprox(get_calc(calc), target; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 3
        @test calc_is_settable(calc.var) == true
    end

    # TA
    @testset "TA" begin
        sc = make_sc()
        calc = OrbitCalc(sc, TA())
        # Test calc is set correctly (matches Keplerian ta)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(scalarize(get_calc(calc)), sc_kep.ta; atol=1e-12)
        # Update via set_calc! and verify spacecraft state updated
        set_calc!(calc, pi/6)
        sc_kep = get_state(sc, Keplerian())
        @test isapprox(sc_kep.ta, pi/6; atol=1e-12)
        @test isapprox(scalarize(get_calc(calc)), pi/6; atol=1e-12)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
    end

    # IncomingAsymptoteFull
    @testset "IncomingAsymptote" begin
        sc = make_sc()
        calc = OrbitCalc(sc, IncomingAsymptote())
        # Test calc is set correctly (matches IncomingAsymptote vector)
        sc_in = get_state(sc, IncomingAsymptote())
        expected = [sc_in.rp, sc_in.c3, sc_in.rla, sc_in.dla, sc_in.bpa, sc_in.ta]
        @test isapprox(get_calc(calc), expected; atol=1e-8)
        # Update via set_calc! and verify spacecraft state updated
        target = [ 30000.0, 5.0, 0.39269908169872414, 0.3141592653589793, 0.7853981633974483, 0.001]
        set_calc!(calc, target)
        sc_in = get_state(sc, IncomingAsymptote())
        expected2 = [sc_in.rp, sc_in.c3, sc_in.rla, sc_in.dla, sc_in.bpa, sc_in.ta]
        @test isapprox(expected2, target; atol=1e-8)
        @test isapprox(get_calc(calc), target; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 6
        @test calc_is_settable(calc.var) == true
    end

    @testset "Keplerian (full state vector) OrbitCalc" begin
        sc = make_sc()
        calc = OrbitCalc(sc, Keplerian())
        # Test: get_calc returns the full Keplerian state vector
        kvec = get_calc(calc)
        @test length(kvec) == 6
        sc_kep = get_state(sc, Keplerian())
        expected = [sc_kep.sma, sc_kep.ecc, sc_kep.inc, sc_kep.raan, sc_kep.aop, sc_kep.ta]
        @test isapprox(kvec, expected; rtol=1e-13)

        # Test: set_calc! updates the spacecraft Keplerian state
        target = [8000.0, 0.05, 0.2, 1.1, 0.7, 0.3]
        set_calc!(calc, target)
        sc_kep2 = get_state(sc, Keplerian())
        expected2 = [sc_kep2.sma, sc_kep2.ecc, sc_kep2.inc, sc_kep2.raan, sc_kep2.aop, sc_kep2.ta]
        @test isapprox(expected2, target; rtol=1e-13)
        @test isapprox(get_calc(calc), target; rtol=1e-13)

        # Trait checks
        @test calc_numvars(Keplerian()) == 6
        @test calc_is_settable(Keplerian()) == true
        @test calc_input_statetag(Keplerian()) == Keplerian()
    end

    # VelocityVector
    @testset "VelocityVector" begin
        sc = make_sc()
        calc = OrbitCalc(sc, VelocityVector())
        # Test calc is set correctly (matches Cartesian velocity vector)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(get_calc(calc), sc_cart.posvel[4:6]; atol=1e-12)
        # Update via set_calc! and verify spacecraft state updated
        target = [0.0, 7.5, 1.0]
        set_calc!(calc, target)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(sc_cart.posvel[4:6], target; atol=1e-12)
        @test isapprox(get_calc(calc), target; atol=1e-12)
        # Test trait methods
        @test calc_numvars(calc.var) == 3
        @test calc_is_settable(calc.var) == true
    end

    # PosX
    @testset "PosX" begin
        sc = make_sc()
        calc = OrbitCalc(sc, PosX())
        # Read: matches Cartesian x component
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(scalarize(get_calc(calc)), sc_cart.posvel[1]; atol=1e-8)
        # Write: update x and verify spacecraft state and calc reflect it
        target_x = 7100.0
        set_calc!(calc, target_x)
        sc_cart = get_state(sc, Cartesian())
        @test isapprox(sc_cart.posvel[1], target_x; atol=1e-8)
        @test isapprox(scalarize(get_calc(calc)), target_x; atol=1e-8)
        # Test trait methods
        @test calc_numvars(calc.var) == 1
        @test calc_is_settable(calc.var) == true
    end

    

end
nothing