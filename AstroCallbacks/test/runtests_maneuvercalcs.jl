using Test
using LinearAlgebra
using AstroFun
using AstroStates
using AstroEpochs

# Helper to normalize get_calc return to scalar if it’s a length-1 vector
#scalarize(x) = x
#scalarize(x::AbstractVector) = (length(x) == 1 ? x[1] : error("expected length-1 vector"))

# Fresh spacecraft for isolation
#function make_sc()
#    Spacecraft(
#        state = SphericalRADECState(CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])),
#        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
#    )
#end

@testset "DeltaVVector (ManeuverCalc)" begin
    sc = make_sc()
    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
    calc = ManeuverCalc(man, sc, DeltaVVector())

    # Initial get matches maneuver fields
    @test all(isapprox.(get_calc(calc), [0.1, 0.2, 0.3]; atol=1e-8))

    # Update via set_calc! and verify maneuver updated (spacecraft not asserted here)
    set_calc!(calc, [0.2, 0.3, 0.4])
    @test all(isapprox.(get_calc(calc), [0.2, 0.3, 0.4]; atol=1e-8))
    @test man.element1 ≈ 0.2 && man.element2 ≈ 0.3 && man.element3 ≈ 0.4
    # Test trait methods
    @test calc_numvars(calc.var) == 3
    @test calc_is_settable(calc.var) == true
end

@testset "DeltaVMag (ManeuverCalc)" begin
    # Use existing helper if available
     sc = make_sc()

    # Maneuver with known components
    man = ImpulsiveManeuver(axes = VNB(), element1 = 0.1, element2 = 0.2, element3 = 0.3)
    calc = ManeuverCalc(man, sc, DeltaVMag())

    # get_calc returns a length-1 vector with ||Δv||
    vals = get_calc(calc)
    @test length(vals) == 1
    @test isapprox(vals[1], sqrt(0.1^2 + 0.2^2 + 0.3^2); atol=1e-12)

    # Not settable: set_calc! should throw with clear error message
    @test_throws "Variable DeltaVMag is not settable" set_calc!(calc, [0.0])

    # Test trait methods
    @test calc_numvars(calc.var) == 1
    @test calc_is_settable(calc.var) == false
end

# Rationale: ManeuverCalc get_calc for DeltaV variables reads from the ImpulsiveManeuver.
@testset "ManeuverCalc variable evaluation" begin
    m = ImpulsiveManeuver(axes=Inertial(), Isp=300.0, element1=0.01, element2=0.02, element3=-0.03)
    sc = Spacecraft(
        state=CartesianState([7000.0,0,0, 0,7.5,0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    mc_vec = ManeuverCalc(m, sc, DeltaVVector())
    mc_mag = ManeuverCalc(m, sc, DeltaVMag())

    @test get_calc(mc_vec) ≈ [0.01, 0.02, -0.03]
    @test scalarize(get_calc(mc_mag)) ≈ sqrt(0.01^2 + 0.02^2 + 0.03^2)
end

nothing