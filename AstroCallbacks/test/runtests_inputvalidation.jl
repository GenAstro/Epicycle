using Test
using LinearAlgebra

using AstroCallbacks
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroFrames
using AstroModels

# Rationale: Unknown OrbitVar must throw via calc_input_statetag fallback.
struct UnknownVar <: AbstractOrbitVar end
@testset "OrbitCalc unknown variable dispatch" begin
    sc = Spacecraft(
        state = CartesianState([1,2,3, 4,5,6]),
        time  = Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
        oc = OrbitCalc(sc, UnknownVar())
    try
        get_calc(oc)
        @test false  # should not reach
    catch e
        msg = sprint(showerror, e)
        @test occursin("There is no AbstractOrbitVar defined named UnknownVar()", msg)
        @test occursin("subtypes(AbstractOrbitVar)", msg)
    end
end

struct NoMuOrigin <: AbstractPoint end
@testset "convert_orbitcalc_state: missing μ triggers error" begin
    cs_nomu = CoordinateSystem(NoMuOrigin(), Inertial())
    cart = CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0])
    os = OrbitState(to_vector(cart), Cartesian())
    @test_throws ErrorException AstroCallbacks.convert_orbitcalc_state(os, cs_nomu, Keplerian())
end

# Rationale: convert_orbitcalc_state: Cartesian→Keplerian succeeds 
# when μ available from earth.
@testset "convert_orbitcalc_state: Cartesian → Keplerian with μ" begin
    cs = CoordinateSystem(earth, Inertial())
    cart = CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0])
    os = OrbitState(to_vector(cart), Cartesian())
    kep = AstroCallbacks.convert_orbitcalc_state(os, cs, Keplerian())
    @test kep isa KeplerianState
end

# Rationale: Constraint positional/keyword constructors validate lengths and promote eltype.
@testset "Constraint constructors and type promotion" begin
    sc = Spacecraft(
        state=CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    oc = OrbitCalc(sc, VelocityVector())

    # Positional constructor: 3-var vector constraint, BigFloat promotion
    lb = big.([-1.0, -1.0, -1.0]); ub = big.([1.0, 1.0, 1.0]); sca = big.([1.0, 1.0, 1.0])
    c = Constraint(oc, lb, ub, sca)
    @test eltype(c.lower_bounds) === BigFloat
    @test c.numvars == 3

    # Keyword constructor infers numvars from variable
    ck = Constraint(calc=oc, lower_bounds=[-1.0,-1.0,-1.0], upper_bounds=[1.0,1.0,1.0], scale=[1.0,1.0,1.0])
    @test ck.numvars == 3

    # Length mismatches throw
    @test_throws ArgumentError Constraint(oc, [-1.0], [1.0,1.0], [1.0,1.0], 2)
    #@test_throws ArgumentError Constraint(calc=oc, lower_bounds=[-1.0, -1.0], upper_bounds=[1.0,1.0], scale=[1.0], numvars=2)
end

nothing