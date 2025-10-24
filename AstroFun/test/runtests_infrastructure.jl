using Test

using AstroFun
using AstroMan
using AstroCoords
using AstroUniverse
using AstroModels
using AstroStates

@testset "OrbitCalc PositionVector: get/set round-trip" begin
    sc = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time  = Time("2020-09-21T12:23:12", TAI(), ISOT()),
    )
    oc = OrbitCalc(sc, PositionVector())

    # get_calc returns current position
    p = get_calc(oc)
    @test p == [7000.0, 300.0, 0.0]

    # set_calc! updates spacecraft state position only
    set_calc!(oc, [7001.0, 299.0, -1.0])
    @test get_calc(oc) == [7001.0, 299.0, -1.0]
    @test to_vector(sc.state)[4:6] == [0.0, 7.5, 1.0]
end

# Rationale: _subjects_from_calc returns the right tuple of subjects per Calc type.
@testset "subjects_from_calc" begin
    sc = Spacecraft(
        state=CartesianState([7000.0,0,0, 0,7.5,0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()))
  
    oc = OrbitCalc(sc, PositionVector())
    bc = BodyCalc(earth, GravParam())
    mc = ManeuverCalc(ImpulsiveManeuver(), sc, DeltaVMag())

    @test AstroFun._subjects_from_calc(oc) == (sc,)
    @test AstroFun._subjects_from_calc(bc) == (earth,)
    @test AstroFun._subjects_from_calc(mc) == (mc.man, sc)
end

# Rationale: ManeuverCalc get_calc for DeltaV variables reads from the ImpulsiveManeuver.
@testset "ManeuverCalc variable evaluation" begin
    m = ImpulsiveManeuver(axes=Inertial(), Isp=300.0, 
          element1=0.01, element2=0.02, element3=-0.03)
    sc = Spacecraft(
        state=CartesianState([7000.0,0,0, 0,7.5,0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    mc_vec = ManeuverCalc(m, sc, DeltaVVector())
    mc_mag = ManeuverCalc(m, sc, DeltaVMag())

    @test get_calc(mc_vec) ≈ [0.01, 0.02, -0.03]
    @test first(get_calc(mc_mag)) ≈ sqrt(0.01^2 + 0.02^2 + 0.03^2)
end

# Rationale: _subjects_from_calc returns the right tuple of subjects per Calc type.
@testset "subjects_from_calc" begin
    sc = Spacecraft(
        state=CartesianState([7000.0,0,0, 0,7.5,0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    oc = OrbitCalc(sc, PositionVector())
    bc = BodyCalc(earth, GravParam())
    mc = ManeuverCalc(ImpulsiveManeuver(), sc, DeltaVMag())

    @test AstroFun._subjects_from_calc(oc) == (sc,)
    @test AstroFun._subjects_from_calc(bc) == (earth,)
    @test AstroFun._subjects_from_calc(mc) == (mc.man, sc)
end

# Rationale: func_eval always returns a Vector, for scalar and vector-valued calcs.
@testset "func_eval output normalization" begin
    sc = Spacecraft(
        state=CartesianState([7000.0,300.0,0.0, 0.0,7.5,1.0]),
        time=Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    c_vec = Constraint(calc=OrbitCalc(sc, PositionVector()),
                       lower_bounds=[-1.0,-1.0,-1.0], upper_bounds=[1.0,1.0,1.0], scale=[1.0,1.0,1.0])
    v = func_eval(c_vec)
    @test v == [7000.0, 300.0, 0.0]

    c_sca = Constraint(calc=BodyCalc(earth, GravParam()),
                       lower_bounds=[0.0], upper_bounds=[1e7], scale=[1.0])
    s = func_eval(c_sca)
    @test s isa Vector
    @test length(s) == 1
end

# Rationale: Default trait calc_numvars(::AbstractCalcVariable) returns 1 for variables without overrides.
struct DummyVarNum <: AstroFun.AbstractCalcVariable end
@testset "calc_numvars default trait" begin
    @test AstroFun.calc_numvars(DummyVarNum()) == 1
end

# Rationale: convert_orbitcalc_state fast-path returns the identical object when target type matches.
@testset "convert_orbitcalc_state fast-path" begin
    cs = CoordinateSystem(earth, Inertial())
    st = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])
    out = AstroFun.convert_orbitcalc_state(st, cs, Cartesian())
    @test out === st
end

# Rationale: Default trait calc_is_settable(::AbstractCalcVariable) is false without an override.
struct DummyVarSettable <: AstroFun.AbstractCalcVariable end
@testset "calc_is_settable default trait" begin
    @test AstroFun.calc_is_settable(DummyVarSettable()) == false
end

# Rationale: Fallback set_calc! on non-settable variable throws a clear error mentioning the variable type.
@testset "set_calc! fallback throws for non-settable variable" begin
    m = ImpulsiveManeuver(axes=Inertial(), Isp=300.0, element1=0.01, element2=0.02, element3=-0.03)
    sc = Spacecraft(
        state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0]),
        time  = Time("2020-01-01T00:00:00", TAI(), ISOT()),
    )
    mc = ManeuverCalc(m, sc, DeltaVMag())
    try
        AstroFun.set_calc!(mc, 0.1)
        @test false
    catch e
        msg = sprint(showerror, e)
        @test msg == "Variable DeltaVMag is not settable."
    end
end

# Rationale: func_eval errors on unsupported calc return type (neither Number nor AbstractVector).
struct FakeCalc <: AstroFun.AbstractCalc end
AstroFun.get_calc(::FakeCalc) = [1.0 2.0; 3.0 4.0]  # Matrix triggers error path
@testset "func_eval unsupported calc return type" begin
    c = Constraint(calc=FakeCalc(), lower_bounds=[0.0], upper_bounds=[1.0], scale=[1.0])
    try
        AstroFun.func_eval(c)
        @test false
    catch e
        msg = sprint(showerror, e)
        @test occursin("unsupported calc return type", msg)
        @test occursin("Matrix", msg)
    end
end



# Rationale: calc_numvars default trait returns 1 for any variable without an override.
struct __DummyVarNum__ <: AstroFun.AbstractCalcVariable end
@testset "calc_numvars default trait = 1" begin
    @test AstroFun.calc_numvars(__DummyVarNum__()) == 1
end

# Rationale: to_concrete_state on an already-concrete AbstractState must be a pass-through.
@testset "to_concrete_state passthrough for concrete state" begin
    st = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])
    @test AstroFun.to_concrete_state(st) === st
end

# Rationale: convert_orbitcalc_state fast-path returns the identical instance when target type matches.
@testset "convert_orbitcalc_state fast-path (no-op) returns same instance" begin
    cs = CoordinateSystem(earth, Inertial())
    st = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])
    out = AstroFun.convert_orbitcalc_state(st, cs, Cartesian())
    @test out === st
end

# Rationale: calc_is_settable default trait is false unless a variable explicitly overrides it.
struct __DummyVarSettable__ <: AstroFun.AbstractCalcVariable end
@testset "calc_is_settable default trait = false" begin
    @test AstroFun.calc_is_settable(__DummyVarSettable__()) == false
end

@eval AstroFun begin
    __cov_fastpath__(st, cs) = convert_orbitcalc_state(st, cs, Cartesian())
end
@testset "convert_orbitcalc_state fast-path (no-op) returns same instance" begin
    cs = CoordinateSystem(earth, Inertial())
    st = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 1.0])
    out = AstroFun.__cov_fastpath__(st, cs)
    @test out === st
end

nothing