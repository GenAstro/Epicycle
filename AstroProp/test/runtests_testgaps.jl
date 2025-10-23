using Test
using LinearAlgebra
using SciMLBase

using AstroEpochs
using AstroStates
using AstroUniverse
using AstroModels: Spacecraft, to_posvel, set_posvel!
using AstroFun: OrbitCalc, PosMag
using AstroProp

# Helpers
make_sat() = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
)

forces_earth_only() = ForceModel(PointMassGravity(earth, ()))
integ_fast() = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)

# 1) find_center branches: return nothing and conflict error
@testset "find_center branches" begin
    # A tiny dummy OrbitODE for exercising the 'no centers' branch
    struct _DummyForce <: AstroProp.OrbitODE end
    fm_none = ForceModel((_DummyForce(),))
    @test fm_none.center === nothing  # covers: return nothing

    # Conflicting centers error
    pm1 = PointMassGravity(earth, ())
    pm2 = PointMassGravity(moon,  ())
    @test_throws ErrorException ForceModel((pm1, pm2))  # covers: error("Multiple conflicting...")
end

# 2) Internal fallbacks/errors in orbit_propagator.jl


# 4) Legacy StopAt* backward-prop terminate branches (integrator.dt < 0)
@testset "Legacy StopAt* terminate on backward propagation" begin
    sat1 = make_sat()
    sat2 = make_sat()
    sat3 = make_sat()
    prop_forces = forces_earth_only()
    integ = integ_fast()

    sol1 = propagate(DynSys(spacecraft=[sat1], forces=prop_forces), integ, StopAtSeconds(-10.0); direction = :backward)
    @test sol1.retcode in (ReturnCode.Success, ReturnCode.Terminated)

    sol2 = propagate(DynSys(spacecraft=[sat2], forces=prop_forces), integ, StopAtDays(-1/86400); direction = :backward)
    @test sol2.retcode in (ReturnCode.Success, ReturnCode.Terminated)

    sol3 = propagate(DynSys(spacecraft=[sat3], forces=prop_forces), integ, StopAtRadius(sat3, 7100.0); direction = :backward)
    @test sol3.retcode in (ReturnCode.Success, ReturnCode.Terminated)
end

nothing