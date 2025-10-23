using Test
using LinearAlgebra
using AstroEpochs, AstroStates, AstroCoords, AstroUniverse
using AstroModels: Spacecraft, to_posvel
using AstroFun: OrbitCalc, PosMag, TA, PosX, get_calc
using AstroProp: OrbitPropagator, StopAt, PointMassGravity, ForceModel, IntegratorConfig, propagate

# Helper: fresh spacecraft for each test
make_sat() = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
    name  = "SC-StopAt",
    coord_sys = CoordinateSystem(earth, ICRFAxes()),
)

# Forces + integrator (same for all tests)
gravity = PointMassGravity(earth, (moon, sun))
forces  = ForceModel(gravity)
integ   = IntegratorConfig(Tsit5(); dt=10.0, reltol=1e-9, abstol=1e-9)
prop    = OrbitPropagator(forces, integ)

# Tolerances for event satisfaction
const POS_TOL = 1e-3      # km
const ANG_TOL = 1e-6      # rad

#@testset "StopAt with OrbitCalc variables" begin
    @testset "Stop at radius magnitude target (any direction)" begin
        sat = make_sat()
        target_r = 7000.0
        sol = propagate(prop, sat, StopAt(sat, PosMag(), target_r; direction = 0))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)

        rmag = norm(to_posvel(sat)[1:3])
        @test isapprox(rmag, target_r; atol = POS_TOL)
    end

    #@testset "Stop at true anomaly = π (decreasing crossing)" begin
    #    sat = make_sat()
    #    sol = propagate(prop, sat, StopAt(sat, TA(), pi; direction = -1))
    #    @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
    #
    #
#
#        # Validate via OrbitCalc(TA)
#        ta_now = get_calc(OrbitCalc(sat, TA()))
#        # Normalize angle error to nearest branch around π
#        err = atan(sin(ta_now - pi), cos(ta_now - pi))
#        @test abs(err) ≤ ANG_TOL
#    end

    @testset "Stop at x-position crossing (increasing)" begin
        sat = make_sat()
        target_x = 7.5
        sol = propagate(prop, sat, StopAt(sat, PosX(), target_x; direction = +1))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)

        posvel = to_posvel(sat)
        x = posvel[1]
        vx = posvel[4]
        @test isapprox(x, target_x; atol = POS_TOL)
        @test vx > 0                      # increasing crossing
    end
#end

nothing