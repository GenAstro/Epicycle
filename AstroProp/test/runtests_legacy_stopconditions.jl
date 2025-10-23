using Test
using LinearAlgebra
using DifferentialEquations

using AstroBase
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroCoords
using AstroModels
using AstroMan
using AstroFun
using AstroProp

using AstroEpochs, AstroStates
using AstroUniverse
using AstroModels: Spacecraft
using AstroProp
using AstroFun: OrbitCalc, PosMag, get_calc

# Helpers
make_sat_cart() = Spacecraft(
    state = CartesianState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03]),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
)

make_sat_orbit() = Spacecraft(
    state = OrbitState([7000.0, 300.0, 0.0, 0.0, 7.5, 0.03], Cartesian()),
    time  = Time("2015-09-21T12:23:12", TAI(), ISOT()),
)

# Earth-only point-mass gravity for deterministic tests
forces = ForceModel(PointMassGravity(earth, ()))
# Fast integrator/tolerances for tests
integ  = IntegratorConfig(Tsit5(); abstol = 1e-9, reltol = 1e-9, dt = 20.0)

# Utility to accept Success or Terminated from SciML
ok_ret(rc) = rc in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)

# Basic vector access
#to_posvel(sat::Spacecraft) = AstroStates.to_vector(sat.state)

const TIME_SEC_ATOL = 1e-3     # seconds
const TIME_DAY_ATOL = 1e-6     # days
const POS_ATOL      = 1e-3     # km
const RVEL_ATOL     = 1e-6     # km/s
const Z_ATOL        = 1e-6     # km

@testset "Legacy StopAt* interface works" begin
    @testset "StopAtSeconds(3600) advances ~1 hour" begin
        sat = make_sat_cart()
        t0 = sat.time.jd
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtSeconds(3600.0))
        @test ok_ret(sol.retcode)
        Δt_sec = (sat.time.jd - t0) * 86400.0
        @test isapprox(Δt_sec, 3600.0; atol = TIME_SEC_ATOL)
    end

    @testset "StopAtDays(10) advances ~10 days" begin
        sat = make_sat_orbit()
        t0 = sat.time.jd
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtDays(10.0))
        @test ok_ret(sol.retcode)
        Δt_days = (sat.time.jd - t0)
        @test isapprox(Δt_days, 10.0; atol = TIME_DAY_ATOL)
    end

    @testset "StopAtApoapsis hits radial-velocity zero (decreasing)" begin
        sat = make_sat_cart()
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtApoapsis(sat))
        @test ok_ret(sol.retcode)
        r, v = to_posvel(sat)[1:3], to_posvel(sat)[4:6]
        rv = dot(r, v) / norm(r)
        @test isapprox(rv, 0.0; atol = RVEL_ATOL)
        # Slightly advance to confirm decreasing crossing (rv < 0)
        sol2 = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtSeconds(1.0))
        r2, v2 = to_posvel(sat)[1:3], to_posvel(sat)[4:6]
        rv2 = dot(r2, v2) / norm(r2)
        @test rv2 ≤ 0
    end

    @testset "StopAtPeriapsis hits radial-velocity zero (increasing)" begin
        sat = make_sat_cart()
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtPeriapsis(sat))
        @test ok_ret(sol.retcode)
        r, v = to_posvel(sat)[1:3], to_posvel(sat)[4:6]
        rv = dot(r, v) / norm(r)
        @test isapprox(rv, 0.0; atol = RVEL_ATOL)
        # Slightly advance to confirm increasing crossing (rv > 0)
        sol2 = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtSeconds(1.0))
        r2, v2 = to_posvel(sat)[1:3], to_posvel(sat)[4:6]
        rv2 = dot(r2, v2) / norm(r2)
        @test rv2 ≥ 0
    end

    @testset "StopAtAscendingNode crosses z=0 with vz>0" begin
        sat = make_sat_cart()
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtAscendingNode(sat))
        @test ok_ret(sol.retcode)
        x, y, z, vx, vy, vz = to_posvel(sat)
        @test isapprox(z, 0.0; atol = Z_ATOL)
        @test vz ≥ 0
    end

    @testset "StopAtRadius hits requested radius" begin
        sat = make_sat_cart()
        target_r = 7100.0
        sol = propagate(DynSys(spacecraft=[sat], forces=forces), integ, StopAtRadius(sat, target_r))
        @test ok_ret(sol.retcode)
        r = norm(to_posvel(sat)[1:3])
        @test isapprox(r, target_r; atol = POS_ATOL)
    end
end

nothing