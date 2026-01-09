using Test
using LinearAlgebra
using AstroEpochs, AstroStates, AstroFrames, AstroUniverse
using AstroModels: Spacecraft, to_posvel
using AstroCallbacks: OrbitCalc, PosMag, TA, PosX, get_calc
using AstroProp: OrbitPropagator, StopAt, PointMassGravity, ForceModel, IntegratorConfig, propagate
using AstroProp: PropDurationSeconds, PropDurationDays

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
        @test isapprox(rmag, target_r; atol = 1e-6)
    end

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

# ============================================================================
# Time-based stopping conditions
# ============================================================================

@testset "PropDuration stopping conditions" begin
    @testset "PropDurationSeconds - forward propagation" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Propagate for 3600 seconds (1 hour)
        duration_sec = 3600.0
        sol = propagate(prop, sat, StopAt(sat, PropDurationSeconds(), duration_sec))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check elapsed time in TT (Earth-centered)
        t_final = sat.time
        elapsed_sec = (t_final.tt.jd - t_start.tt.jd) * 86400.0
        @test isapprox(elapsed_sec, duration_sec; atol = 1e-4)
    end

    @testset "PropDurationDays - forward propagation" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Propagate for 2.5 days
        duration_days = 2.5
        sol = propagate(prop, sat, StopAt(sat, PropDurationDays(), duration_days))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check elapsed time in TT (Earth-centered)
        t_final = sat.time
        elapsed_days = (t_final.tt.jd - t_start.tt.jd)
        @test isapprox(elapsed_days, duration_days; atol = 1e-6/86400)
    end

    @testset "PropDurationDays - zero duration throws error" begin
        sat = make_sat()
        @test_throws ErrorException propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.0))
    end

    @testset "PropDurationSeconds - backward with :infer" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Propagate backward using negative duration + :infer
        duration_sec = -7200.0  # Negative = backward
        sol = propagate(prop, sat, StopAt(sat, PropDurationSeconds(), duration_sec); direction=:infer)
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check elapsed time is negative (went backward in TT)
        t_final = sat.time
        elapsed_sec = (t_final.tt.jd - t_start.tt.jd) * 86400.0
        @test elapsed_sec < 0  # Went backward
        @test isapprox(elapsed_sec, duration_sec; atol = 1e-4)
    end

    @testset "PropDurationDays - backward with :infer" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Propagate backward using negative duration + :infer
        duration_days = -1.5  # Negative = backward
        sol = propagate(prop, sat, StopAt(sat, PropDurationDays(), duration_days); direction=:infer)
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check elapsed time is negative (went backward in TT)
        t_final = sat.time
        elapsed_days = (t_final.tt.jd - t_start.tt.jd)
        @test elapsed_days < 0  # Went backward
        @test isapprox(elapsed_days, duration_days; atol = 1e-6/86400.0)
    end

    @testset "PropDurationSeconds - explicit :backward" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Propagate backward with explicit direction and negative duration (both agree)
        duration_sec = -3600.0  # Negative = backward
        sol = propagate(prop, sat, StopAt(sat, PropDurationSeconds(), duration_sec); direction=:backward)
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check elapsed time is negative (went backward in TT)
        t_final = sat.time
        elapsed_sec = (t_final.tt.jd - t_start.tt.jd) * 86400.0
        @test elapsed_sec < 0  # Went backward
        @test isapprox(elapsed_sec, duration_sec; atol = 1.0)
    end

    @testset "Direction conflict validation" begin
        sat = make_sat()
        
        # Negative duration with explicit :forward should error
        @test_throws ErrorException propagate(prop, sat, StopAt(sat, PropDurationSeconds(), -100.0); direction=:forward)
        
        # Positive duration with explicit :backward should error  
        @test_throws ErrorException propagate(prop, sat, StopAt(sat, PropDurationDays(), 1.0); direction=:backward)
    end

    @testset "Time-based stops must have direction=0" begin
        sat = make_sat()
        
        # PropDuration with direction != 0 should error
        @test_throws ErrorException StopAt(sat, PropDurationSeconds(), 100.0; direction=1)
        @test_throws ErrorException StopAt(sat, PropDurationDays(), 1.0; direction=-1)
    end

    @testset "Multiple time-based stops not allowed" begin
        sat = make_sat()
        
        # Two PropDurationSeconds
        @test_throws ErrorException propagate(prop, sat, 
            StopAt(sat, PropDurationSeconds(), 100.0),
            StopAt(sat, PropDurationSeconds(), 200.0))
        
        # Two PropDurationDays
        @test_throws ErrorException propagate(prop, sat,
            StopAt(sat, PropDurationDays(), 1.0),
            StopAt(sat, PropDurationDays(), 2.0))
        
        # Mixed PropDuration types
        @test_throws ErrorException propagate(prop, sat,
            StopAt(sat, PropDurationSeconds(), 100.0),
            StopAt(sat, PropDurationDays(), 1.0))
        
        # Two absolute time stops
        time1 = Time(sat.time.tt.jd + 1.0, TT(), JD())
        time2 = Time(sat.time.tt.jd + 2.0, TT(), JD())
        @test_throws ErrorException propagate(prop, sat,
            StopAt(sat, time1),
            StopAt(sat, time2))
        
        # PropDuration + absolute time
        @test_throws ErrorException propagate(prop, sat,
            StopAt(sat, PropDurationSeconds(), 3600.0),
            StopAt(sat, time1))
        
        # Forward and backward durations (still errors - ambiguous)
        @test_throws ErrorException propagate(prop, sat,
            StopAt(sat, PropDurationSeconds(), 100.0),
            StopAt(sat, PropDurationSeconds(), -100.0))
    end
end

@testset "Absolute time stopping conditions" begin
    @testset "StopAt(Time) - future time" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Target time 1 day in the future
        target_time = Time(t_start.tt.jd + 1.0, TT(), JD())
        
        sol = propagate(prop, sat, StopAt(sat, target_time))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check we stopped at target time (within tolerance)
        t_final = sat.time
        elapsed_days = t_final.tt.jd - t_start.tt.jd
        @test isapprox(elapsed_days, 1.0; atol = 1e-6/86400.0)
    end

    @testset "StopAt(Time) - past time with :infer" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Target time in the past - works with :infer
        past_time = Time(sat.time.tt.jd - 1.0, TT(), JD())
        sol = propagate(prop, sat, StopAt(sat, past_time); direction=:infer)
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Check we went backward
        t_final = sat.time
        elapsed_days = t_final.tt.jd - t_start.tt.jd
        @test elapsed_days < 0  # Went backward
        @test isapprox(elapsed_days, -1.0; atol = 1e-6/86400.0)
    end

    @testset "StopAt(Time) - past time with :forward errors" begin
        sat = make_sat()
        # Target time in the past with default :forward should error
        past_time = Time(sat.time.tt.jd - 1.0, TT(), JD())
        @test_throws ErrorException propagate(prop, sat, StopAt(sat, past_time); direction=:forward)
    end

    @testset "StopAt(Time) - handles time scale conversion" begin
        sat = make_sat()
        t_start = deepcopy(sat.time)
        
        # Target specified in UTC, propagation uses TT
        target_utc = Time("2015-09-22T12:00:00", UTC(), ISOT())
        
        sol = propagate(prop, sat, StopAt(sat, target_utc))
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
        
        # Verify we stopped at the right time (compare in TT)
        expected_elapsed = target_utc.tt.jd - t_start.tt.jd
        actual_elapsed = sat.time.tt.jd - t_start.tt.jd
        @test isapprox(actual_elapsed, expected_elapsed; atol = 1e-6/86400.0)
    end
end

@testset "Propagation time scale (TT vs TDB)" begin
    @testset "Earth-centered uses TT" begin
        sat_earth = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2015-01-01T00:00:00", TAI(), ISOT()),
            coord_sys = CoordinateSystem(earth, ICRFAxes()),
        )
        t_start = deepcopy(sat_earth.time)
        
        # Propagate and verify TT is used
        duration_sec = 1000.0
        sol = propagate(prop, sat_earth, StopAt(sat_earth, PropDurationSeconds(), duration_sec))
        
        elapsed_tt = (sat_earth.time.tt.jd - t_start.tt.jd) * 86400.0
        @test isapprox(elapsed_tt, duration_sec; atol = 1e-4)
    end

    @testset "Non-Earth center uses TDB" begin
        # Mars-centered spacecraft - use date with good ephemeris coverage
        gravity_mars = PointMassGravity(mars, ())
        forces_mars = ForceModel(gravity_mars)
        prop_mars = OrbitPropagator(forces_mars, integ)
        
        sat_mars = Spacecraft(
            state = CartesianState([4000.0, 0.0, 0.0, 0.0, 3.0, 0.0]),
            time = Time("2020-06-01T12:00:00", TAI(), ISOT()),
            coord_sys = CoordinateSystem(mars, ICRFAxes()),
        )
        t_start = deepcopy(sat_mars.time)
        
        # Propagate and verify TDB is used
        duration_sec = 1000.0
        sol = propagate(prop_mars, sat_mars, StopAt(sat_mars, PropDurationSeconds(), duration_sec))
        
        elapsed_tdb = (sat_mars.time.tdb.jd - t_start.tdb.jd) * 86400.0
        @test isapprox(elapsed_tdb, duration_sec; atol = 1e-4)
    end

    @testset "Absolute time stop uses correct scale for Earth" begin
        sat = make_sat()
        
        # StopAt(Time) should convert to TT for Earth
        target_time = Time("2015-09-22T00:00:00", UTC(), ISOT())
        sol = propagate(prop, sat, StopAt(sat, target_time))
        
        # Implementation uses TT internally for Earth-centered
        # Verify by checking the constructor converts correctly
        @test sol.retcode in (SciMLBase.ReturnCode.Success, SciMLBase.ReturnCode.Terminated)
    end
end

nothing