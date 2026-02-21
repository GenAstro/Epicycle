# Visual Regression Tests for Graphics Module
# 
# These tests require GLMakie with GPU/display access and should be run locally
# before releases, not in CI.
#
# Usage:
#   julia --project=. test/graphics/visual_regression.jl
#
# First run: Creates reference images in test/reference/
# Future runs: Compares against references, fails if different

using Test
using Epicycle
using ReferenceTests
using FileIO

@testset "Visual Regression Tests" begin
    
    @testset "Basic LEO orbit" begin
        # Simple circular LEO orbit
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "LEO Test"
        )
        
        # Propagate for one orbit
        gravity = PointMassGravity(earth, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 6000.0))
        
        # Create view and render
        view = View3D()
        add_spacecraft!(view, sat)
        display_view(view; size=(800, 600))
        
        # Save/compare reference
        @test_reference "../reference/leo_orbit.png" view._scene
    end
    
    @testset "Hohmann transfer with iterations" begin
        # Create spacecraft
        sat = Spacecraft(
            state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 7.5, 1.0),
            time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            name = "Hohmann"
        )
        
        # Force models and dynamics system
        pm_grav = PointMassGravity(earth, (moon, sun))
        forces = ForceModel(pm_grav)
        integ = IntegratorConfig(DP8(); abstol=1e-11, reltol=1e-11, dt=4000)
        dynsys = DynSys(forces=forces, spacecraft=[sat])
        
        # Create maneuvers
        toi = ImpulsiveManeuver(axes=VNB(), element1=0.1, element2=0.2, element3=0.3)
        moi = ImpulsiveManeuver(axes=VNB(), element1=0.4, element2=0.5, element3=0.6)
        
        # Define solver variables
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "toi",
            lower_bound = [0.0, 0.0, 0.0],
            upper_bound = [2.5, 0.0, 0.0]
        )
        
        var_moi = SolverVariable(
            calc = ManeuverCalc(moi, sat, DeltaVVector()),
            name = "moi",
            lower_bound = [0.0, 0.0, 0.0],
            upper_bound = [3.0, 0.0, 0.0]
        )
        
        # Define constraints
        pos_target = 45000.0
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0]
        )
        
        ecc_con = Constraint(
            calc = OrbitCalc(sat, Ecc()),
            lower_bounds = [0.0],
            upper_bounds = [0.0],
            scale = [1.0]
        )
        
        # Create events
        toi_event = Event(
            name = "TOI",
            event = () -> maneuver(sat, toi),
            vars = [var_toi],
            funcs = []
        )
        
        prop_event = Event(
            name = "Prop to Apoapsis",
            event = () -> propagate(dynsys, integ, StopAtApoapsis(sat))
        )
        
        moi_event = Event(
            name = "MOI",
            event = () -> maneuver(sat, moi),
            vars = [var_moi],
            funcs = [pos_con, ecc_con]
        )
        
        # Build sequence and solve
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event])
        add_events!(seq, moi_event, [prop_event])
        
        result = trajectory_solve(seq; record_iterations=true)
        
        # Create view with iterations shown
        view = View3D()
        add_spacecraft!(view, sat; show_iterations=true)
        display_view(view; size=(800, 600))
        
        # Save/compare reference
        @test_reference "../reference/hohmann_with_iterations.png" view._scene
    end
    
    @testset "Hohmann transfer without iterations" begin
        # Create spacecraft
        sat = Spacecraft(
            state = KeplerianState(7000.0, 0.001, 0.0, 0.0, 7.5, 1.0),
            time = Time("2020-09-21T12:23:12", TAI(), ISOT()),
            name = "Hohmann"
        )
        
        # Force models and dynamics system
        pm_grav = PointMassGravity(earth, (moon, sun))
        forces = ForceModel(pm_grav)
        integ = IntegratorConfig(DP8(); abstol=1e-11, reltol=1e-11, dt=4000)
        dynsys = DynSys(forces=forces, spacecraft=[sat])
        
        # Create maneuvers
        toi = ImpulsiveManeuver(axes=VNB(), element1=0.1, element2=0.2, element3=0.3)
        moi = ImpulsiveManeuver(axes=VNB(), element1=0.4, element2=0.5, element3=0.6)
        
        # Define solver variables
        var_toi = SolverVariable(
            calc = ManeuverCalc(toi, sat, DeltaVVector()),
            name = "toi",
            lower_bound = [0.0, 0.0, 0.0],
            upper_bound = [2.5, 0.0, 0.0]
        )
        
        var_moi = SolverVariable(
            calc = ManeuverCalc(moi, sat, DeltaVVector()),
            name = "moi",
            lower_bound = [0.0, 0.0, 0.0],
            upper_bound = [3.0, 0.0, 0.0]
        )
        
        # Define constraints
        pos_target = 45000.0
        pos_con = Constraint(
            calc = OrbitCalc(sat, PosMag()),
            lower_bounds = [pos_target],
            upper_bounds = [pos_target],
            scale = [1.0]
        )
        
        ecc_con = Constraint(
            calc = OrbitCalc(sat, Ecc()),
            lower_bounds = [0.0],
            upper_bounds = [0.0],
            scale = [1.0]
        )
        
        # Create events
        toi_event = Event(
            name = "TOI",
            event = () -> maneuver(sat, toi),
            vars = [var_toi],
            funcs = []
        )
        
        prop_event = Event(
            name = "Prop to Apoapsis",
            event = () -> propagate(dynsys, integ, StopAtApoapsis(sat))
        )
        
        moi_event = Event(
            name = "MOI",
            event = () -> maneuver(sat, moi),
            vars = [var_moi],
            funcs = [pos_con, ecc_con]
        )
        
        # Build sequence and solve
        seq = Sequence()
        add_events!(seq, prop_event, [toi_event])
        add_events!(seq, moi_event, [prop_event])
        
        result = trajectory_solve(seq; record_iterations=true)
        
        # Create view WITHOUT iterations shown
        view = View3D()
        add_spacecraft!(view, sat; show_iterations=false)
        display_view(view; size=(800, 600))
        
        # Save/compare reference
        @test_reference "../reference/hohmann_no_iterations.png" view._scene
    end
    
    @testset "DeepSpace spacecraft model" begin
        # Simple LEO orbit for DeepSpace model
        sat = Spacecraft(
            state = KeplerianState(30164.0, 0.0001, deg2rad(35.0), deg2rad(5.0), 3.07, 0.0),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "DeepSpace",
            cad_model = CADModel(
                file_path = joinpath(pkgdir(Epicycle), "assets", "DeepSpace1.obj"),
                scale = 1000.0,
                visible = true
            )
        )
        
        # Propagate for one orbit
        gravity = PointMassGravity(earth, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        propagate(prop, sat, StopAt(sat, PropDurationDays(), 0.39))
        
        # Create view and render
        view = View3D()
        add_spacecraft!(view, sat)
        display_view(view; size=(800, 600))
        
        # Save/compare reference
        @test_reference "../reference/deepspace_model.png" view._scene
    end
    
    @testset "Multiple spacecraft rendering" begin
        # Test visual rendering of multiple spacecraft in same view
        # Use different orbital planes so they're clearly visible
        sat1 = Spacecraft(
            state = KeplerianState(7000.0, 0.01, deg2rad(0.0), deg2rad(0.0), 0.0, 0.0),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "Equatorial"
        )
        sat2 = Spacecraft(
            state = KeplerianState(8500.0, 0.01, deg2rad(45.0), deg2rad(90.0), 0.0, 0.0),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "Inclined"
        )
        sat3 = Spacecraft(
            state = KeplerianState(10000.0, 0.01, deg2rad(90.0), deg2rad(180.0), 0.0, 0.0),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            name = "Polar"
        )
        
        # Propagate all three
        gravity = PointMassGravity(earth, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        
        propagate(prop, sat1, StopAt(sat1, PropDurationSeconds(), 5400.0))
        propagate(prop, sat2, StopAt(sat2, PropDurationSeconds(), 6000.0))
        propagate(prop, sat3, StopAt(sat3, PropDurationSeconds(), 6600.0))
        
        # Render all three together
        view = View3D()
        add_spacecraft!(view, sat1)
        add_spacecraft!(view, sat2)
        add_spacecraft!(view, sat3)
        display_view(view; size=(800, 600))
        
        @test_reference "../reference/multiple_spacecraft.png" view._scene
    end
    
    @testset "Body with texture loading failure" begin
        # Test that body renders with gray fallback when texture fails
        fake_body = CelestialBody(
            name = "BadTexture",
            mu = 398600.4418,
            equatorial_radius = 6378.137,
            texture_file = "nonexistent_texture.png"
        )
        
        coord_sys = CoordinateSystem(fake_body, ICRFAxes())
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            coord_sys = coord_sys
        )
        
        gravity = PointMassGravity(fake_body, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 3000.0))
        
        view = View3D(coord_sys=coord_sys)
        add_spacecraft!(view, sat)
        
        # Should warn but still render
        @test_logs (:warn, r"Failed to load texture.*Using solid gray") begin
            display_view(view; size=(800, 600))
        end
        
        @test_reference "../reference/body_no_texture.png" view._scene
    end
    
    @testset "Body without texture specified" begin
        # Test body with empty texture_file renders gray
        no_tex_body = CelestialBody(
            name = "GrayBody",
            mu = 398600.4418,
            equatorial_radius = 6378.137,
            texture_file = ""
        )
        
        coord_sys = CoordinateSystem(no_tex_body, ICRFAxes())
        sat = Spacecraft(
            state = CartesianState([7000.0, 0.0, 0.0, 0.0, 7.5, 0.0]),
            time = Time("2024-01-01T00:00:00", TAI(), ISOT()),
            coord_sys = coord_sys
        )
        
        gravity = PointMassGravity(no_tex_body, ())
        forces = ForceModel(gravity)
        integ = IntegratorConfig(DP8(); abstol=1e-12, reltol=1e-12, dt=60.0)
        prop = OrbitPropagator(forces, integ)
        propagate(prop, sat, StopAt(sat, PropDurationSeconds(), 3000.0))
        
        view = View3D(coord_sys=coord_sys)
        add_spacecraft!(view, sat)
        display_view(view; size=(800, 600))
        
        @test_reference "../reference/body_gray_default.png" view._scene
    end
    
end

println("\n" * "="^70)
println("Visual Regression Tests Complete")
println("="^70)
println("\nIf this was the first run, reference images were created in:")
println("  test/reference/")
println("\nPlease review them visually and commit if they look correct.")
println("\nFuture runs will compare against these references.")
println("="^70)
