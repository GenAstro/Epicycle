# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# target!: a targeting problem written in flight order.
#
# The GEO transfer is the problem in runtests_sequence_geotransfer.jl, which checks the assembled
# event graph against GMAT. Here the same problem is written as a block and solved, and the solved
# burns must meet every target. The error cases check that a malformed block says what is wrong
# and leaves no recorder installed behind it.

using Test
using AstroProp: Vern9

@testset "target! block" begin

    function _geo_setup()
        sat = Spacecraft(
            state = CartesianState([3737.792, -4607.692, -2845.644, 5.411, 5.367, -1.566]),
            time  = Time("2000-01-01T11:59:28.000", UTC(), ISOT()),
            name  = "GeoSat-1")
        prop = OrbitPropagator(ForceModel(PointMassGravity(earth, ())),
                               IntegratorConfig(Vern9(); abstol = 1e-12, reltol = 1e-12, dt = 60.0))
        toi = ImpulsiveManeuver(axes = VNB(), element1 = 2.518, element2 = 0.0,   element3 = 0.0)
        mcc = ImpulsiveManeuver(axes = VNB(), element1 = 0.559, element2 = 0.588, element3 = 0.0)
        moi = ImpulsiveManeuver(axes = VNB(), element1 = 0.282, element2 = 0.0,   element3 = 0.0)
        return sat, prop, toi, mcc, moi
    end

    @testset "verbs act immediately outside a block" begin
        sat, prop, toi, _, _ = _geo_setup()
        r0 = copy(to_posvel(sat))
        maneuver!(sat, toi)
        @test to_posvel(sat)[4:6] != r0[4:6]
        propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), 60.0))
        @test to_posvel(sat)[1:3] != r0[1:3]
        v = Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 0.0, 0.0])
        @test v isa SolverVariable
        @test Constraint(position_magnitude, sat; equals = 1.0) isa Constraint
    end

    @testset "malformed blocks" begin
        sat, prop, toi, _, _ = _geo_setup()
        stop = StopAt(sat, PropDurationSeconds(), 60.0)

        # A constraint checks the step before it, so it cannot come first.
        @test_throws ArgumentError target!() do
            Constraint(position_magnitude, sat; equals = 42000.0)
        end
        @test_throws ArgumentError target!(() -> nothing)
        @test_throws ArgumentError target!() do
            propagate!(prop, sat, stop)
            Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 0.0, 0.0])
        end
        @test_throws ArgumentError target!() do
            target!() do
                propagate!(prop, sat, stop)
            end
        end

        # Every failure above left the recorders off, so the verbs act again.
        @test AstroSolve._BLOCK[] === nothing
        @test AstroProp._RECORDER[] === nothing
        @test AstroManeuvers._RECORDER[] === nothing
        r0 = copy(to_posvel(sat))
        propagate!(prop, sat, stop)
        @test to_posvel(sat) != r0
    end

    @testset "an exception inside the block clears the recorders" begin
        @test_throws ErrorException target!() do
            error("thrown by the block")
        end
        @test AstroSolve._BLOCK[] === nothing
        @test AstroProp._RECORDER[] === nothing
        @test AstroManeuvers._RECORDER[] === nothing
    end

    @testset "GEO transfer meets every target" begin
        sat, prop, toi, mcc, moi = _geo_setup()
        z_crossing = StopAt(position_z, sat, EarthMJ2000Eq; equals = 0.0)
        apoapsis   = StopAt(position_dot_velocity, sat; equals = 0.0, direction = -1)
        perigee    = StopAt(position_dot_velocity, sat; equals = 0.0, direction =  1)

        result = target!(method = Optimize(derivatives = :fd, print_level = 0)) do
            propagate!(prop, sat, z_crossing)

            Vary(delta_v, toi; lower_bound = [0.0, 0.0, 0.0], upper_bound = [8.0, 0.0, 0.0])
            maneuver!(sat, toi)

            propagate!(prop, sat, apoapsis)
            Constraint(position_magnitude, sat; equals = 85000.0)

            propagate!(prop, sat, perigee)
            propagate!(prop, sat, z_crossing)

            Vary(delta_v, mcc; lower_bound = [-1.0, -1.0, -0.001],
                               upper_bound = [ 4.0,  1.0,  0.001])
            maneuver!(sat, mcc)

            propagate!(prop, sat, perigee)
            Constraint(inclination, sat, EarthMJ2000Eq; equals = deg2rad(2.0))
            Constraint(position_magnitude, sat; equals = 42195.0)

            Vary(delta_v, moi; lower_bound = [-1.0, -0.001, -0.001],
                               upper_bound = [ 4.0,  0.001,  0.001])
            maneuver!(sat, moi)
            Constraint(semi_major_axis, sat; equals = 42166.90)
        end

        @test result.info in (:Solve_Succeeded, :Solved_To_Acceptable_Level)
        @test isapprox(result.constraints, [85000.0, deg2rad(2.0), 42195.0, 42166.90];
                       rtol = 1e-6)
        # The TOI the assembled graph converges to for the same problem.
        @test isapprox(delta_v(toi)[1], 2.819829; atol = 1e-5)
    end

end
nothing
