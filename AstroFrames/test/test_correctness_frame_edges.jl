# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# The cases found in the 2026-09-30 production review, one testset each.

using AstroFrames
using AstroUniverse
using AstroEpochs
using LinearAlgebra
using ForwardDiff
using Test
import SatelliteToolboxTransformations as STB

@testset "Frame edge cases" begin

    @testset "Earth-fixed axes are right on a leap-second day" begin
        # AstroEpochs spreads a leap-second day over 86401 s; the IERS tables and UT1 assume an
        # ordinary day. Fed the stretched date, ICRF → ITRF was 474 m off at LEO by 23:00 on
        # 2016-12-31. Checked against SatelliteToolboxTransformations given the ordinary date.
        original = frame_theory()
        try
            set_frame_theory!(IAU2006())
            for (h, m, s) in ((12, 0, 0), (23, 0, 0), (23, 59, 59))
                t = Time("2016-12-31T$(lpad(h, 2, '0')):$(lpad(m, 2, '0')):$(lpad(s, 2, '0'))", UTC(), ISOT())
                R = Matrix(axes_rotation(ICRF(), ITRF(), t))[1:3, 1:3]
                Rs = Matrix(STB.r_eci_to_ecef(STB.DCM, STB.GCRF(), STB.ITRF(),
                                              STB.date_to_jd(2016, 12, 31, h, m, s), eop(IAU2006())))
                @test opnorm(R - Rs) * 6778 < 0.02e-3        # km; under 2 cm at LEO
            end
        finally
            set_frame_theory!(original)
        end
    end

    @testset "the Earth rotation keeps the two-part precision of Time" begin
        # A single Float64 Julian date resolves 40 µs; a 10 µs step must still turn the Earth.
        t0 = Time(2460400.0, 0.1, TDB(), JD())
        t1 = Time(2460400.0, 0.1 + 10e-6 / 86400, TDB(), JD())
        R0 = Matrix(axes_rotation(ICRF(), ITRF(), t0))[1:3, 1:3]
        R1 = Matrix(axes_rotation(ICRF(), ITRF(), t1))[1:3, 1:3]
        @test isapprox(opnorm(R1 - R0), 7.2921e-5 * 10e-6; rtol = 1e-2)
        # and the epoch derivative is still there
        g = ForwardDiff.derivative(x -> axes_rotation(ICRF(), ITRF(), Time(2460400.0, x, TDB(), JD()))[1, 2], 0.1)
        @test isfinite(g) && g != 0
    end

    @testset "an EOP read allocates little" begin
        e = AstroFrames._scales(Time("2024-03-01T07:13:00", UTC(), ISOT()))
        axes_rotation(ICRF(), ITRF(), e)
        # 13.6 KB before the field names were made compile-time
        @test (@allocated axes_rotation(ICRF(), ITRF(), e)) < 2000
    end

    @testset "rotations between fixed celestial axes do not convert the epoch" begin
        t = Time("2024-03-01T07:13:00", UTC(), ISOT())
        axes_rotation(ICRF(), MJ2000Eq(), t)
        @test (@allocated axes_rotation(ICRF(), MJ2000Eq(), t)) == 0
        @test axes_rotation(ICRF(), MJ2000Ec(), t) == axes_rotation(ICRF(), MJ2000Ec(), 2451545.0)
        # An orbit-relative frame with a dummy epoch no longer warns about pre-1972 leap seconds
        x = [7000.0, 0, 0, 0, 7.5, 0]
        @test_logs axes_rotation(ICRF(), RIC(), 0.0, (; reference_state = x))
    end

    @testset "coordinate systems are values" begin
        a = CoordinateSystem(earth, ICRF())
        b = CoordinateSystem(earth, ICRF())
        @test a === b
        @test a == CoordinateSystem(deepcopy(earth), ICRF())
        @test hash(a) == hash(CoordinateSystem(deepcopy(earth), ICRF()))
        @test a != CoordinateSystem(moon, ICRF())
        @test a != CoordinateSystem(earth, MJ2000Eq())
        # The exported frames cannot be changed for everyone
        @test_throws ErrorException (EarthICRF.origin = sun)
        # Coordinates built from the same numbers are equal
        t = Time("2024-03-01T07:13:00", UTC(), ISOT())
        @test Coordinate([7000.0, 0, 0, 0, 7.5, 0], a, t) == Coordinate([7000.0, 0, 0, 0, 7.5, 0], b, t)
    end

    @testset "orbit-relative parameters are checked" begin
        e = 2451545.0
        @test_throws ArgumentError axes_rotation(ICRF(), RIC(), e, (; reference_state = [7000.0, 0, 0]))
        @test_throws ArgumentError axes_rotation(ICRF(), VNB(), e,
            (; reference_state = [7000.0, 0, 0, 0, 7.5, 0], reference_accel = [0.0, 0.0]))
        err = try
            axes_rotation(ICRF(), VNB(), e, (; reference_state = [7000.0, 0, 0, 0, 7.5, 0],
                                              reference_acceleration = [0.0, 0.0, 1e-6]))
        catch x
            x
        end
        @test err isa ArgumentError && occursin("reference_acceleration", err.msg)
        # A missing reference names the calls a user writes
        err = try axes_rotation(ICRF(), RIC(), e) catch x; x end
        @test err isa ArgumentError && occursin("CoordinateSystem(chief, RIC())", err.msg)
        # A reference with no plane warns, and the axes are NaN
        @test_logs (:warn, r"zero angular momentum") match_mode = :any begin
            M = axes_rotation(ICRF(), RIC(), e, (; reference_state = [7000.0, 0, 0, 1.0, 0, 0]))
            @test any(isnan, M)
        end
    end

    @testset "epochs of other number types" begin
        @test axes_rotation(ICRF(), ITRF(), Float32(2.4604e6)) ≈ axes_rotation(ICRF(), ITRF(), 2.4604e6)
        # A dual epoch cannot reach the ephemeris; the message says so
        err = try
            origin_translation(earth, moon, ICRF(), ForwardDiff.Dual(2451545.0, 1.0))
        catch x
            x
        end
        @test err isa ArgumentError && occursin("not differentiable", err.msg)
    end
end
