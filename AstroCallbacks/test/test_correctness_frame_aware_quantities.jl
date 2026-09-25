# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# The frame-aware quantity readers.
#
# These are what a user reaches for most, and they were exported untested. Two
# properties matter and neither is obvious from the code: the frame argument
# changes the answer, and it changes it by the right amount.
#
# They take any subject, so every check runs twice — once through a
# `Spacecraft` and once through an `AstroFrames.Coordinate`, which is the same
# call with no vehicle involved.
# =============================================================================

using Test
using LinearAlgebra
using AstroCallbacks
using AstroFrames
using AstroModels
using AstroStates
using AstroEpochs
using AstroUniverse

const _EPOCH_FQ = Time(2458849.5, 0.0, :tdb, :jd)
const _VEC_FQ   = [7000.0, 1200.0, 400.0, -1.3, 6.9, 2.1]
const _EQ_FQ    = CoordinateSystem(earth, MJ2000Eq())

_sc_fq() = Spacecraft(state = CartesianState(_VEC_FQ), time = _EPOCH_FQ, coord_sys = _EQ_FQ)
_co_fq() = Coordinate(_VEC_FQ, _EQ_FQ, _EPOCH_FQ)

@testset "position_vector" begin

    @testset "reads the state back in its own frame" begin
        for subject in (_sc_fq(), _co_fq())
            @test position_vector(subject, _EQ_FQ) ≈ _VEC_FQ[1:3] rtol = 1e-12
        end
    end

    @testset "the frame argument does real work" begin
        # A rotation preserves length and changes components. If the argument
        # were ignored these would be equal, and every non-native-frame caller
        # would be silently wrong.
        for subject in (_sc_fq(), _co_fq())
            r_eq   = position_vector(subject, _EQ_FQ)
            r_itrf = position_vector(subject, CoordinateSystem(earth, ITRF()))
            @test norm(r_itrf) ≈ norm(r_eq) rtol = 1e-12
            @test !isapprox(r_itrf, r_eq; rtol = 1e-6)
        end
    end

    @testset "an origin change moves it by the origin separation" begin
        for subject in (_sc_fq(), _co_fq())
            r_moon = position_vector(subject, CoordinateSystem(moon, MJ2000Eq()))
            @test 350_000 < norm(r_moon) < 420_000
        end
    end

    @testset "a spacecraft and a bare coordinate agree exactly" begin
        # The claim the subject interface makes: same call, same answer, one
        # implementation. A drift here means two code paths exist.
        for cs in (_EQ_FQ, CoordinateSystem(earth, ITRF()),
                   CoordinateSystem(earth, MJ2000Ec()), CoordinateSystem(moon, MJ2000Eq()))
            @test position_vector(_sc_fq(), cs) ≈ position_vector(_co_fq(), cs) rtol = 1e-14
        end
    end
end

@testset "velocity_vector" begin

    @testset "reads the state back in its own frame" begin
        for subject in (_sc_fq(), _co_fq())
            @test velocity_vector(subject, _EQ_FQ) ≈ _VEC_FQ[4:6] rtol = 1e-12
        end
    end

    @testset "a rotating target frame contributes omega cross r" begin
        # The term that makes velocity frame-dependent in a stronger sense than
        # position. At 7100 km and Earth's rate it is about 0.5 km/s, so an
        # implementation that rotated velocity without it would land far
        # outside this band.
        for subject in (_sc_fq(), _co_fq())
            v_eq   = velocity_vector(subject, _EQ_FQ)
            v_itrf = velocity_vector(subject, CoordinateSystem(earth, ITRF()))
            @test 0.2 < abs(norm(v_itrf) - norm(v_eq)) < 0.8
        end
    end

    @testset "a spacecraft and a bare coordinate agree exactly" begin
        for cs in (_EQ_FQ, CoordinateSystem(earth, ITRF()))
            @test velocity_vector(_sc_fq(), cs) ≈ velocity_vector(_co_fq(), cs) rtol = 1e-14
        end
    end
end

@testset "raan" begin

    @testset "the number changes with the frame, which is the whole point" begin
        # RAAN is measured from the frame's X axis in the frame's equatorial
        # plane. Reading it without naming a frame is how it ends up referred
        # to the wrong one, and nothing about the number says so.
        for subject in (_sc_fq(), _co_fq())
            Ω_eq   = raan(subject, _EQ_FQ)
            Ω_tod  = raan(subject, CoordinateSystem(earth, TODEq()))
            Ω_itrf = raan(subject, CoordinateSystem(earth, ITRF()))

            for Ω in (Ω_eq, Ω_tod, Ω_itrf)
                @test isfinite(Ω)
                @test 0 <= Ω <= 2π
            end

            # Precession over two decades is degrees, not arcseconds.
            @test rad2deg(abs(Ω_eq - Ω_tod)) > 0.1
            # Earth rotation puts the terrestrial value somewhere else entirely.
            @test !isapprox(Ω_itrf, Ω_eq; atol = deg2rad(1.0))
        end
    end

    @testset "a spacecraft and a bare coordinate agree exactly" begin
        for cs in (_EQ_FQ, CoordinateSystem(earth, TODEq()), CoordinateSystem(earth, ITRF()))
            @test raan(_sc_fq(), cs) ≈ raan(_co_fq(), cs) rtol = 1e-12
        end
    end
end

@testset "orbit-relative frames take their reference through params" begin
    # A frame defined by another object's orbit cannot be resolved downward, so
    # the reader has to carry `params` through as well.
    chief = [7100.0, 1200.0, 400.0, -1.28, 6.95, 2.08]
    ric   = CoordinateSystem(earth, RIC())

    for subject in (_sc_fq(), _co_fq())
        r = position_vector(subject, ric, (; reference_state = chief))
        @test length(r) == 3
        @test all(isfinite, r)
        # The two orbits are ~100 km apart, and RIC is centred on the same
        # origin, so the magnitude is the deputy's radius, not the separation.
        @test 7000 < norm(r) < 7300
    end

    @testset "and a missing reference fails by name" begin
        e = try
            position_vector(_co_fq(), ric)
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("reference_state", e.msg)
    end
end
