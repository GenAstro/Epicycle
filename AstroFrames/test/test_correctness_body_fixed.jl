# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Body-fixed axes: the tag, and the rotation it produces.
#
# `I-1`/`I-2`/`I-3` in `test_correctness_invariants.jl` already check block
# form, orthonormality and round trip across every body. This file covers what
# those cannot: that the tag resolves to the right body, that the rotation
# turns at the rate the body actually turns at, and that the tag prints
# usefully.
# =============================================================================

using AstroFrames
using AstroUniverse
using LinearAlgebra
using Test

const _JD_BF = 2458849.5

"""Rotation period implied by the rate block of a state transform, in hours."""
function _period_hours(M)
    R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]
    Ω = Ṙ * R'
    ω = sqrt(Ω[3,2]^2 + Ω[1,3]^2 + Ω[2,1]^2)
    return 2π / ω / 3600
end

@testset "body-fixed axes resolve to a body" begin

    @testset "the sentinel takes the body from the origin" begin
        # `CelestialBodyFixed()` carries no body; the coordinate system fills
        # it in. This is the form a user writes, so it is the form under test.
        for (body, naif) in ((sun, 10), (mercury, 199), (venus, 299), (mars, 499),
                             (jupiter, 599), (saturn, 699), (uranus, 799),
                             (neptune, 899), (pluto, 999))
            cs = CoordinateSystem(body, CelestialBodyFixed())
            @test cs.axes === CelestialBodyFixed{naif}()
            @test AstroFrames.naifid(cs.axes) == naif
        end
    end

    @testset "an explicit body agrees with the sentinel" begin
        @test CoordinateSystem(mars, CelestialBodyFixed()).axes ===
              CelestialBodyFixed{499}()
    end

    @testset "the tag says which body it is" begin
        # A user reading a coordinate system back has to be able to tell.
        @test sprint(show, CelestialBodyFixed{499}()) == "CelestialBodyFixed(Mars)"
        @test sprint(show, CelestialBodyFixed{10}())  == "CelestialBodyFixed(Sun)"

        # A body outside the shipped table still identifies itself, by number.
        @test sprint(show, CelestialBodyFixed{401}()) == "CelestialBodyFixed(NAIF-401)"

        # And the unresolved sentinel does not pretend to be a body.
        @test sprint(show, CelestialBodyFixed()) == "CelestialBodyFixed(unresolved)"

        @test sprint(show, MIME"text/plain"(), CelestialBodyFixed{499}()) ==
              "CelestialBodyFixed(Mars)"
    end
end

@testset "body-fixed axes turn at the body's rate" begin
    # The strongest check available without a truth row: the rate block implies
    # a rotation period, and every one of these is independently published to
    # four or more figures. It tests the pole construction, the chain rule in
    # `pole_axes_rotation`, and the orientation constants together — and it is
    # what exposed the Neptune convention trap below.
    #
    # Values are the commonly published sidereal rotation periods. Tolerance is
    # loose (0.1%) because the sources round and some are System III
    # conventions; a transcription error shows up far outside it.

    @testset "Sun and the inner planets" begin
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{10}(),  _JD_BF)) ≈ 25.38 * 24 rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{199}(), _JD_BF)) ≈ 58.646 * 24 rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{299}(), _JD_BF)) ≈ 243.025 * 24 rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{499}(), _JD_BF)) ≈ 24.6229 rtol = 1e-3
    end

    @testset "the giants and Pluto" begin
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{599}(), _JD_BF)) ≈ 9.9250 rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{699}(), _JD_BF)) ≈ 10.656 rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{799}(), _JD_BF)) ≈ 17.24  rtol = 1e-3
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{999}(), _JD_BF)) ≈ 6.3872 * 24 rtol = 1e-3
    end

    @testset "Neptune, whose published period is not the IAU rotation rate" begin
        # Worth its own testset because comparing Neptune against the commonly
        # quoted 16.11 h makes correct code look wrong — which is exactly what
        # happened when this file was written.
        #
        # 16.11 h is the rotation of Neptune's MAGNETIC FIELD, from Voyager 2
        # radio measurements, and it is what most general sources still give.
        # The IAU rotational element is
        #
        #     W = 249.978 + 541.1397757 d − 0.48 sin N
        #
        # (Archinal et al., 2018, Table 1), and 541.1397757°/day is 15.967 h —
        # the atmospheric determination that superseded the magnetic one. They
        # differ by 0.9%, which is far outside any tolerance and reads as a
        # transcription error if you do not know the two conventions exist.
        #
        # So this checks against the IAU constant itself rather than against a
        # published period, which is the only defensible reference for a value
        # the IAU defines.
        expected = 360 / 541.1397757 * 24
        @test _period_hours(axes_rotation(ICRF(), CelestialBodyFixed{899}(), _JD_BF)) ≈ expected rtol = 1e-4

        # And it is emphatically not the magnetic period, so nobody "corrects"
        # it back later.
        @test !isapprox(_period_hours(axes_rotation(ICRF(), CelestialBodyFixed{899}(), _JD_BF)),
                        16.11; rtol = 1e-3)
    end

    @testset "the pole is where the orientation model puts it" begin
        # Row 3 of R is the body's north pole in ICRF. Checked against the
        # model directly, so a mix-up between bodies cannot pass.
        for naif in (10, 199, 299, 499, 599, 699, 799, 899, 999)
            R = axes_rotation(ICRF(), CelestialBodyFixed{naif}(), _JD_BF)[1:3, 1:3]
            o = iau2015_orientation(naif, _JD_BF)
            pole = [cos(o.dec_pole)cos(o.ra_pole), cos(o.dec_pole)sin(o.ra_pole), sin(o.dec_pole)]
            @test R[3, :] ≈ pole atol = 1e-12
        end
    end
end

@testset "body-fixed axes route" begin
    # Parametric over NAIF ID, so they cannot be enumerated in the edge list —
    # their routes come from `hub_axes`. Before that existed, every one of
    # these raised.

    @testset "reached from the Earth frames under either theory" begin
        original = frame_theory()
        try
            set_frame_theory!(FK5())
            for source in (MJ2000Eq(), MODEq(), TODEq(), ITRF())
                M = axes_rotation(source, CelestialBodyFixed{499}(), _JD_BF)
                @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
            end

            set_frame_theory!(IAU2006())
            for source in (GCRF(), CIRS(), TIRS(), ITRF())
                M = axes_rotation(source, CelestialBodyFixed{499}(), _JD_BF)
                @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
            end
        finally
            set_frame_theory!(original)
        end
    end

    @testset "the routed answer equals the hand-composed one" begin
        M_routed = axes_rotation(ITRF(), CelestialBodyFixed{499}(), _JD_BF)
        M_hand   = axes_rotation(ICRF(), CelestialBodyFixed{499}(), _JD_BF) *
                   axes_rotation(ITRF(), ICRF(), _JD_BF)
        @test M_routed ≈ M_hand
    end

    @testset "one body to another" begin
        M = axes_rotation(CelestialBodyFixed{499}(), CelestialBodyFixed{599}(), _JD_BF)
        @test M ≈ axes_rotation(ICRF(), CelestialBodyFixed{599}(), _JD_BF) *
                  axes_rotation(CelestialBodyFixed{499}(), ICRF(), _JD_BF)
        @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-13
    end
end
