# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# The rate block of the Earth chain.
#
# The ERFA truth rows validate `R` for every Earth edge to microarcseconds and
# say nothing about `Ṙ`, which every velocity transformation uses. This file
# covers it, by finite difference — an independent check, since the derivative
# is taken of the rotation rather than read from the same source.
#
# It sorts the edges into three kinds, because they are not equally correct and
# reporting them together would hide that.
#
# **Static.** The frame bias and the J2000 obliquity really are constant, so
# `Ṙ = 0` is exact and is asserted as exact.
#
# **Earth rotation.** `CIRS → TIRS` and `TODEq → PEF` carry the sidereal
# rotation and model its rate. These agree with a finite difference to 1.5e-6
# relative, which is the difference's own truncation at h = 10 s.
#
# **Modelled as static, but not.** Precession, nutation, polar motion and the
# obliquity of date all vary, and Epicycle differentiates none of them: those
# edges return `Ṙ = 0`. This matches GMAT and Vallado, which also differentiate
# only the sidereal rotation, and it is a deliberate approximation rather than
# an oversight — but it was undocumented and untested, so it could have changed
# in either direction without anyone noticing.
#
# The tests below pin it from both sides. They assert the rate is exactly zero
# *and* bound how much is being left out, so that implementing the true rates
# fails this file loudly and forces the documentation to be updated with it.
#
# Measured at 2020-01-01, as a fraction of Earth's rotation rate:
#
#     precession       1.06e-7        polar motion     1.43e-9
#     nutation         4.32e-8        obliquity rate   9.86e-10
#     precession-nutation (CIO)  4.25e-8
#
# Reaching the user as a velocity error in the composed ICRF → ITRF:
#
#     LEO  0.76 mm/s      GPS  2.88 mm/s      GEO  4.56 mm/s
#
# Small in absolute terms. Worth knowing at GEO, where the ITRF velocity is
# itself near zero and this is a percent-level share of it.
# =============================================================================

using AstroFrames
using AstroUniverse
using LinearAlgebra
using Test

const _RATE_EPOCHS = (2451545.0007428695, 2455197.5007660175,
                      2458849.5008007395, 2460676.50080074)

const _OMEGA_EARTH = 7.292115e-5      # rad/s, IERS nominal

"""Central difference of the rotation block, in per-second units."""
function _numerical_rate(source, target, jd; h_seconds = 10.0)
    h_days = h_seconds / 86_400
    forward  = Matrix(axes_rotation(source, target, jd + h_days))[1:3, 1:3]
    backward = Matrix(axes_rotation(source, target, jd - h_days))[1:3, 1:3]
    return (forward - backward) ./ (2 * h_seconds)
end

"""Magnitude of the angular velocity implied by a rotation and its rate."""
_spin(R, Ṙ) = (Ω = Ṙ * R'; norm((Ω[3, 2], Ω[1, 3], Ω[2, 1])))

"""Run `f` under `theory`, then put the original theory back."""
function _with_theory(f, theory)
    original = frame_theory()
    try
        set_frame_theory!(theory)
        f()
    finally
        set_frame_theory!(original)
    end
end

@testset "static edges have exactly no rate" begin
    # Not "nearly zero" — these rotations do not depend on the epoch at all, so
    # anything nonzero here is a bug rather than a small error.
    _with_theory(FK5()) do
        for jd in _RATE_EPOCHS
            for (source, target) in ((ICRF(), MJ2000Eq()),        # frame bias
                                     (MJ2000Eq(), MJ2000Ec()))    # J2000 obliquity
                M = Matrix(axes_rotation(source, target, jd))
                @test all(iszero, M[4:6, 1:3])
                @test _numerical_rate(source, target, jd) == zeros(3, 3)
            end
        end
    end
end

@testset "the Earth-rotation edges carry a correct rate" begin
    # The only edges that model a rate. Checked against a finite difference,
    # which is independent of how the rate is derived. Measured 1.25e-6 to
    # 1.5e-6; the bound is the difference's truncation with headroom, not a
    # tolerance for a modelling error.
    for (theory, source, target) in ((IAU2006(), CIRS(), TIRS()),
                                     (FK5(),     TODEq(), PEF()))
        _with_theory(theory) do
            for jd in _RATE_EPOCHS
                M = Matrix(axes_rotation(source, target, jd))
                analytic = M[4:6, 1:3]
                numerical = _numerical_rate(source, target, jd)

                @test norm(numerical - analytic) / norm(analytic) < 1e-5

                # And it is Earth's rotation, not some other rate of the right
                # size — this is what catches a sidereal/solar day mix-up.
                @test abs(_spin(M[1:3, 1:3], analytic) / _OMEGA_EARTH - 1) < 1e-6
            end
        end
    end
end

@testset "the composed chain turns at Earth's rate" begin
    # End to end, both theories. Agreement with the IERS nominal is 1.5e-8,
    # which is UT1 drift rather than error.
    for theory in (IAU2006(), FK5())
        _with_theory(theory) do
            for jd in _RATE_EPOCHS
                M = Matrix(axes_rotation(ICRF(), ITRF(), jd))
                @test abs(_spin(M[1:3, 1:3], M[4:6, 1:3]) / _OMEGA_EARTH - 1) < 1e-7
            end
        end
    end
end

@testset "precession, nutation and polar motion are modelled as static" begin
    # A deliberate approximation, asserted from both sides.
    #
    # The first assertion pins the behaviour: these edges return no rate. The
    # second bounds what that leaves out, so the approximation stays small even
    # if the underlying models change.
    #
    # If someone implements the true rates, the first assertion fails. That is
    # the intent — the change should be deliberate and should carry a
    # documentation update with it, not slip in.
    edges = ((IAU2006(), GCRF(),     CIRS(),  "precession-nutation"),
             (IAU2006(), TIRS(),     ITRF(),  "polar motion"),
             (FK5(),     MJ2000Eq(), MODEq(), "precession"),
             (FK5(),     MODEq(),    TODEq(), "nutation"),
             (FK5(),     PEF(),      ITRF(),  "polar motion"),
             (FK5(),     MODEq(),    MODEc(), "obliquity of date"))

    for (theory, source, target, what) in edges
        _with_theory(theory) do
            @testset "$(what)" begin
                for jd in _RATE_EPOCHS
                    M = Matrix(axes_rotation(source, target, jd))
                    @test all(iszero, M[4:6, 1:3])

                    # What is being neglected, as a share of Earth's rotation.
                    # Largest is precession at 1.06e-7.
                    R = M[1:3, 1:3]
                    @test _spin(R, _numerical_rate(source, target, jd)) / _OMEGA_EARTH < 2e-7
                end
            end
        end
    end
end

@testset "the velocity this costs a user is bounded" begin
    # The number that actually matters: transform a state to ITRF with our
    # rate block, and with a fully numerical one, and difference the velocity.
    #
    # Measured 0.76 mm/s at LEO, 2.9 at GPS, 4.6 at GEO. The bound is 1 cm/s,
    # which is loose enough to survive an EOP update and tight enough that
    # losing the sidereal rate entirely — the failure that matters — trips it
    # by four orders of magnitude.
    _with_theory(IAU2006()) do
        for jd in _RATE_EPOCHS
            M = Matrix(axes_rotation(ICRF(), ITRF(), jd))
            numerical = _numerical_rate(ICRF(), ITRF(), jd)

            for (radius, speed) in ((7000.0, 7.546), (26_600.0, 3.873), (42_164.0, 3.075))
                r = [radius, 0.0, 0.0]
                v = [0.0, speed, 0.0]
                ours = M[4:6, 1:3] * r + M[1:3, 1:3] * v
                full = numerical * r + M[1:3, 1:3] * v
                @test norm(ours - full) < 1e-5      # km/s, i.e. 1 cm/s
            end
        end
    end
end

@testset "the rate block keeps its shape" begin
    # [R 0; Ṙ R]. The lower-right block being R rather than another copy of Ṙ
    # is what makes the matrix map velocity correctly.
    for theory in (IAU2006(), FK5())
        _with_theory(theory) do
            for jd in _RATE_EPOCHS
                for (source, target) in ((ICRF(), ITRF()), (GCRF(), TIRS()),
                                         (MJ2000Eq(), PEF()), (ICRF(), MJ2000Ec()))
                    M = Matrix(axes_rotation(source, target, jd))
                    @test size(M) == (6, 6)
                    @test all(iszero, M[1:3, 4:6])
                    @test M[4:6, 4:6] == M[1:3, 1:3]

                    # The reverse edge transposes both blocks.
                    B = Matrix(axes_rotation(target, source, jd))
                    @test norm(B[1:3, 1:3] - M[1:3, 1:3]') < 1e-14
                    @test norm(B[4:6, 1:3] - M[4:6, 1:3]') < 1e-14
                end
            end
        end
    end
end
