# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# What the Earth chain does outside the EOP table.
#
# **This file records a known limitation. It is not an endorsement of it.**
#
# Earth orientation parameters are measured, not modelled: UT1−UTC and polar
# motion come from an IERS table that ends where the observations and the
# short-term predictions end. Today that is 1973-01-01 to about 2027-08.
#
# Ask for an epoch outside that and the table's end values are held constant.
# You get a rotation that is orthonormal, has determinant 1, composes
# correctly, round-trips, and is wrong — by an amount that grows the further
# out you go. Every structural check in this suite passes on it.
#
# The two directions are not alike, and the difference is the point.
#
# **Before 1973** you are warned, though not about EOP: the leap-second table
# starts in 1972, so an earlier epoch emits "Leapsecond of date ... not
# available, returning 0" — four times per call. Noisy, but you find out.
#
# **After 2027-08 nothing is said at all.** That is the direction anyone doing
# mission design goes, and it is silent.
#
# Magnitude: UT1−UTC drifts by roughly half a second per year once the
# predictions run out, and a second of UT1 is 15 arcseconds of Earth rotation,
# so the error accumulates at something like 250 m/year at LEO. A propagation
# to 2035 is a few kilometres out in Earth-fixed coordinates and says nothing
# about it.
#
# **Recommendation, not yet implemented.** Warn once per session when an epoch
# falls outside the loaded table, naming the span and the epoch asked for. Not
# an error — extending past the table is legitimate for mission design, and
# throwing would break studies that are entirely reasonable. But it should be
# a decision the user makes knowingly rather than one they never learn about.
#
# The tests below pin the current behaviour so that implementing the warning
# fails this file and forces it to be updated deliberately.
# =============================================================================

using AstroFrames
using AstroUniverse
using LinearAlgebra
using Test

# Just inside the table at each end, and outside it in both directions. They
# are kept apart because only one direction is silent.
const _EOP_INSIDE = (2442000.0, 2451545.0, 2458849.5, 2461000.0)
const _EOP_BEFORE = (2440000.0,)                                   # 1968
const _EOP_AFTER  = (2462500.0, 2465000.0, 2470000.0, 2488000.0)   # 2030-2100
const _EOP_OUTSIDE = (_EOP_BEFORE..., _EOP_AFTER...)

@testset "inside the table, polar motion varies with epoch" begin
    # The control for the test below. Polar motion is what TIRS → ITRF applies
    # and it is pure EOP, so inside the table this edge moves from year to
    # year. If this stopped being true the next test would pass for the wrong
    # reason.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())
        rotations = [Matrix(axes_rotation(TIRS(), ITRF(), jd))[1:3, 1:3]
                     for jd in _EOP_INSIDE]
        # Year to year, polar motion moves the matrix by 3e-7 to 6e-7.
        for i in 2:length(rotations)
            @test maximum(abs.(rotations[i] - rotations[1])) > 1e-7
        end
    finally
        set_frame_theory!(original)
    end
end

@testset "outside the table, the answer is silently frozen" begin
    # Two epochs a decade apart, both past the end. Polar motion is identical
    # between them, because both read the table's last row. Nothing warns.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())

        # Fifty years apart, both past the end. Measured 1.1e-10, against the
        # 3e-7 to 6e-7 that a real half-century of polar motion produces — so
        # the table has stopped contributing.
        #
        # Not exactly zero, because this edge also carries the TIO locator s',
        # which is a polynomial in time and keeps moving after the measured
        # data runs out. That residual is what is left here.
        far = Matrix(axes_rotation(TIRS(), ITRF(), 2470000.0))[1:3, 1:3]
        further = Matrix(axes_rotation(TIRS(), ITRF(), 2488000.0))[1:3, 1:3]
        @test maximum(abs.(far - further)) < 1e-9

        # And it looks perfectly healthy, which is the problem.
        #
        # The bound is 1e-12 rather than machine epsilon because orthonormality
        # decays slowly with distance from J2000 — measured 2.5e-16 near the
        # table, 1.5e-13 at year 2100, from the Earth-rotation angle being
        # reduced from a quantity that has grown to ~100 turns. That is 0.03
        # microarcseconds, four orders below anything the frames claim, and it
        # is unrelated to the EOP fault this file is about.
        for jd in _EOP_OUTSIDE
            M = Matrix(axes_rotation(ICRF(), ITRF(), jd))
            R = M[1:3, 1:3]
            @test norm(R' * R - I) < 1e-12
            @test isapprox(det(R), 1; atol = 1e-12)
            @test all(isfinite, M)

            # It round-trips too. No invariant in this suite can see the fault.
            B = Matrix(axes_rotation(ITRF(), ICRF(), jd))
            @test norm(Matrix(B * M) - I) < 1e-11
        end
    finally
        set_frame_theory!(original)
    end
end

@testset "the table covers the span we claim it does" begin
    # Pins the documented range. Polar motion still varying just inside each
    # end, and frozen outside, locates the boundary without reaching into the
    # table's internals — which belong to SatelliteToolboxTransformations and
    # are not ours to depend on.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())
        # 1e-7 sits between the two populations: real polar motion moves this
        # edge by 3e-7 or more, the frozen tail by 1.1e-10 at worst.
        moves(a, b) = maximum(abs.(
            Matrix(axes_rotation(TIRS(), ITRF(), a))[1:3, 1:3] -
            Matrix(axes_rotation(TIRS(), ITRF(), b))[1:3, 1:3])) > 1e-7

        # 1973-01-02 against 1974: inside, so it varies.
        @test moves(2441685.5, 2442050.0)
        # 2027-06 against 2027-08: inside, so it varies.
        @test moves(2461550.0, 2461640.0)

        # 1970 against 1972: both before the table, so frozen.
        @test !moves(2440587.5, 2441317.5)
        # 2028 against 2029: both after, so frozen.
        @test !moves(2461800.0, 2462200.0)
    finally
        set_frame_theory!(original)
    end
end

@testset "past the end of the table, nothing is said" begin
    # Recording today's behaviour explicitly. When the warning recommended in
    # this file's header is implemented, this test fails — which is the point.
    # Update it then, and update the header with it.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())
        axes_rotation(ICRF(), ITRF(), 2451545.0)   # drain any once-per-pair notice

        for jd in _EOP_AFTER
            logs, _ = Test.collect_test_logs() do
                axes_rotation(ICRF(), ITRF(), jd)
            end
            @test isempty(logs)
        end
    finally
        set_frame_theory!(original)
    end
end

@testset "before the table, the leap-second table does speak up" begin
    # Not an EOP warning — it comes from the leap-second table, which starts in
    # 1972 — but it does mean the past direction is not silent. Recorded so the
    # asymmetry with the future direction stays visible.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())
        for jd in _EOP_BEFORE
            logs, _ = Test.collect_test_logs() do
                axes_rotation(ICRF(), ITRF(), jd)
            end
            @test !isempty(logs)
            @test any(l -> occursin("leap-second data before 1972", string(l.message)), logs)
        end
    finally
        set_frame_theory!(original)
    end
end
