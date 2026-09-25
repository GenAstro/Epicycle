# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# RIC, LVLH and VNB against GMAT — external truth for the orbit-relative family.
#
# Everything else these frames are held to is analytic: orthonormality, the
# right rotation rate for a circular orbit, axes pointing where the docstring
# says. All of that was true of RIC while its second component was cross-track
# and its third pointed against the motion, which is what makes an independent
# implementation worth having.
#
# GMAT has no frames under these names, but `ObjectReferenced` expresses all
# three, and the mapping is stated here because it is the part that can be got
# wrong:
#
#     RIC    XAxis = R,  ZAxis = N     → Y = Z × X = N × R, in-track
#     VNB    XAxis = V,  YAxis = N     → Z = X × Y = V × N
#     LVLH   ZAxis = -R, YAxis = -N    → X = Y × Z = N × R
#
# The reported values are self-consistent in a way that confirms the mapping
# independently: LVLH X equals RIC Y, LVLH Z is the negative of RIC X, and
# VNB Y equals RIC Z.
#
# GMAT builds these about an unaccelerated reference, and so does Epicycle
# unless given an acceleration, so the two agree on all three frames including
# the rate block.
#
# **The comparison is exact, not approximate.** These axes are built from the
# reference state alone — not from an epoch, and not from whichever inertial
# frame that state is labelled with. Handing GMAT and Epicycle the same six
# numbers removes the frame bias, the time scale and the ephemeris from the
# comparison entirely. What is left is the rotation.
#
# The reference orbit is eccentric and inclined on purpose. On a circular orbit
# in-track coincides with velocity and these three frames are easy to confuse.
#
# Truth source: GMAT R2026a. The script and report are internal and not shipped.
# =============================================================================

using AstroFrames
using AstroUniverse
using LinearAlgebra
using Test

# The two states handed to both tools, km and km/s.
const _GMAT_REF = [8000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]
const _GMAT_TGT = [1000.0, -2000.0, 3000.0, 1.0, 2.0, 3.0]

# GMAT's Earth μ, used for the reference acceleration so the two tools are
# working from the same force model as well as the same state.
const _GMAT_MU = 398600.4415

# Reported by GMAT R2026a: the target expressed in each set of axes, so
# position is R·r and velocity is R·v + Ṙ·r.
const _GMAT_TRUTH = (
    ric  = [1444.63023702923, -1921.919212390072, 2866.926929537409,
            0.4736921610548745, 0.9386513912868919, 2.483391554348458],
    vnb  = [-1260.123838323872, 2866.926929537409, 2047.637148710952,
            2.594372608313854, 2.483391554348458, 1.049760523653619],
    lvlh = [-1921.919212390072, -2866.926929537409, -1444.63023702923,
            0.9386513912868919, -2.483391554348458, -0.4736921610548745],
)

_gmat_params() = (; reference_state = _GMAT_REF,
                    reference_accel = (-_GMAT_MU / norm(_GMAT_REF[1:3])^3) .* _GMAT_REF[1:3])

# Measured agreement is 2.3e-13 km on a 2900 km component — 8e-17 relative,
# which is double-precision round-off in GMAT's printed output.
const _TOL_GMAT = 1e-12

@testset "RIC and LVLH match GMAT exactly" begin
    # Position and velocity together, so the rate block is under test as well
    # as the rotation.
    for (axes, truth, name) in ((RIC(),  _GMAT_TRUTH.ric,  "RIC"),
                                (LVLH(), _GMAT_TRUTH.lvlh, "LVLH"))
        @testset "$(name)" begin
            ours = Matrix(axes_rotation(ICRF(), axes, 2458849.5, _gmat_params())) * _GMAT_TGT
            @test norm(ours[1:3] - truth[1:3]) / norm(truth[1:3]) < _TOL_GMAT
            @test norm(ours[4:6] - truth[4:6]) / norm(truth[4:6]) < _TOL_GMAT
        end
    end
end

@testset "GMAT's own numbers confirm the axis mapping" begin
    # Independent of Epicycle: if ObjectReferenced had been set up to mean
    # something other than intended, these identities would not hold, and the
    # comparison above would be against the wrong frame.
    ric, lvlh, vnb = _GMAT_TRUTH.ric, _GMAT_TRUTH.lvlh, _GMAT_TRUTH.vnb

    @test lvlh[1] ≈  ric[2]      # LVLH x̂ is the in-track direction
    @test lvlh[2] ≈ -ric[3]      # LVLH ŷ is −ĥ
    @test lvlh[3] ≈ -ric[1]      # LVLH ẑ is nadir
    @test vnb[2]  ≈  ric[3]      # VNB N̂ and RIC Ĉ are both the orbit normal
end

@testset "VNB matches GMAT, and both say the same thing" begin
    # Both tools build VNB about an unaccelerated reference unless told
    # otherwise, so with no acceleration supplied they agree exactly — rotation
    # and rate together.
    #
    # This is not GMAT getting the rate wrong. VNB's first axis follows the
    # velocity, which turns only when the velocity changes, so the frame's rate
    # is a statement about the force environment rather than about the
    # geometry. With no acceleration there is nothing to turn it, and a
    # non-rotating frame is the correct answer to the question asked — the one
    # you want when resolving a vector into along-track, normal and binormal
    # components at an instant, which is most of what VNB is for.
    unaccelerated = Matrix(axes_rotation(ICRF(), VNB(), 2458849.5,
                                         (; reference_state = _GMAT_REF))) * _GMAT_TGT
    truth = _GMAT_TRUTH.vnb

    @test norm(unaccelerated[1:3] - truth[1:3]) / norm(truth[1:3]) < _TOL_GMAT
    @test norm(unaccelerated[4:6] - truth[4:6]) / norm(truth[4:6]) < _TOL_GMAT

    # Where Epicycle goes further is that it will take a real acceleration and
    # give the frame's true rate. The axes are unchanged; only the rate moves,
    # and on this geometry it moves the velocity by 1.6 km/s.
    M = Matrix(axes_rotation(ICRF(), VNB(), 2458849.5, _gmat_params()))
    accelerated = M * _GMAT_TGT

    @test norm(accelerated[1:3] - truth[1:3]) / norm(truth[1:3]) < _TOL_GMAT
    @test norm(M[4:6, 1:3]) > 1e-4
    @test norm(accelerated[4:6] - truth[4:6]) > 1.0
end

@testset "RIC's in-track component is the one GMAT calls in-track" begin
    # The regression that motivated this file. Before the fix, RIC returned
    # (radial, cross-track, −in-track); component 2 below would have held
    # GMAT's component 3, and component 3 the negative of GMAT's component 2.
    ours = Matrix(axes_rotation(ICRF(), RIC(), 2458849.5, _gmat_params())) * _GMAT_TGT
    truth = _GMAT_TRUTH.ric

    @test sign(ours[2]) == sign(truth[2])
    @test !isapprox(ours[2], truth[3]; rtol = 1e-6)     # not the old swap
    @test !isapprox(ours[3], -truth[2]; rtol = 1e-6)    # not the old sign flip
end
