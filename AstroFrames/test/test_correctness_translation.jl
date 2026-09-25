# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Origin translation and route composition.
#
# The second primitive, and the routing that a translation immediately needs:
# the ephemeris arrives in ICRF, so expressing an offset in any other axes
# requires composing the chain.
#
# Requires SPICE kernels (an ephemeris) and IERS EOP.
# =============================================================================

using AstroFrames
using AstroUniverse
using AstroEpochs: Time
using LinearAlgebra
using Test

const _JD = 2458849.5

@testset "origin translation" begin

    @testset "Earth ↔ Moon offset is physically right" begin
        Δ = origin_translation(earth, moon, ICRF(), _JD)
        @test length(Δ) == 6
        # The lunar orbit spans roughly 356 000–407 000 km, and the Moon moves
        # at about 1 km/s. Bounds rather than a fixed value: this is a check
        # that the ephemeris is being read correctly, not a truth row.
        @test 350_000 < norm(Δ[1:3]) < 410_000
        @test 0.9 < norm(Δ[4:6]) < 1.1
    end

    @testset "direction: source origin relative to target" begin
        # The contract is the state of the SOURCE origin seen from the TARGET.
        # Getting it backwards is a sign error of full magnitude that no
        # structural check would catch, so it is pinned against the raw
        # ephemeris, whose own direction is documented the other way round.
        Δ  = origin_translation(earth, moon, ICRF(), _JD)
        mre = translate_state(earth, moon, _JD)      # Moon seen from Earth
        @test all(isapprox.(Δ, -mre; rtol = 1e-12))

        # Stated as the degenerate case that fixes the convention: a point at
        # the source origin has zero coordinates there, so about the target
        # origin its coordinates are exactly the offset.
        state_at_earth_centre = zeros(6)
        @test all(isapprox.(axes_rotation(ICRF(), ICRF(), _JD) * state_at_earth_centre .+ Δ,
                            Δ; rtol = 1e-12))
    end

    @testset "both epoch forms are supported and agree" begin
        # The `Time` overload is the one a user reaches for, because a `Time`
        # carries its own scale and cannot be passed as the wrong one.
        Δ_jd   = origin_translation(earth, moon, ICRF(), _JD)
        Δ_time = origin_translation(earth, moon, ICRF(), Time(_JD, 0.0, :tdb, :jd))
        @test all(isapprox.(Δ_jd, Δ_time; rtol = 1e-14))

        # And the scale is honoured rather than ignored: the same number read
        # as UTC is a different instant, so the Moon is somewhere else.
        Δ_utc = origin_translation(earth, moon, ICRF(), Time(_JD, 0.0, :utc, :jd))
        @test !all(isapprox.(Δ_jd, Δ_utc; rtol = 1e-12))
    end

    @testset "the axes argument does real work" begin
        Δ_icrf = origin_translation(earth, moon, ICRF(), _JD)
        Δ_itrf = origin_translation(earth, moon, ITRF(), _JD)

        # A rotation preserves length but changes components. If the axes
        # argument were ignored these would be equal, and the offset would be
        # silently wrong by a rotation for every non-ICRF caller.
        @test isapprox(norm(Δ_itrf[1:3]), norm(Δ_icrf[1:3]); rtol = 1e-12)
        @test !all(isapprox.(Δ_itrf, Δ_icrf; rtol = 1e-9))
    end

    @testset "the frame bias is applied for FK5 axes" begin
        # The 16 km-at-1-AU trap, at lunar range. Asking for the offset in an
        # FK5 frame must rotate through the bias rather than hand back the
        # ICRF-aligned ephemeris unchanged.
        Δ_icrf = origin_translation(earth, moon, ICRF(), _JD)
        Δ_fk5  = origin_translation(earth, moon, MJ2000Eq(), _JD)
        Δkm    = norm(Δ_fk5[1:3] - Δ_icrf[1:3])

        # 23.15 mas over ~404 000 km is ~45 m. Bounded on both sides: zero
        # would mean the bias was skipped, large would mean it was misapplied.
        @test 0.02 < Δkm < 0.10
    end
end

@testset "route composition" begin

    @testset "a routed pair equals the hand-composed chain" begin
        original = frame_theory()
        try
            set_frame_theory!(IAU2006())
            @test axes_rotation(ICRF(), ITRF(), _JD) ≈
                  axes_rotation(TIRS(), ITRF(), _JD) *
                  axes_rotation(CIRS(), TIRS(), _JD) *
                  axes_rotation(GCRF(), CIRS(), _JD) *
                  axes_rotation(ICRF(), GCRF(), _JD)

            set_frame_theory!(FK5())
            @test axes_rotation(ICRF(), ITRF(), _JD) ≈
                  axes_rotation(PEF(), ITRF(), _JD) *
                  axes_rotation(TODEq(), PEF(), _JD) *
                  axes_rotation(MODEq(), TODEq(), _JD) *
                  axes_rotation(MJ2000Eq(), MODEq(), _JD) *
                  axes_rotation(ICRF(), MJ2000Eq(), _JD)
        finally
            set_frame_theory!(original)
        end
    end

    @testset "the theory selects the route" begin
        original = frame_theory()
        try
            set_frame_theory!(IAU2006()); M6 = axes_rotation(ICRF(), ITRF(), _JD)
            set_frame_theory!(FK5());     M5 = axes_rotation(ICRF(), ITRF(), _JD)

            # Same endpoints, different chain, therefore different numbers —
            # which is the theory setting doing its job, not an error.
            @test !(M5 ≈ M6)

            # And the gap is the frame-bias scale: the FK5 chain refers its
            # precession to MJ2000Eq, itself ~23 mas off ICRF.
            R = M5[1:3,1:3] * M6[1:3,1:3]'
            v = 0.5 .* (R[3,2]-R[2,3], R[1,3]-R[3,1], R[2,1]-R[1,2])
            @test 20.0 < rad2deg(asin(sqrt(sum(abs2, v)))) * 3600 * 1000 < 30.0
        finally
            set_frame_theory!(original)
        end
    end

    @testset "routes round-trip and stay orthonormal" begin
        for theory in (FK5(), IAU2006())
            original = frame_theory()
            try
                set_frame_theory!(theory)
                M = axes_rotation(ICRF(), ITRF(), _JD)
                @test norm(M * axes_rotation(ITRF(), ICRF(), _JD) - I) < 1e-14
                @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-14
                # The spin edge is inside the route, so the composed rate is
                # still Earth's — the route must not lose it.
                Ω = -M[4:6,1:3] * M[1:3,1:3]'
                @test Ω[2,1] ≈ 7.292115146706979e-5 rtol = 1e-7
            finally
                set_frame_theory!(original)
            end
        end
    end

    @testset "an undeclared pair still fails loudly" begin
        # Routing must not turn a genuinely unsupported pair into a silent
        # wrong answer. `Inertial` has no edges, so no route exists under any
        # theory and the loud error stands.
        #
        # MODEq → CIRS is deliberately NOT the example here: it is reachable
        # across the two chains, so under FR-FRAME-4 it warns and proceeds.
        err = try
            axes_rotation(Inertial(), ITRF(), _JD)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Inertial", err.msg)
    end
end
