# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Orbit-relative axes, and the `params` mechanism they exercise.
#
# Needs no EOP and no kernels: the reference orbit is supplied by the caller,
# which is the whole point of `params`.
# =============================================================================

using AstroFrames
using LinearAlgebra
using Test

const _MU = 398600.4418        # km³/s², Earth
const _JD_OR = 2458849.5

"""Two-body state and acceleration for a circular orbit of radius `a`."""
function _circular(a)
    r = [a, 0.0, 0.0]
    v = [0.0, sqrt(_MU / a), 0.0]
    return vcat(r, v), (-_MU / a^3) .* r
end

"""An inclined, eccentric state, so no axis is accidentally aligned."""
function _generic()
    r = [-4550.0, 2220.0, 4980.0]
    v = [-3.10, -6.60, 0.12]
    a = (-_MU / norm(r)^3) .* r
    return vcat(r, v), a
end

"""Rotation rate vector and matrix read out of a transform block `[R 0; Rdot R]`."""
function _rate_vector(M)
    R, Rdot = M[1:3, 1:3], M[4:6, 1:3]
    W = Rdot * R'
    return [W[3,2], W[1,3], W[2,1]], W
end

@testset "orbit-relative axes (VNB)" begin

    @testset "for a circular orbit the frame turns at the orbital rate" begin
        # The definitive check. A circular orbit's velocity direction rotates
        # at exactly the mean motion, so the analytic Ṙ must reproduce
        # √(μ/a³). Nothing about orthonormality or a round trip would catch a
        # wrong rate; this does, and to full precision.
        for a in (7000.0, 26_600.0, 42_164.0)
            state, accel = _circular(a)
            p = (; reference_state = state, reference_accel = accel)
            M = axes_rotation(ICRF(), VNB(), _JD_OR, p)
            R, Ṙ = M[1:3, 1:3], M[4:6, 1:3]

            Ω = Ṙ * R'
            @test norm(Ω + Ω') < 1e-18              # a rotation rate is skew
            ω = sqrt(sum(abs2, (Ω[3,2], Ω[1,3], Ω[2,1])))
            @test ω ≈ sqrt(_MU / a^3) rtol = 1e-12
        end
    end

    @testset "geometry: V along velocity, N along r × v" begin
        state, accel = _circular(7000.0)
        p = (; reference_state = state, reference_accel = accel)
        R = axes_rotation(ICRF(), VNB(), _JD_OR, p)[1:3, 1:3]

        r = state[1:3]; v = state[4:6]
        @test R[1, :] ≈ v / norm(v)                       # V̂
        @test R[2, :] ≈ cross(r, v) / norm(cross(r, v))   # N̂
        @test R[3, :] ≈ cross(R[1, :], R[2, :])           # B̂ completes the set
        @test norm(R' * R - I) < 1e-15
        @test det(R) ≈ 1                                   # right-handed
    end

    @testset "round trip" begin
        state, accel = _circular(7000.0)
        p = (; reference_state = state, reference_accel = accel)
        @test norm(axes_rotation(ICRF(), VNB(), _JD_OR, p) *
                   axes_rotation(VNB(), ICRF(), _JD_OR, p) - I) < 1e-14
    end

    @testset "acceleration is used, not ignored" begin
        # If `reference_accel` were quietly dropped the frame would report as
        # non-rotating, which is the failure the requirement exists to prevent.
        # Two different accelerations must give two different rates.
        state, accel = _circular(7000.0)
        M_true = axes_rotation(ICRF(), VNB(), _JD_OR,
                               (; reference_state = state, reference_accel = accel))
        M_zero = axes_rotation(ICRF(), VNB(), _JD_OR,
                               (; reference_state = state, reference_accel = zeros(3)))
        @test M_true[1:3, 1:3] ≈ M_zero[1:3, 1:3]     # same orientation
        @test !(M_true[4:6, 1:3] ≈ M_zero[4:6, 1:3])  # different rate
    end

    @testset "a missing reference orbit fails by name, not as MethodError" begin
        for axes in (RIC(), LVLH(), VNB())
            e = try; axes_rotation(ICRF(), axes, _JD_OR); nothing; catch e; e; end
            @test e isa ArgumentError
            @test occursin("reference_state", e.msg)
            @test occursin("does not propagate", e.msg)
        end
    end

    @testset "acceleration is optional, and its absence means something" begin
        # VNB turns only under acceleration, so omitting it is not a missing
        # input — it states that the reference is not being accelerated, and
        # the frame that follows is exactly non-rotating. That is the frame you
        # want when resolving a vector into components at an instant, which is
        # most of what VNB is for.
        state, accel = _circular(7000.0)

        unaccelerated = axes_rotation(ICRF(), VNB(), _JD_OR, (; reference_state = state))
        accelerated   = axes_rotation(ICRF(), VNB(), _JD_OR,
                                      (; reference_state = state, reference_accel = accel))

        @test all(iszero, unaccelerated[4:6, 1:3])
        @test norm(accelerated[4:6, 1:3]) > 1e-6

        # Only the rate differs. The axes themselves do not depend on it.
        @test unaccelerated[1:3, 1:3] == accelerated[1:3, 1:3]

        # RIC and LVLH turn from the geometry, so they are unaffected either
        # way under two-body motion: their rate needs r × a, which is zero for
        # any radial acceleration.
        for axes in (RIC(), LVLH())
            without = axes_rotation(ICRF(), axes, _JD_OR, (; reference_state = state))
            with    = axes_rotation(ICRF(), axes, _JD_OR,
                                    (; reference_state = state, reference_accel = accel))
            @test without == with
            @test norm(without[4:6, 1:3]) > 1e-6      # and they do turn
        end
    end

    @testset "params are ignored by edges that do not need them" begin
        # So a caller can pass params uniformly without knowing which edges
        # along a route consume them.
        p = (; reference_state = _circular(7000.0)[1], reference_accel = zeros(3))
        @test axes_rotation(ICRF(), MJ2000Eq(), _JD_OR, p) ==
              axes_rotation(ICRF(), MJ2000Eq(), _JD_OR)
    end
end

@testset "orbit-relative axes (RIC)" begin

    @testset "geometry: rows are radial, in-track, cross-track in that order" begin
        state, _ = _generic()
        r = state[1:3]; v = state[4:6]
        R = axes_rotation(ICRF(), RIC(), _JD_OR, (; reference_state = state))[1:3, 1:3]

        R̂ = r / norm(r)
        Ĉ = cross(r, v) / norm(cross(r, v))
        Î = cross(Ĉ, R̂)

        @test R[1, :] ≈ R̂
        @test R[2, :] ≈ Î
        @test R[3, :] ≈ Ĉ
        @test norm(R' * R - I) < 1e-15
        @test det(R) ≈ 1

    end

    @testset "in-track is not the velocity direction off a circular orbit" begin
        # Pinning this keeps RIC from being quietly rebuilt as a velocity
        # frame, which is the standard way these two get conflated.
        #
        # It needs an orbit with a real flight-path angle to show. `_generic`
        # is nearly circular — 0.056°, so in-track and velocity agree to
        # 0.001 — which is why this uses its own state: radial velocity of
        # 1.5 km/s against 6.5 transverse, about 13°.
        r = [8000.0, 0.0, 0.0]
        v = [1.5, 6.5, 0.0]
        R = axes_rotation(ICRF(), RIC(), _JD_OR,
                          (; reference_state = vcat(r, v)))[1:3, 1:3]

        v̂ = v / norm(v)
        @test !isapprox(R[2, :], v̂; atol = 1e-2)

        # It is the transverse direction: perpendicular to the radius, in the
        # orbit plane, and on the same side as the motion.
        @test abs(dot(R[2, :], r)) < 1e-10
        @test dot(R[2, :], v) > 0
        @test rad2deg(acos(clamp(dot(R[2, :], v̂), -1, 1))) > 5.0
    end

    @testset "in-track points along the motion, not against it" begin
        # This is the check that was missing. The frame was built as
        # (R̂, Ĉ, R̂ × Ĉ), which is orthonormal with determinant +1 and passes
        # every structural test — while putting cross-track in the in-track
        # slot and pointing the third axis *against* the velocity.
        #
        # A circular orbit is what exposes it, because there in-track is
        # exactly the velocity direction and the sign has nowhere to hide.
        for a in (7000.0, 26_600.0, 42_164.0)
            state, _ = _circular(a)
            r = state[1:3]; v = state[4:6]
            R = axes_rotation(ICRF(), RIC(), _JD_OR, (; reference_state = state))[1:3, 1:3]

            @test R[2, :] ≈ v / norm(v)                     # +velocity, not −
            @test dot(R[2, :], v) > 0
            @test R[3, :] ≈ cross(r, v) / norm(cross(r, v))

            # Stated the way a user meets it: a point one kilometre ahead of
            # the reference reads +1 km in the second component.
            ahead = v / norm(v)
            @test R * ahead ≈ [0.0, 1.0, 0.0] atol = 1e-12
        end
    end

    @testset "a radial-primary frame needs no acceleration" begin
        # The practical difference from VNB, and the reason both exist.
        state, accel = _circular(7000.0)
        M_no_accel = axes_rotation(ICRF(), RIC(), _JD_OR, (; reference_state = state))
        M_accel    = axes_rotation(ICRF(), RIC(), _JD_OR,
                                   (; reference_state = state, reference_accel = accel))

        # Two-body acceleration is radial, so r × a = 0 and the orbit plane is
        # fixed: supplying it changes nothing. That is exactly why omitting it
        # is safe here and not safe for VNB.
        @test M_no_accel ≈ M_accel
    end

    @testset "for a circular orbit the frame turns at the orbital rate" begin
        for a in (7000.0, 26_600.0, 42_164.0)
            state, _ = _circular(a)
            M = axes_rotation(ICRF(), RIC(), _JD_OR, (; reference_state = state))
            w, W = _rate_vector(M)
            @test norm(W + W') < 1e-18
            @test norm(w) ≈ sqrt(_MU / a^3) rtol = 1e-12
        end
    end

    @testset "rate is h/r² for an eccentric orbit too" begin
        # The general statement, of which the circular case is one point: a
        # radial-primary frame turns with the true anomaly.
        state, _ = _generic()
        r = state[1:3]; v = state[4:6]
        M = axes_rotation(ICRF(), RIC(), _JD_OR, (; reference_state = state))
        w, _ = _rate_vector(M)
        @test norm(w) ≈ norm(cross(r, v)) / dot(r, r) rtol = 1e-12
    end

    @testset "out-of-plane acceleration turns the orbit plane" begin
        # And is the only thing that does. With acceleration supplied the rate
        # picks up a component along the radius — the plane rotating.
        state, _ = _generic()
        M = axes_rotation(ICRF(), RIC(), _JD_OR,
                          (; reference_state = state, reference_accel = [0.0, 0.0, 1e-6]))
        w, _ = _rate_vector(M)
        rhat = state[1:3] / norm(state[1:3])
        @test abs(dot(w, rhat)) > 1e-12
    end

    @testset "round trip" begin
        state, _ = _generic()
        p = (; reference_state = state)
        @test norm(axes_rotation(ICRF(), RIC(), _JD_OR, p) *
                   axes_rotation(RIC(), ICRF(), _JD_OR, p) - I) < 1e-14
    end
end

@testset "orbit-relative axes (LVLH)" begin

    @testset "the stated sign convention, pinned" begin
        # LVLH signs vary by source, so the docstring commits to one and this
        # test is what makes that commitment real: z nadir, y = -h, x = y × z.
        state, _ = _generic()
        r = state[1:3]; v = state[4:6]
        hhat = cross(r, v) / norm(cross(r, v))
        rhat = r / norm(r)
        R = axes_rotation(ICRF(), LVLH(), _JD_OR, (; reference_state = state))[1:3, 1:3]

        @test R[3, :] ≈ -rhat                      # z is nadir
        @test R[2, :] ≈ -hhat                      # y is the negative orbit normal
        @test R[1, :] ≈ cross(R[2, :], R[3, :])
        @test norm(R' * R - I) < 1e-15
        @test det(R) ≈ 1

        # x must be roughly along-track, not anti-track. A sign slip in the
        # triad flips it by 180 degrees and every other check here still passes.
        @test dot(R[1, :], v) > 0
    end

    @testset "LVLH and RIC describe the same motion" begin
        # Same reference orbit, so the two frames differ by a fixed rotation
        # and must share a rotation rate.
        state, _ = _generic()
        p = (; reference_state = state)
        w_ric,  _ = _rate_vector(axes_rotation(ICRF(), RIC(),  _JD_OR, p))
        w_lvlh, _ = _rate_vector(axes_rotation(ICRF(), LVLH(), _JD_OR, p))
        @test norm(w_ric) ≈ norm(w_lvlh) rtol = 1e-12
    end

    @testset "for a circular orbit the frame turns at the orbital rate" begin
        for a in (7000.0, 42_164.0)
            state, _ = _circular(a)
            M = axes_rotation(ICRF(), LVLH(), _JD_OR, (; reference_state = state))
            w, W = _rate_vector(M)
            @test norm(W + W') < 1e-18
            @test norm(w) ≈ sqrt(_MU / a^3) rtol = 1e-12
        end
    end

    @testset "round trip" begin
        state, _ = _generic()
        p = (; reference_state = state)
        @test norm(axes_rotation(ICRF(), LVLH(), _JD_OR, p) *
                   axes_rotation(LVLH(), ICRF(), _JD_OR, p) - I) < 1e-14
    end

    @testset "missing params fail by name" begin
        e = try; axes_rotation(ICRF(), LVLH(), _JD_OR); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("reference_state", e.msg)
        @test occursin("LVLH", e.msg)
    end
end

@testset "orbit-relative frames route" begin
    # The point of threading `params` through composition. Before this only the
    # direct ICRF-VNB pair worked, and ITRF to VNB raised.
    #
    # Needs EOP, unlike the rest of this file: the route runs through the Earth
    # chain to reach ICRF.

    state, accel = _generic()
    p = (; reference_state = state, reference_accel = accel)

    @testset "a multi-edge route reaches an orbit-relative frame" begin
        for target in (RIC(), LVLH(), VNB())
            M_routed = axes_rotation(ITRF(), target, _JD_OR, p)

            # Equal to the hand-composed chain: the parameters must survive
            # every edge and be consumed by the last one.
            M_hand = axes_rotation(ICRF(), target, _JD_OR, p) *
                     axes_rotation(ITRF(), ICRF(), _JD_OR)
            @test M_routed ≈ M_hand
            @test norm(M_routed[1:3,1:3]' * M_routed[1:3,1:3] - I) < 1e-13
        end
    end

    @testset "routing carries the rate, not just the orientation" begin
        # ITRF turns with the Earth and RIC with the orbit, so the composed
        # rate is the difference. Dropping Rdot anywhere in the chain would
        # leave the orientation right and this wrong.
        M = axes_rotation(ITRF(), RIC(), _JD_OR, p)
        w, W = _rate_vector(M)
        @test norm(W + W') < 1e-18
        @test 1e-5 < norm(w) < 1e-3
    end

    @testset "one orbit-relative frame to another" begin
        M = axes_rotation(RIC(), VNB(), _JD_OR, p)
        @test M ≈ axes_rotation(ICRF(), VNB(), _JD_OR, p) *
                  axes_rotation(RIC(), ICRF(), _JD_OR, p)
        @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-14
    end

    @testset "Robustness: the reverse direction fails by name too" begin
        # `VNB -> ICRF` is as unusable without a reference orbit as
        # `ICRF -> VNB`, and only the forward direction was checked.
        for F in (RIC(), LVLH(), VNB())
            e = try; axes_rotation(F, ICRF(), _JD_OR); nothing; catch e; e; end
            @test e isa ArgumentError
            @test occursin("reference_state", e.msg)
        end
    end

    @testset "out of theory, with parameters, warns and proceeds" begin
        # The parameterless router had this path tested; the one that carries
        # parameters did not, and it is a separate branch.
        original = frame_theory()
        try
            set_frame_theory!(IAU2006())
            M = axes_rotation(MODEq(), RIC(), _JD_OR, p)   # MODEq is FK5-only
            @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
        finally
            set_frame_theory!(original)
        end
    end

    @testset "a route into VNB carries the acceleration through" begin
        # The long way round must reach the same frame as the direct edge —
        # params have to survive every hop, not just the last one.
        accel = (-398600.4418 / norm(state[1:3])^3) .* state[1:3]

        for p in ((; reference_state = state),
                  (; reference_state = state, reference_accel = accel))
            routed = axes_rotation(ITRF(), VNB(), _JD_OR, p)
            direct = axes_rotation(ICRF(), VNB(), _JD_OR, p) *
                     axes_rotation(ITRF(), ICRF(), _JD_OR)
            @test routed ≈ direct
        end

        # And the two differ from each other, so the check above is not
        # comparing zero with zero.
        @test !(axes_rotation(ITRF(), VNB(), _JD_OR, (; reference_state = state)) ≈
                axes_rotation(ITRF(), VNB(), _JD_OR,
                              (; reference_state = state, reference_accel = accel)))
    end
end
