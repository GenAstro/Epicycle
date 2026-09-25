# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# A body this package has never heard of.
#
# Phobos stands in for the case the design exists to serve: a body whose
# orientation is not in our tables, whose model the user supplies, and whose
# orientation parameters are being estimated. Nothing below is written by
# AstroFrames or AstroUniverse — it is the file a user would write, run as a
# test so it cannot quietly stop working.
#
# Written against the public interface only. If any of this needs an
# underscore name, the interface is not finished.
# =============================================================================

using AstroFrames
using AstroUniverse
using ForwardDiff
using LinearAlgebra
using Test

const _JD_UB = 2458849.5

# --- What a user writes -----------------------------------------------------

phobos = CelestialBody("Phobos", 7.087e-4, 11.2667, 0.0, 401)

"""A pole that drifts and a body that spins about it — the usual shape."""
struct SimpleSpin{T<:Real} <: AbstractOrientationModel
    pole_ra::T          # deg
    pole_dec::T         # deg
    pm0::T              # deg, prime meridian at J2000
    spin_rate::T        # deg/day
end
SimpleSpin(a, d, w, s) = SimpleSpin(promote(a, d, w, s)...)

function AstroUniverse.body_axes_rotation(m::SimpleSpin, naifid, jd_tdb)
    d  = jd_tdb - 2451545.0
    Ẇ  = deg2rad(m.spin_rate)                      # rad/day
    return pole_axes_rotation(deg2rad(m.pole_ra), deg2rad(m.pole_dec),
                              deg2rad(m.pm0) + Ẇ * d,
                              zero(Ẇ), zero(Ẇ), Ẇ / 86_400)
end

struct PhobosFixed <: AstroFrames.AbstractAxes end

AstroFrames.hub_axes(::Type{PhobosFixed}) = ICRF()
AstroFrames.valid_origin(::PhobosFixed, o) = o.naifid == 401

function AstroFrames.axes_rotation(::ICRF, ::PhobosFixed, e::EpochScales)
    return body_axes_rotation(orientation_model(phobos), 401, epoch_tdb(e))
end
AstroFrames.axes_rotation(a::PhobosFixed, ::ICRF, e::EpochScales) =
    inv(axes_rotation(ICRF(), a, e))

# IAU 2015 values for Phobos, stripped of the libration terms. Good enough to
# be a body; not a truth row, and not claimed as one.
const _PHOBOS_MODEL = SimpleSpin(317.68, 52.90, 35.06, 1128.8445850)

set_orientation!(phobos, _PHOBOS_MODEL)

# --- What must then be true -------------------------------------------------

@testset "a user-defined body" begin

    @testset "registering a model is all it takes" begin
        @test orientation_model(phobos) === _PHOBOS_MODEL
        @test has_orientation_model(phobos)

        # Before registration there was nothing, and the error said so rather
        # than failing somewhere further down.
        unknown = CelestialBody("Chariklo", 1.0e-10, 100.0, 0.0, 2010199)
        @test !has_orientation_model(unknown)
        e = try; orientation_model(unknown); nothing; catch e; e; end
        @test e isa ArgumentError
        @test occursin("set_orientation!", e.msg)
    end

    @testset "the spin rate is what was asked for" begin
        M = axes_rotation(ICRF(), PhobosFixed(), _JD_UB)
        Ω = M[4:6, 1:3] * M[1:3, 1:3]'
        @test norm(Ω + Ω') < 1e-18                          # a rate is skew
        ω = sqrt(Ω[3,2]^2 + Ω[1,3]^2 + Ω[2,1]^2)
        @test 2π / ω / 3600 ≈ 360 / 1128.8445850 * 24 rtol = 1e-10   # 7.653 h
    end

    @testset "the pole points where it was told to" begin
        # Row 3 of R is the body's north pole in ICRF, by construction.
        R  = axes_rotation(ICRF(), PhobosFixed(), _JD_UB)[1:3, 1:3]
        α, δ = deg2rad(317.68), deg2rad(52.90)
        @test R[3, :] ≈ [cos(δ)cos(α), cos(δ)sin(α), sin(δ)] atol = 1e-12
        @test norm(R' * R - I) < 1e-14
        @test det(R) ≈ 1
    end

    @testset "orthonormal and reversible" begin
        M = axes_rotation(ICRF(), PhobosFixed(), _JD_UB)
        @test norm(M * axes_rotation(PhobosFixed(), ICRF(), _JD_UB) - I) < 1e-12
    end
end

@testset "a user-defined body routes" begin
    # Declaring `hub_axes` is the whole of what makes this work. The user wrote
    # one edge, to ICRF, and never touched the graph.

    @testset "reached from frames the user never mentioned" begin
        original = frame_theory()
        try
            set_frame_theory!(FK5())
            for source in (MJ2000Eq(), MODEq(), TODEq(), PEF(), ITRF(), TEME(), MJ2000Ec())
                M = axes_rotation(source, PhobosFixed(), _JD_UB)
                @test size(M) == (6, 6)
                @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
            end

            set_frame_theory!(IAU2006())
            for source in (GCRF(), CIRS(), TIRS(), ITRF())
                M = axes_rotation(source, PhobosFixed(), _JD_UB)
                @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
            end
        finally
            set_frame_theory!(original)
        end
    end

    @testset "out of theory it warns and proceeds, as everywhere else" begin
        # A frame hanging off a hub inherits the hub's reachability, so this
        # has to behave like any other cross-theory pair (FR-FRAME-4) rather
        # than being the one place the rule quietly does not apply. It did
        # throw here until the hub route got its own fallback.
        original = frame_theory()
        try
            set_frame_theory!(IAU2006())
            M = axes_rotation(MODEq(), PhobosFixed(), _JD_UB)   # MODEq is FK5-only
            @test norm(M[1:3,1:3]' * M[1:3,1:3] - I) < 1e-12
        finally
            set_frame_theory!(original)
        end
    end

    @testset "the routed answer equals the hand-composed one" begin
        original = frame_theory()
        try
            for theory in (FK5(), IAU2006())
                set_frame_theory!(theory)
                M_routed = axes_rotation(ITRF(), PhobosFixed(), _JD_UB)
                M_hand   = axes_rotation(ICRF(), PhobosFixed(), _JD_UB) *
                           axes_rotation(ITRF(), ICRF(), _JD_UB)
                @test M_routed ≈ M_hand
            end
        finally
            set_frame_theory!(original)
        end
    end

    @testset "to another frame that hangs off the same hub" begin
        # Leaf to leaf, through the hub: Phobos-fixed to a spacecraft's
        # orbit frame, neither of which is in the graph.
        state = [-4550.0, 2220.0, 4980.0, -3.10, -6.60, 0.12]
        M = axes_rotation(PhobosFixed(), RIC(), _JD_UB, (; reference_state = state))
        @test M ≈ axes_rotation(ICRF(), RIC(), _JD_UB, (; reference_state = state)) *
                  axes_rotation(PhobosFixed(), ICRF(), _JD_UB)
    end

    @testset "and back the other way" begin
        M = axes_rotation(PhobosFixed(), ITRF(), _JD_UB)
        @test norm(M * axes_rotation(ITRF(), PhobosFixed(), _JD_UB) - I) < 1e-12
    end

    @testset "the coordinate system builds, and refuses elsewhere" begin
        cs = CoordinateSystem(phobos, PhobosFixed())
        @test cs.axes isa PhobosFixed
        @test_throws ArgumentError CoordinateSystem(earth, PhobosFixed())
    end
end

@testset "orientation parameters are readable and writable" begin
    # AstroFrames does not know which of these are being solved for. It only
    # has to let whatever does read them and write them back.

    @testset "read by name, off the body" begin
        p = orientation_parameters(phobos)
        @test keys(p) == (:pole_ra, :pole_dec, :pm0, :spin_rate)
        @test p.pole_ra ≈ 317.68
    end

    @testset "writing changes the transform" begin
        before = axes_rotation(ICRF(), PhobosFixed(), _JD_UB)
        set_orientation_parameters!(phobos, (; pole_ra = 318.68))
        after = axes_rotation(ICRF(), PhobosFixed(), _JD_UB)

        @test orientation_parameters(phobos).pole_ra ≈ 318.68
        @test orientation_parameters(phobos).spin_rate ≈ 1128.8445850  # untouched
        @test !(before ≈ after)

        set_orientation!(phobos, _PHOBOS_MODEL)
        @test axes_rotation(ICRF(), PhobosFixed(), _JD_UB) ≈ before
    end

    @testset "a name that is not a parameter says so" begin
        e = try
            set_orientation_parameters!(phobos, (; oblateness = 1.0))
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("oblateness", e.msg)
        @test occursin("pole_ra", e.msg)      # names what it does have
    end

    @testset "a SPICE model has nothing to solve for, and says why" begin
        @test isempty(orientation_parameters(SpiceOrientation("IAU_MARS")))
        e = try
            set_orientation_parameters(SpiceOrientation("IAU_MARS"), (; pole_ra = 1.0))
            nothing
        catch e; e; end
        @test e isa ArgumentError
        @test occursin("kernel", e.msg)
    end
end

@testset "a user-defined body differentiates" begin
    # The reason a custom model has to exist at all: an estimated orientation
    # is a solve-for, and SPICE cannot carry a derivative.

    @testset "with respect to the pole being estimated" begin
        function R11(q)
            m = SimpleSpin(q[1], q[2], 35.06, 1128.8445850)
            return body_axes_rotation(m, 401, _JD_UB)[1, 1]
        end
        g = ForwardDiff.gradient(R11, [317.68, 52.90])
        @test all(isfinite, g)
        @test norm(g) > 1e-6          # the pole actually moves the frame

        # Against a difference quotient, so this is not just "it ran".
        h = 1e-6
        fd = (R11([317.68 + h, 52.90]) - R11([317.68 - h, 52.90])) / 2h
        @test g[1] ≈ fd rtol = 1e-5
    end

    @testset "with respect to the spin rate" begin
        f = s -> body_axes_rotation(SimpleSpin(317.68, 52.90, 35.06, s), 401, _JD_UB)[1, 1]
        d = ForwardDiff.derivative(f, 1128.8445850)
        @test isfinite(d) && abs(d) > 1e-6
    end

    @testset "and through the whole routed chain" begin
        # Not just the model in isolation: the derivative survives five Earth
        # edges, the EOP interpolation in them, and the composition.
        function routed(α)
            set_orientation!(phobos, SimpleSpin(α, 52.90, 35.06, 1128.8445850))
            return axes_rotation(ITRF(), PhobosFixed(), _JD_UB)[1, 1]
        end
        d = ForwardDiff.derivative(routed, 317.68)
        set_orientation!(phobos, _PHOBOS_MODEL)

        @test isfinite(d)
        h = 1e-6
        fd = (routed(317.68 + h) - routed(317.68 - h)) / 2h
        set_orientation!(phobos, _PHOBOS_MODEL)
        @test d ≈ fd rtol = 1e-5
    end

    @testset "through the estimator's actual write path" begin
        # Not by rebuilding the model by hand, but by writing a solved-for
        # value onto the body the way whatever manages an estimation would.
        # The model's fields are all one type, so writing a single derivative
        # value has to widen it — that failed until `set_orientation_parameters`
        # promoted the numeric fields.
        function by_write(α)
            set_orientation_parameters!(phobos, (; pole_ra = α))
            return axes_rotation(ICRF(), PhobosFixed(), _JD_UB)[1, 1]
        end
        d = ForwardDiff.derivative(by_write, 317.68)
        set_orientation!(phobos, _PHOBOS_MODEL)

        @test isfinite(d) && abs(d) > 1e-6
        h = 1e-6
        fd = (by_write(317.68 + h) - by_write(317.68 - h)) / 2h
        set_orientation!(phobos, _PHOBOS_MODEL)
        @test d ≈ fd rtol = 1e-5
    end

    @testset "the SPICE model cannot, which is the point" begin
        # Recorded so nobody builds an estimator on a SPICE-backed frame and
        # finds out late. This is a property of SPICE, not a defect here.
        f = x -> body_axes_rotation(SpiceOrientation("IAU_MARS"), 499, x)[1, 1]
        @test_throws MethodError ForwardDiff.derivative(f, _JD_UB)
    end
end
