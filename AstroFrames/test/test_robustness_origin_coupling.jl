# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: MIT

# =============================================================================
# Refusing coordinate systems that are not coordinate systems.
#
# ITRF axes about the Sun is not a frame, it is a mistake, and the only place
# it can be caught cheaply is where the coordinate system is built. `I-11` in
# `test_correctness_invariants.jl` checks that a few such pairs are refused;
# this file checks that *every* restricted frame is refused, and that the
# message tells the user what to pass instead.
#
# The messages matter as much as the refusal: a
# rejection that does not say what would have worked sends the user to the
# source.
# =============================================================================

using AstroFrames
using AstroUniverse
using EpicycleBase
using Test

"""The ArgumentError raised by `f`, or `nothing` if it did not raise."""
_raised(f) = try; f(); nothing; catch e; e; end

# Every Earth-restricted frame, not a sample of them. Adding a frame to the
# family without adding it here leaves a hole this loop closes.
const _EARTH_ONLY = (GCRF(), CIRS(), TIRS(), ITRF(), MODEq(), TODEq(),
                     MODEc(), TODEc(), PEF(), TEME())

"""A point carrying no NAIF ID, which nothing shipped provides."""
struct _AnonymousPoint <: EpicycleBase.AbstractPoint end

@testset "origin coupling: Earth frames" begin

    @testset "accepted at Earth" begin
        for axes in _EARTH_ONLY
            @test CoordinateSystem(earth, axes) isa CoordinateSystem
        end
    end

    @testset "refused elsewhere, and the message names the fix" begin
        for axes in _EARTH_ONLY, origin in (sun, mars, moon)
            e = _raised(() -> CoordinateSystem(origin, axes))
            @test e isa ArgumentError
            # Names the frame, the origin it got, and what it wanted.
            @test occursin(string(nameof(typeof(axes))), e.msg)
            @test occursin("earth", e.msg)
        end
    end
end

@testset "origin coupling: Moon frames" begin

    @testset "accepted at the Moon" begin
        for axes in (MoonPA(), MoonME())
            @test CoordinateSystem(moon, axes) isa CoordinateSystem
        end
    end

    @testset "refused elsewhere" begin
        for axes in (MoonPA(), MoonME()), origin in (earth, sun, mars)
            e = _raised(() -> CoordinateSystem(origin, axes))
            @test e isa ArgumentError
            @test occursin("moon", e.msg)
        end
    end
end

@testset "origin coupling: body-fixed axes" begin

    @testset "the body must be the origin" begin
        # The trap this catches is Mars-fixed axes about Venus: both are valid
        # frames, the pair is not, and nothing downstream would notice.
        for (body, other) in ((mars, venus), (jupiter, saturn), (sun, pluto))
            @test CoordinateSystem(body, CelestialBodyFixed(body)) isa CoordinateSystem
            e = _raised(() -> CoordinateSystem(other, CelestialBodyFixed(body)))
            @test e isa ArgumentError
            @test occursin("NAIF", e.msg)
        end
    end

    @testset "a body with no orientation model is refused where it is asked for" begin
        # Not at transform time, several calls down.
        rock = CelestialBody("Unoriented Rock", 1e-9, 1.0, 0.0, 2888888)
        e = _raised(() -> CoordinateSystem(rock, CelestialBodyFixed()))
        @test e isa ArgumentError
        @test occursin("set_orientation!", e.msg)
    end

    @testset "an origin with no NAIF ID at all" begin
        # `CelestialBodyFixed()` has nothing to resolve against.
        e = _raised(() -> CoordinateSystem(_AnonymousPoint(), CelestialBodyFixed()))
        @test e isa ArgumentError
        @test occursin("NAIF", e.msg)
    end
end

@testset "the unresolved sentinel refuses to transform" begin
    # `CelestialBodyFixed()` names no body, so it cannot rotate anything. It
    # must say that rather than fail as a `MethodError` inside the routing.
    for (a, b) in ((ICRF(), CelestialBodyFixed()),
                   (CelestialBodyFixed(), ICRF()),
                   (ITRF(), CelestialBodyFixed()),
                   (CelestialBodyFixed(), ITRF()),
                   (CelestialBodyFixed(), CelestialBodyFixed()))
        e = _raised(() -> axes_rotation(a, b, 2458849.5))
        @test e isa ArgumentError
        @test occursin("Unresolved", e.msg)
        @test occursin("CoordinateSystem", e.msg)   # names the way out
    end
end

@testset "frames with no origin restriction accept any origin" begin
    # The default. An orbit-relative or inertial frame is meaningful anywhere,
    # and must not be caught by the coupling check.
    for axes in (ICRF(), MJ2000Eq(), MJ2000Ec(), Inertial(), RIC(), LVLH(), VNB())
        for origin in (earth, sun, mars, moon)
            @test CoordinateSystem(origin, axes) isa CoordinateSystem
        end
    end
end
