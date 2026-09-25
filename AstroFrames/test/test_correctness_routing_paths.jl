# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# Routing paths and a failure message that the rest of the suite does not reach.
#
# Truth: composition. Every route below is checked against the product of rotations that are
# themselves verified elsewhere in this suite against ERFA or SPICE, so a route that composes in
# the wrong order or through the wrong frame disagrees with the product.
#
#   - an extension frame whose hub is not ICRF, which is the only way to reach the leaf-to-leaf
#     route through two different hubs, since every shipped leaf frame uses ICRF
#   - a transform that crosses from the FK5 chain to the IAU 2006 chain, whose existing test
#     checks only that the result is orthonormal
#   - the same crossing through the method that takes orbit-relative parameters, a second copy of
#     the routing logic
#   - a subject that holds its state in the `OrbitState` container, as a Spacecraft does
#   - the message raised when a lunar frame cannot be evaluated

using AstroFrames
using AstroUniverse
using AstroEpochs
using AstroStates
using LinearAlgebra
using Test

const _RP_JD  = 2458849.5
const _RP_TOL = 1e-13

# A frame fixed to Earth's rotating axes by a constant 30° turn about the pole. Its hub is ITRF,
# which is interior to the Earth chains, so reaching it from a lunar frame routes MoonME → ICRF
# (MoonME's hub), ICRF → ITRF (the chain between the two hubs), ITRF → OffsetFixed.
struct _OffsetFixed <: AstroFrames.AbstractAxes end

AstroFrames.hub_axes(::Type{_OffsetFixed}) = ITRF()

function _rot_z6(θ)
    c, s = cos(θ), sin(θ)
    R = [c s 0.0; -s c 0.0; 0.0 0.0 1.0]
    return [R zeros(3, 3); zeros(3, 3) R]
end

AstroFrames.axes_rotation(::ITRF, ::_OffsetFixed, e::EpochScales) = _rot_z6(deg2rad(30.0))
AstroFrames.axes_rotation(::_OffsetFixed, ::ITRF, e::EpochScales) = _rot_z6(deg2rad(-30.0))

@testset "an extension frame on a non-ICRF hub routes through both hubs" begin
    M = axes_rotation(MoonME(), _OffsetFixed(), _RP_JD)
    truth = axes_rotation(ITRF(), _OffsetFixed(), _RP_JD) * axes_rotation(MoonME(), ITRF(), _RP_JD)
    @test M ≈ truth atol = _RP_TOL

    back = axes_rotation(_OffsetFixed(), MoonME(), _RP_JD)
    @test back ≈ axes_rotation(ITRF(), MoonME(), _RP_JD) *
                 axes_rotation(_OffsetFixed(), ITRF(), _RP_JD) atol = _RP_TOL

    # Going and coming back is the identity on the state.
    x = [7000.0, 1000.0, -2000.0, 1.0, 7.0, 0.5]
    @test back * (M * x) ≈ x rtol = 1e-12
end

@testset "a transform across the two Earth theories is the product of its halves" begin
    # Under IAU 2006, MODEq is reachable only through ICRF, which both chains share.
    original = frame_theory()
    try
        set_frame_theory!(IAU2006())
        empty!(AstroFrames._OUT_OF_THEORY_WARNED)
        M = @test_logs (:warn,) axes_rotation(MODEq(), CIRS(), _RP_JD)
        truth = axes_rotation(ICRF(), CIRS(), _RP_JD) * axes_rotation(MODEq(), ICRF(), _RP_JD)
        @test M ≈ truth atol = _RP_TOL

        # The method that takes orbit-relative parameters carries its own copy of the routing.
        empty!(AstroFrames._OUT_OF_THEORY_WARNED)
        Mp = @test_logs (:warn,) axes_rotation(MODEq(), CIRS(), _RP_JD, NamedTuple())
        @test Mp ≈ M atol = _RP_TOL
    finally
        set_frame_theory!(original)
    end
end

# A subject of the kind AstroModels' Spacecraft is: its state in the tagged container.
struct _ContainerSubject
    state::OrbitState
    frame::CoordinateSystem
    epoch::Time
end
AstroFrames.state_of(s::_ContainerSubject) = s.state
AstroFrames.frame_of(s::_ContainerSubject) = s.frame
AstroFrames.epoch_of(s::_ContainerSubject) = s.epoch

@testset "a subject holding an OrbitState converts as its Cartesian state does" begin
    epoch = Time(_RP_JD, 0.0, TDB(), JD())
    cs    = CoordinateSystem(earth, MJ2000Eq())
    cart  = CartesianState([7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9])
    kep   = KeplerianState(cart, earth.mu)

    from_kep  = _ContainerSubject(OrbitState(to_vector(kep), Keplerian()), cs, epoch)
    from_cart = _ContainerSubject(OrbitState(to_vector(cart), Cartesian()), cs, epoch)

    target = CoordinateSystem(earth, ITRF())
    x_kep  = to_vector(CartesianState(Coordinate(from_kep, target)))
    x_cart = to_vector(CartesianState(Coordinate(from_cart, target)))
    @test x_kep ≈ x_cart rtol = 1e-10
    @test norm(x_kep[1:3]) ≈ norm(to_vector(cart)[1:3]) rtol = 1e-10
end

@testset "a lunar frame that cannot be evaluated says what to load" begin
    cause = ErrorException("SPICE(FRAMEDATANOTFOUND)")

    # A frame no kernel defines: the message is about the frame kernel.
    e1 = AstroFrames._lunar_frame_error("MOON_NOT_A_FRAME", 0.0, cause)
    @test e1 isa ArgumentError
    @test occursin("MoonME axes need the lunar frame kernels", e1.msg)
    @test occursin("moon_de440_250416.tf", e1.msg)
    @test occursin("FRAMEDATANOTFOUND", e1.msg)

    # A frame that is defined but has no data at the epoch: the message gives the epoch as a date.
    e2 = AstroFrames._lunar_frame_error("MOON_PA", 0.0, cause)
    @test occursin("MoonPA axes are defined", e2.msg)
    @test occursin("2000-01-01T12:00:00.000 TDB", e2.msg)
    @test occursin("moon_pa_de440_200625.bpc", e2.msg)

    # And the lookup raises that error, not SPICE's.
    @test_throws ArgumentError AstroFrames._lunar_sxform("J2000", "MOON_NOT_A_FRAME", 0.0)
end

nothing
