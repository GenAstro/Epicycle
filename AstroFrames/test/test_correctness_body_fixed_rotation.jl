# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: MIT

# `body_fixed_rotation(model, naifid, epoch)`: the rotation into the fixed axes an orientation
# model defines, for any model, Earth's frame theories included. It must agree with the axes the
# registered models give, since both are the same rotation reached two ways.

using AstroFrames
using AstroUniverse
using AstroEpochs
using LinearAlgebra
using Test

@testset "body_fixed_rotation — any orientation model" begin
    t = Time("2024-01-01T00:00:00", UTC(), ISOT())
    original = frame_theory()
    try
        @testset "Earth: the frame theory's ITRF chain, whichever theory is active" begin
            for th in (IAU2006(), FK5())
                set_frame_theory!(th)
                @test body_fixed_rotation(th, 399, t) == axes_rotation(ICRF(), ITRF(), t)
                # Body-fixed axes for the Earth are ITRF, now that it has an orientation model.
                @test axes_rotation(ICRF(), CelestialBodyFixed{399}(), t) ==
                      axes_rotation(ICRF(), ITRF(), t)
            end
            # A theory other than the active one is evaluated as asked, not as set.
            set_frame_theory!(IAU2006())
            fk5 = body_fixed_rotation(FK5(), 399, t)
            @test fk5 != axes_rotation(ICRF(), ITRF(), t)
            set_frame_theory!(FK5())
            @test fk5 == axes_rotation(ICRF(), ITRF(), t)
        end

        @testset "the Moon: LunarPA and LunarME are the MoonPA and MoonME axes" begin
            @test body_fixed_rotation(LunarPA(), 301, t) == axes_rotation(ICRF(), MoonPA(), t)
            @test body_fixed_rotation(LunarME(), 301, t) == axes_rotation(ICRF(), MoonME(), t)
            @test axes_rotation(ICRF(), CelestialBodyFixed{301}(), t) ==
                  axes_rotation(ICRF(), MoonPA(), t)
        end

        @testset "Earth and the Moon now have body-fixed coordinate systems" begin
            # They raised before AstroUniverse 0.4, having no orientation model.
            cs = CoordinateSystem(earth, CelestialBodyFixed())
            @test cs.axes === CelestialBodyFixed{399}()
            @test axes_rotation(ICRF(), cs.axes, t) == axes_rotation(ICRF(), ITRF(), t)
            cs = CoordinateSystem(moon, CelestialBodyFixed())
            @test cs.axes === CelestialBodyFixed{301}()
            @test axes_rotation(ICRF(), cs.axes, t) == axes_rotation(ICRF(), MoonPA(), t)
        end

        @testset "a planet: the model's own rotation, at the epoch's TDB" begin
            @test body_fixed_rotation(IAU1991(), 499, t) == body_axes_rotation(IAU1991(), 499, t.tdb.jd)
            @test body_fixed_rotation(IAU2015(), 499, t) == axes_rotation(ICRF(), CelestialBodyFixed{499}(), t)
            # A TDB Julian date is accepted as the epoch too.
            @test body_fixed_rotation(IAU1991(), 499, t.tdb.jd) == body_fixed_rotation(IAU1991(), 499, t)
        end
    finally
        set_frame_theory!(original)
    end
end

using ForwardDiff

@testset "a Julian date carrying a derivative works as a Time does" begin
    # A dual Julian date used to recurse in `_scales` until the stack overflowed.
    t0 = Time(2460000.5, 0.0, TDB(), JD())
    viajd(x)   = axes_rotation(ICRF(), ITRF(), x)[1, 2]
    viatime(x) = axes_rotation(ICRF(), ITRF(), Time(x, zero(x), TDB(), JD()))[1, 2]
    d = ForwardDiff.derivative(viajd, 2460000.5)
    @test d == ForwardDiff.derivative(viatime, 2460000.5)
    # The truth is the rate block, per second, times the seconds in a day. A finite
    # difference on a date near 2.46e6 would resolve only 40 µs. The Earth's block
    # leaves out the precession, nutation and polar-motion rates that differentiation
    # includes; they differ by 2.5e-6 of the spin here.
    @test d ≈ axes_rotation(ICRF(), ITRF(), t0)[4, 2] * 86400 rtol = 1e-5
    mars(x) = body_fixed_rotation(IAU1991(), 499, x)[1, 1]
    @test ForwardDiff.derivative(mars, 2460000.5) ≈
          body_fixed_rotation(IAU1991(), 499, t0)[4, 1] * 86400 rtol = 1e-9
    # SPICE-backed axes cannot differentiate, and say so rather than blame the kernels.
    e = try; ForwardDiff.derivative(x -> axes_rotation(ICRF(), MoonPA(), x)[1, 1], 2460000.5); nothing
        catch e; e; end
    @test e isa ArgumentError && occursin("not differentiable", e.msg)
end

@testset "a frame theory gives the Earth's axes and no other body's" begin
    t = Time("2024-01-01T00:00:00", UTC(), ISOT())
    for th in (IAU2006(), FK5())
        e = try; body_fixed_rotation(th, 499, t); nothing; catch e; e; end
        @test e isa ArgumentError && occursin("ITRF", e.msg) && occursin("499", e.msg)
    end
end

@testset "an epoch of the wrong type is named as such" begin
    e = try; body_fixed_rotation(IAU2006(), 399, "2024-01-01"); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("String", e.msg) && occursin("Time", e.msg)
    # An integer Julian date is still fine.
    @test body_fixed_rotation(IAU2015(), 499, 2451545) == body_fixed_rotation(IAU2015(), 499, 2451545.0)
end
