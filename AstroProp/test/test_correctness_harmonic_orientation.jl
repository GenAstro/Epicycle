# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Which body-fixed axes HarmonicGravity evaluates a field in: the field's own (`field_orientation`),
# which by default is the body's orientation model, and for the Earth the frame theory. The choice
# is made at construction, and no force takes axes of its own. Also the seams a field gives tides:
# its reference radius and tide system, and a tide model resolved against them.

using AstroProp
using AstroProp: field_orientation
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroUniverse
using LinearAlgebra: norm
using StaticArrays: SVector
using ForwardDiff
using Test

function _zonal_accel(force, t)
    x  = [6878.137, 100.0, 2000.0, 0.0, 7.5, 0.0]
    xd = zeros(6)
    sc = Spacecraft(state = CartesianState(x), time = t)
    return AstroProp.accel_eval!(force, t, x, xd, sc, nothing)[4:6]
end

@testset "HarmonicGravity — which axes a field is evaluated in" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    original = frame_theory()
    try
        @testset "the Earth's default is the frame theory at construction" begin
            set_frame_theory!(IAU2006())
            @test field_orientation(Zonal(), earth) === IAU2006()
            g06 = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
            @test g06.orientation === IAU2006()

            set_frame_theory!(FK5())
            gfk = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
            @test gfk.orientation === FK5()

            # A force built earlier keeps its axes, and the two differ by the theories' pole.
            @test g06.orientation === IAU2006()
            a06, afk = _zonal_accel(g06, t), _zonal_accel(gfk, t)
            @test a06 != afk
            @test norm(a06 - afk) / norm(a06) < 1e-6
        end

        @testset "a force takes no axes of its own" begin
            @test_throws MethodError HarmonicGravity(earth; degree = 5, order = 0, model = Zonal(),
                                                     orientation = IAU2006())
        end

        @testset "Zonal is the Earth's field only" begin
            e = try
                HarmonicGravity(mars; degree = 2, order = 0, model = Zonal())
                nothing
            catch e; e; end
            @test e isa ArgumentError
            @test occursin("Earth", e.msg) && occursin("Mars", e.msg)
        end
    finally
        set_frame_theory!(original)
    end
end

# A field of a user's own, declaring its axes: the point-mass term alone, which reads the same in
# any axes, so the rotation there and back is all that can change it.
struct _PointField <: AstroProp.AbstractGeopotential end
AstroProp.max_degree(::_PointField) = 0
AstroProp.max_order(::_PointField)  = 0
AstroProp.geopotential_data(::_PointField, body, degree, order) = body.mu * 1.0e9      # m³/s²
AstroProp.geopotential_accel(::_PointField, μ, r, tsec, degree, order) = -μ .* r ./ norm(r)^3
AstroProp.field_orientation(::_PointField, body) = IAU1991()

# The same field, declaring axes that are not Mars's.
struct _EarthAxesField <: AstroProp.AbstractGeopotential end
AstroProp.max_degree(::_EarthAxesField) = 0
AstroProp.max_order(::_EarthAxesField)  = 0
AstroProp.geopotential_data(::_EarthAxesField, body, degree, order) = body.mu * 1.0e9
AstroProp.geopotential_accel(::_EarthAxesField, μ, r, tsec, degree, order) = -μ .* r ./ norm(r)^3
AstroProp.field_orientation(::_EarthAxesField, body) = IAU2006()

@testset "HarmonicGravity — a user's field declares its axes" begin
    t = Time("2024-03-15T00:00:00", UTC(), ISOT())
    g = HarmonicGravity(mars; degree = 0, order = 0, model = _PointField())
    @test g.orientation === IAU1991()
    x  = [3794.2, 300.0, 800.0, 0.0, 0.0, 0.0]
    a  = AstroProp.accel_eval!(g, t, x, zeros(6), Spacecraft(state = CartesianState(x), time = t), nothing)[4:6]
    @test a ≈ -mars.mu .* x[1:3] ./ norm(x[1:3])^3 rtol = 1e-14
end

@testset "HarmonicGravity — a field declaring another body's axes is refused when built" begin
    e = try
        HarmonicGravity(mars; degree = 0, order = 0, model = _EarthAxesField())
        nothing
    catch e; e; end
    @test e isa ArgumentError && occursin("NAIF 499", e.msg)
end

# A tide model of a user's own: a constant acceleration, so what HarmonicGravity adds is exact.
struct _ConstantTide <: AstroProp.AbstractTideModel
    a::Float64
end
struct _ResolvedConstantTide
    a::Float64
    radius::Float64
    system::Symbol
end
AstroProp.resolve_tides(m::_ConstantTide, body, model, data, axes) =
    _ResolvedConstantTide(m.a, AstroProp.field_radius(model, data), AstroProp.tide_system(model, data))
AstroProp.tide_accel(m::_ResolvedConstantTide, r, R, t, params) = SVector{3}(m.a, 0.0, 0.0)

@testset "HarmonicGravity — tides are resolved against the field and added to it" begin
    @test AstroProp.tide_system(_PointField(), nothing) === :unknown
    data = AstroProp.geopotential_data(Zonal(), earth, 5, 0)
    @test AstroProp.field_radius(Zonal(), data) == 6378.1363
    @test AstroProp.tide_system(Zonal(), data) === :tide_free

    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    plain = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
    tidal = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal(), tides = _ConstantTide(1e-9))
    @test plain.tides === nothing
    @test tidal.tides.radius == 6378.1363 && tidal.tides.system === :tide_free
    @test _zonal_accel(tidal, t) - _zonal_accel(plain, t) ≈ [1e-9, 0.0, 0.0] atol = 1e-16   # round-off of the 8e-3 field term
end

@testset "HarmonicGravity — derivatives through time on the new path" begin
    # Through the epoch, as an estimator of the epoch or a free final time would take them.
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    g = HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
    x = [6878.137, 100.0, 2000.0, 0.0, 7.5, 0.0]
    sc = Spacecraft(state = CartesianState(x), time = t)
    f(dt) = AstroProp.accel_eval!(g, t + dt, x, zeros(typeof(dt), 6), sc, nothing)[4]
    # A zonal field is symmetric about the ITRF z axis, which polar motion carries once a day
    # around the celestial pole, so the acceleration has a small daily term and a central
    # difference errs by (ωh)²/6. Richardson's combination of two cancels that term.
    D(h) = (f(h) - f(-h)) / 2h                     # h in days; a Time keeps full precision
    # At h = 2e-3 the differences are near rounding: the extrapolation scattered by 1e-4 across
    # steps. At 5e-3 it agrees with ForwardDiff to 2e-7, and truncation is still cancelled.
    h = 5e-3
    @test ForwardDiff.derivative(f, 0.0) ≈ (4D(h / 2) - D(h)) / 3 rtol = 1e-6
end
