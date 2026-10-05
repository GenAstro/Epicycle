# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Which axes AtmosphericDrag places the atmosphere in, the body's orientation model, which for the
# Earth is the frame theory at construction; which atmospheres it accepts, a model of its own
# body's; and the wind, which is the body's own spin.

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroUniverse
using Test
using LinearAlgebra: cross, norm

@testset "AtmosphericDrag — the body's axes" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    x = [6578.137, 100.0, 50.0, 0.0, 7.7, 0.1]
    original = frame_theory()
    try
        set_frame_theory!(FK5())
        @test AtmosphericDrag(earth).orientation === FK5()
        set_frame_theory!(IAU2006())
        drag = AtmosphericDrag(earth)
        @test drag.orientation === IAU2006()

        # A force takes no axes of its own.
        @test_throws MethodError AtmosphericDrag(earth; model = Exponential(), orientation = FK5())

        # The density sees the rotation the force chose.
        R = body_fixed_rotation(IAU2006(), 399, t)[1:3, 1:3]
        @test AstroProp.density(Exponential(), t.utc.jd, x, R) ==
              AstroProp._exponential_density(AstroProp._geodetic(t.utc.jd, x, R)[3])
    finally
        set_frame_theory!(original)
    end
end

# An atmosphere of a user's own, of Mars: constant density, so the drag is exact.
struct _MarsConstant <: AbstractDensityModel end
AstroProp.atmosphere_body(::_MarsConstant) = 499
AstroProp.density(::_MarsConstant, jd, x̄, axes) = 1.0e-12

@testset "AtmosphericDrag — a model of the body's own atmosphere" begin
    @test atmosphere_body(Exponential()) == 399
    # Another body's atmosphere is refused, naming both bodies.
    for body in (mars, moon)
        e = try; AtmosphericDrag(body; model = Exponential()); nothing; catch e; e; end
        @test e isa ArgumentError && occursin("Earth", e.msg) && occursin(body.name, e.msg)
    end
    e = try; AtmosphericDrag(earth; model = _MarsConstant()); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("NAIF 499", e.msg) && occursin("Earth", e.msg)

    # Mars's own atmosphere works, turning with Mars's axes.
    drag = AtmosphericDrag(mars; model = _MarsConstant())
    @test drag.orientation === orientation_model(mars)
    t  = Time("2024-03-15T00:00:00", UTC(), ISOT())
    x  = [3700.0, 0.0, 0.0, 0.0, 3.4, 0.0]
    sc = Spacecraft(state = CartesianState(x), time = t, coord_sys = CoordinateSystem(mars, ICRF()),
                    mass = 100.0, drag = SphericalDrag(c_d = 2.0, drag_area = 1.0))
    a  = AstroProp.accel_eval!(drag, t, x, zeros(6), sc, nothing)[4:6]
    R, Ṙ = AstroProp._rotation_blocks(body_fixed_rotation(orientation_model(mars), 499, t))
    v_rel = R' * (R * x[4:6] + Ṙ * x[1:3])
    @test a ≈ -0.5 * 1.0e-12 * (2.0 * 1.0 / 100.0) * norm(v_rel) * v_rel * 1.0e3 rtol = 1e-12
end

@testset "AtmosphericDrag — the density seam takes the 3×3 rotation" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    M = body_fixed_rotation(IAU2006(), 399, t)
    e = try; AstroProp.density(Exponential(), t.utc.jd, [6578.137, 0.0, 0.0], M); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("3×3", e.msg) && occursin("6×6", e.msg)
end

@testset "AtmosphericDrag — the wind is the body's own spin" begin
    # v_rel = Rᵀ(R v + Ṙ r). In the frame theory's ITRF the spin is about the Earth's own pole, 0.12°
    # of precession from ICRF z by 2020, at the sidereal rate corrected for the length of day: so
    # the wind is v − ω × r with ω along that pole, to the length-of-day correction.
    t  = Time("2020-10-20T12:00:00", UTC(), ISOT())
    r  = [6578.137, 100.0, 50.0]
    v  = [0.0, 7.7, 0.1]
    R, Ṙ = AstroProp._rotation_blocks(body_fixed_rotation(IAU2006(), 399, t))
    pole = R' * [0.0, 0.0, 1.0]
    ω = AstroProp.EARTH_ANGULAR_SPEED * pole
    wind = R' * (R * v + Ṙ * r)
    @test norm(wind - (v - cross(ω, r))) < 1e-6 * norm(cross(ω, r))
    # About the ICRF z axis instead, the wind would be wrong by ω (p̂ − ẑ) × r, the pole's tilt.
    ωz = [0.0, 0.0, AstroProp.EARTH_ANGULAR_SPEED]
    tilt = norm(cross(ω - ωz, r))
    @test tilt > 0
    @test norm(wind - (v - cross(ωz, r))) ≈ tilt rtol = 0.05
end
