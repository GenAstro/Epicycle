# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Which Earth-fixed axes AtmosphericDrag places the atmosphere in: the `orientation` keyword, else
# the body's orientation model, which for the Earth is the frame theory at construction. The same
# rule as HarmonicGravity.

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroUniverse
using Test
using LinearAlgebra: cross, norm

@isdefined(GmatEarthAxes) || include("gmat_earth_axes.jl")

@testset "AtmosphericDrag — which Earth axes the atmosphere is in" begin
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    x = [6578.137, 100.0, 50.0, 0.0, 7.7, 0.1]
    original = frame_theory()
    try
        set_frame_theory!(FK5())
        @test AtmosphericDrag(earth).orientation === FK5()
        set_frame_theory!(IAU2006())
        drag = AtmosphericDrag(earth)
        @test drag.orientation === IAU2006()
        @test AtmosphericDrag(earth; orientation = GmatEarthAxes()).orientation === GmatEarthAxes()

        # The density sees the rotation the force chose.
        R = body_fixed_rotation(IAU2006(), 399, t)[1:3, 1:3]
        @test AstroProp.density(Exponential(), t.utc.jd, x, R) ==
              AstroProp._exponential_density(AstroProp._geodetic(t.utc.jd, x, R)[3])

        # A direct call with an EOP table takes SatelliteToolbox's FK5 route from mean J2000, which
        # is what GmatEarthAxes returns, so the two agree.
        G = body_fixed_rotation(GmatEarthAxes(), 399, t)[1:3, 1:3]
        eop = AstroProp.fetch_iers_eop()
        @test AstroProp.density(Exponential(), t.utc.jd, x, eop) ≈
              AstroProp.density(Exponential(), t.utc.jd, x, G) rtol = 1e-12
    finally
        set_frame_theory!(original)
    end
end

@testset "AtmosphericDrag — the Earth only, in the Earth's axes" begin
    # The density models and the geodetic altitude under them are the Earth's.
    for body in (mars, moon)
        e = try; AtmosphericDrag(body); nothing; catch e; e; end
        @test e isa ArgumentError && occursin("Earth", e.msg) && occursin(body.name, e.msg)
    end
    # Axes that are not the Earth's are refused when the force is built.
    e = try; AtmosphericDrag(earth; orientation = IAU2015()); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("NAIF 399", e.msg)
    # The density seam takes the 3×3 rotation, and says so when handed the 6×6.
    t = Time("2020-10-20T12:00:00", UTC(), ISOT())
    M = body_fixed_rotation(IAU2006(), 399, t)
    e = try; AstroProp.density(Exponential(), t.utc.jd, [6578.137, 0.0, 0.0], M); nothing; catch e; e; end
    @test e isa ArgumentError && occursin("3×3", e.msg) && occursin("6×6", e.msg)
end

@testset "AtmosphericDrag — the wind is the axes' own spin" begin
    # v_rel = Rᵀ(R v + Ṙ r). With GMAT's axes, whose rate block is a constant spin about the
    # inertial z axis, that is exactly v − ω × r, the form the force used before.
    t  = Time("2020-10-20T12:00:00", UTC(), ISOT())
    r  = [6578.137, 100.0, 50.0]
    v  = [0.0, 7.7, 0.1]
    R, Ṙ = AstroProp._rotation_blocks(body_fixed_rotation(GmatEarthAxes(), 399, t))
    ω = [0.0, 0.0, AstroProp.EARTH_ANGULAR_SPEED]
    @test R' * (R * v + Ṙ * r) ≈ v - cross(ω, r) rtol = 1e-14
    # In the frame theory's ITRF the pole is the Earth's own, 0.28° of precession from ICRF z
    # by 2020, and the rate includes the length of day: close to the constant spin, not equal.
    R, Ṙ = AstroProp._rotation_blocks(body_fixed_rotation(IAU2006(), 399, t))
    w = R' * (R * v + Ṙ * r) - (v - cross(ω, r))
    @test 0 < norm(w) < 0.01 * norm(cross(ω, r))
end
