# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Ground stations, spacecraft geometry, and a spacecraft's state in another frame.
#
# Truth, by section:
#   geodetic_to_body_fixed   analytic. The equator, the poles and the reference meridian have
#                            closed-form positions on an ellipsoid, and altitude adds along the
#                            ellipsoid normal.
#   get_state                analytic. A body-fixed station has zero body-fixed velocity, and an
#                            inertial rotation keeps |r| and gives |v| = ω·ρ, with ρ the distance
#                            from the spin axis, up to polar motion and length-of-day variation.
#   is_visible               analytic. Positions built at a stated elevation above the station's
#                            horizon, either side of the cutoff.
#   geometry, frames         by construction: stored values, printed forms, and a conversion that
#                            round-trips.

using Test
using LinearAlgebra

using AstroEpochs
using AstroStates
using AstroFrames
using AstroUniverse
using AstroModels
using AstroModels: geodetic_to_body_fixed

const _GS_EPOCH = Time("2020-06-01T00:00:00", UTC(), ISOT())

@testset "GroundStation — geodetic to body-fixed on the ellipsoid" begin
    a = earth.equatorial_radius
    f = earth.flattening
    b = a * (1 - f)

    @test geodetic_to_body_fixed(earth, 0.0, 0.0, 0.0, Ellipsoid())   ≈ [a, 0, 0]  atol = 1e-9
    @test geodetic_to_body_fixed(earth, 0.0, 90.0, 0.0, Ellipsoid())  ≈ [0, a, 0]  atol = 1e-9
    @test geodetic_to_body_fixed(earth, 0.0, 180.0, 0.0, Ellipsoid()) ≈ [-a, 0, 0] atol = 1e-9
    @test geodetic_to_body_fixed(earth, 90.0, 0.0, 0.0, Ellipsoid())  ≈ [0, 0, b]  atol = 1e-9
    @test geodetic_to_body_fixed(earth, -90.0, 0.0, 0.0, Ellipsoid()) ≈ [0, 0, -b] atol = 1e-9

    # Altitude adds along the normal, which at the pole and the equator is radial.
    @test geodetic_to_body_fixed(earth, 90.0, 0.0, 2.5, Ellipsoid()) ≈ [0, 0, b + 2.5] atol = 1e-9
    @test geodetic_to_body_fixed(earth, 0.0, 0.0, 2.5, Ellipsoid())  ≈ [a + 2.5, 0, 0] atol = 1e-9

    # At 45° the point sits on the ellipsoid: x²/a² + z²/b² = 1.
    r = geodetic_to_body_fixed(earth, 45.0, 30.0, 0.0, Ellipsoid())
    ρ = hypot(r[1], r[2])
    @test ρ^2 / a^2 + r[3]^2 / b^2 ≈ 1.0 atol = 1e-12
    @test atan(r[2], r[1]) ≈ deg2rad(30.0) atol = 1e-12
    # The geodetic latitude is the angle of the surface normal, (x/a², z/b²), not of the radius.
    @test atan(r[3] / b^2, ρ / a^2) ≈ deg2rad(45.0) atol = 1e-12
    @test atan(r[3], ρ) < deg2rad(45.0)                 # geocentric latitude is lower
end

@testset "GroundStation — construction and validation" begin
    gs = GroundStation(name = "DSS-14", body = earth, latitude = 35, longitude = -117,
                       altitude = 1, min_elevation = 10)
    @test gs.name == "DSS-14"
    @test gs.body === earth
    @test gs.reference isa Ellipsoid
    @test gs.latitude isa Float64 && gs.latitude == 35.0
    @test gs.longitude == -117.0 && gs.altitude == 1.0
    @test gs.min_elevation_deg == 10.0
    @test GroundStation(name = "x", body = earth, latitude = 0, longitude = 0,
                        altitude = 0).min_elevation_deg == -90.0

    # Latitude at the poles is allowed; beyond is not. The elevation cutoff excludes the zenith.
    @test GroundStation(name = "pole", body = earth, latitude = -90.0, longitude = 0.0,
                        altitude = 0.0) isa GroundStation
    @test_throws ArgumentError GroundStation(name = "x", body = earth, latitude = 90.5,
                                             longitude = 0.0, altitude = 0.0)
    @test_throws ArgumentError GroundStation(name = "x", body = earth, latitude = 0.0,
                                             longitude = 0.0, altitude = 0.0, min_elevation = 90.0)
    @test_throws ArgumentError GroundStation(name = "x", body = earth, latitude = 0.0,
                                             longitude = 0.0, altitude = 0.0, min_elevation = -91.0)
end

@testset "GroundStation — state in body-fixed and inertial axes" begin
    gs = GroundStation(name = "Madrid", body = earth, latitude = 40.43, longitude = -4.25,
                       altitude = 0.83)
    r_bf = geodetic_to_body_fixed(earth, gs.latitude, gs.longitude, gs.altitude, Ellipsoid())

    r_itrf, v_itrf = get_state(gs, _GS_EPOCH; axes = ITRF())
    @test r_itrf ≈ r_bf atol = 1e-6
    @test norm(v_itrf) < 1e-12

    r_gcrf, v_gcrf = get_state(gs, _GS_EPOCH)
    ω = 7.292115e-5                                   # rad/s, mean Earth rotation rate
    @test norm(r_gcrf) ≈ norm(r_bf) rtol = 1e-12
    @test norm(v_gcrf) ≈ ω * hypot(r_bf[1], r_bf[2]) rtol = 1e-4
    @test abs(dot(r_gcrf, v_gcrf)) / (norm(r_gcrf) * norm(v_gcrf)) < 1e-4
    @test r_gcrf != r_bf                                # the axes really did change
end

@testset "GroundStation — visibility either side of the elevation cutoff" begin
    gs = GroundStation(name = "Canberra", body = earth, latitude = -35.40, longitude = 148.98,
                       altitude = 0.69, min_elevation = 10.0)
    r_gs, _ = get_state(gs, _GS_EPOCH)
    up = r_gs / norm(r_gs)
    east = normalize(cross([0.0, 0.0, 1.0], up))       # any direction in the local horizontal

    at_elevation(deg; range = 2000.0) =
        r_gs .+ range .* (cosd(deg) .* east .+ sind(deg) .* up)

    @test is_visible(gs, at_elevation(90.0), _GS_EPOCH)
    @test is_visible(gs, at_elevation(10.5), _GS_EPOCH)
    @test !is_visible(gs, at_elevation(9.5), _GS_EPOCH)
    @test !is_visible(gs, at_elevation(-30.0), _GS_EPOCH)
    @test !is_visible(gs, -r_gs .* 1.1, _GS_EPOCH)       # the far side of the Earth
    @test !is_visible(gs, collect(r_gs), _GS_EPOCH)      # zero range has no direction

    # The default cutoff of -90° sees everything with a direction.
    anywhere = GroundStation(name = "any", body = earth, latitude = -35.40, longitude = 148.98,
                             altitude = 0.69)
    @test is_visible(anywhere, -r_gs .* 1.1, _GS_EPOCH)
end

@testset "GroundStation — elevation is measured from the geodetic vertical" begin
    # At 45° latitude the geodetic and geocentric verticals differ by 0.19°, in the meridian
    # plane. Positions are built in body-fixed axes due north and due south of the station at a
    # stated geodetic elevation, 0.02° either side of the cutoff, which a geocentric vertical
    # gets wrong in one direction or the other. The vertical is the analytic ellipsoid normal
    # (cos φ cos λ, cos φ sin λ, sin φ); the positions go to GCRF through AstroFrames.
    φ, λ = 45.0, 10.0
    gs = GroundStation(name = "mid-latitude", body = earth, latitude = φ, longitude = λ,
                       altitude = 0.0, min_elevation = 10.0)
    r_gs = geodetic_to_body_fixed(earth, φ, λ, 0.0, Ellipsoid())
    up    = [cosd(φ) * cosd(λ), cosd(φ) * sind(λ), sind(φ)]
    north = [-sind(φ) * cosd(λ), -sind(φ) * sind(λ), cosd(φ)]

    function in_gcrf(el_deg, horizontal)
        r_bf = r_gs .+ 1500.0 .* (cosd(el_deg) .* horizontal .+ sind(el_deg) .* up)
        c = Coordinate(vcat(r_bf, zeros(3)), CoordinateSystem(earth, ITRF()), _GS_EPOCH)
        return to_vector(CartesianState(Coordinate(c, CoordinateSystem(earth, GCRF()))))[1:3]
    end

    for horizontal in (north, -north)
        @test  is_visible(gs, in_gcrf(10.02, horizontal), _GS_EPOCH)
        @test !is_visible(gs, in_gcrf(9.98,  horizontal), _GS_EPOCH)
    end

    # The geocentric vertical, for contrast: it disagrees with the geodetic one here.
    geocentric_el(r) = asind(dot(r - r_gs, r_gs) / (norm(r - r_gs) * norm(r_gs)))
    r_n = r_gs .+ 1500.0 .* (cosd(10.02) .* north .+ sind(10.02) .* up)
    @test geocentric_el(r_n) < 10.0
end

@testset "SphericalDrag and SphericalSRP — construction and display" begin
    d = SphericalDrag(c_d = 2, drag_area = 3.5)
    @test d isa AbstractDragGeometry
    @test d.c_d === 2.0 && d.drag_area === 3.5
    @test sprint(show, d) == "SphericalDrag(c_d = 2.0, drag_area = 3.5 m²)"
    @test sprint(show, MIME"text/plain"(), d) == "SphericalDrag:\n  c_d  = 2.0\n  area = 3.5 m²"

    s = SphericalSRP(c_r = 1.8, srp_area = 10)
    @test s isa AbstractSRPGeometry
    @test s.c_r === 1.8 && s.srp_area === 10.0
    @test sprint(show, s) == "SphericalSRP(c_r = 1.8, srp_area = 10.0 m²)"
    @test sprint(show, MIME"text/plain"(), s) == "SphericalSRP:\n  c_r  = 1.8\n  area = 10.0 m²"

    sc = Spacecraft(drag = d, srp = s)
    @test sc.drag === d && sc.srp === s
end

@testset "CartesianState(sc, cs) — a spacecraft's state in another frame" begin
    pv = [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9]
    sc = Spacecraft(state = CartesianState(pv), time = _GS_EPOCH)

    @test AstroFrames.state_of(sc) === sc.state
    @test AstroFrames.frame_of(sc) === sc.coord_sys
    @test AstroFrames.epoch_of(sc) === sc.time

    # Same frame: the state comes back unchanged.
    same = CartesianState(sc, sc.coord_sys)
    @test to_vector(same) ≈ pv rtol = 1e-12

    # To Earth-fixed and back through a second spacecraft: a rotation, so |r| is kept, and the
    # round trip recovers the state.
    itrf  = CoordinateSystem(earth, ITRF())
    fixed = CartesianState(sc, itrf)
    @test norm(to_vector(fixed)[1:3]) ≈ norm(pv[1:3]) rtol = 1e-12
    back_sc = Spacecraft(state = fixed, time = _GS_EPOCH, coord_sys = itrf)
    @test to_vector(CartesianState(back_sc, sc.coord_sys)) ≈ pv rtol = 1e-9

    # The params form reaches the same implementation.
    @test to_vector(CartesianState(sc, itrf, NamedTuple())) == to_vector(fixed)
end

nothing
