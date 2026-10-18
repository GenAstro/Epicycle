# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Named shortcuts, body shape and rotation, and Sun and Moon geometry (#11, sections 2–4).
#
# Truth, strongest first:
# - External reference: GMAT R2026a, from truth/body_parameters_truth.script, for the apsis
#   radii, semi-latus rectum, C3, beta angle, and on the Earth the geodetic altitude, latitude
#   and longitude, sidereal time and local sidereal time. Three orbits at 20 Mar 2024 12:00 UTC.
# - Independent inverse: AstroModels' closed-form geodetic-to-body-fixed conversion, run back
#   through the geodetic readers, on Mars, whose ellipsoid is oblate.
# - Geometry with a known answer: due-north and due-east velocities, a subsolar point, a
#   spacecraft behind the Earth or the Moon, the angles of a triangle summing to π.
# - Differentiability: ForwardDiff against central differences.
# =============================================================================

using Test
using ForwardDiff
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using EpicycleBase: label, unit, is_cyclic, is_settable, set_quantity!
using LinearAlgebra: cross, dot, norm, normalize

const _BP_EPOCH = Time("2024-03-20T12:00:00", UTC(), ISOT())
const _BP_JD    = _BP_EPOCH.tdb.jd

const _BP_STATES = Dict(
    :leo     => [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9],
    :molniya => [4637.031328726552, 178.53697947901998, -5679.055240387162,
                 6.252424682730313, 6.928411997008257, 2.5730558589825416],
    :polar   => [1731.6965430023222, -2999.386395771405, -7427.270441381914,
                 -3.0611710106575205, 5.302103721115798, -2.5783940176845728],
    :hyp     => [-3701.281373968369, 7485.70648182266, 4011.578013161604,
                 -9.268115532432553, -1.7246408988602304, 4.13619682261008])

_bp_sat(name) = Spacecraft(state = CartesianState(_BP_STATES[name]), time = _BP_EPOCH,
                           coord_sys = EarthMJ2000Eq)

# A lunar orbiter and a Mars-fixed point: body-fixed frames that need no Earth orientation data.
_bp_lunar() = Spacecraft(state = CartesianState([1900.0, 100.0, 300.0, -0.1, 1.6, 0.3]),
                         time = _BP_EPOCH, coord_sys = CoordinateSystem(moon, ICRF()))
const _BP_MARS_FIXED = CoordinateSystem(mars, CelestialBodyFixed())

# GMAT output, angles in degrees.
const _BP_TRUTH = Dict(
    :leo     => (ra = 7144.525448883372, rp = 6985.629545395033, p = 7064.184091930428,
                 c3 = -56.4184104082937, beta = -9.327968910882271, alt = 730.9375422751282,
                 lat = 9.911516037658375, lon = 4.253356622782521, mha = 358.5115570961498,
                 lst = 2.764913718932291),
    :molniya => (ra = 46283.99999999989, rp = 6915.999999999998, p = 12033.84,
                 c3 = -14.98497901503762, beta = 35.06184334952253, alt = 968.5250417859525,
                 lat = -50.77646507872066, lon = 4.000170227630063, mha = 358.5115570961498,
                 lst = 2.511727323779837),
    :polar   => (ra = 8250.0, rp = 6750.0, p = 7424.999999999999, c3 = -53.14672557333333,
                 beta = 60.03377086877112, alt = 1834.506842108826, lat = -65.04905050163741,
                 lon = -57.94996266998206, mha = 358.5115570961498, lst = 300.5615944261677))

_bp_angle_deg(a, b) = abs(rem(a - b + 180, 360) - 180)      # difference across the wrap

@testset "agree with an external reference (GMAT): inertial" begin
    for name in (:leo, :molniya, :polar)
        @testset "$name" begin
            sat, t = _bp_sat(name), _BP_TRUTH[name]
            @test apoapsis_radius(sat)  ≈ t.ra rtol = 1e-12
            @test periapsis_radius(sat) ≈ t.rp rtol = 1e-12
            @test semilatus_rectum(sat) ≈ t.p  rtol = 1e-12
            @test c3(sat)               ≈ t.c3 rtol = 1e-12
            # The reference reads the Sun from its default DE405 ephemeris; the planetary kernel
            # here is newer. The two place the Sun about 13 km apart, 5e-6° as seen from Earth.
            @test rad2deg(beta_angle(sat)) ≈ t.beta atol = 2e-5
        end
    end
end

@testset "agree with an external reference (GMAT): Earth-fixed" begin
    for name in (:leo, :molniya, :polar)
        @testset "$name" begin
            sat, t = _bp_sat(name), _BP_TRUTH[name]
            @test altitude(sat) ≈ t.alt atol = 1e-3                          # km
            @test rad2deg(geodetic_latitude(sat)) ≈ t.lat atol = 1e-5
            @test _bp_angle_deg(rad2deg(longitude(sat)), t.lon) < 1e-4
            @test _bp_angle_deg(rad2deg(meridian_hour_angle(sat)), t.mha) < 1e-4
            @test _bp_angle_deg(rad2deg(local_sidereal_time(sat)), t.lst) < 1e-4
        end
    end
end

@testset "named shortcuts satisfy the conic identities" begin
    for name in (:leo, :molniya, :polar, :hyp)
        @testset "$name" begin
            sat = _bp_sat(name)
            a, e = semi_major_axis(sat), eccentricity(sat)
            @test periapsis_radius(sat) ≈ a * (1 - e) rtol = 1e-12
            @test semilatus_rectum(sat) ≈ a * (1 - e^2) rtol = 1e-12
            @test c3(sat) ≈ 2 * specific_energy(sat) rtol = 1e-14
            @test true_longitude(sat) ≈
                  mod(raan(sat) + argument_of_periapsis(sat) + true_anomaly(sat), 2π) atol = 1e-12
            e < 1 && @test apoapsis_radius(sat) ≈ a * (1 + e) rtol = 1e-12
        end
    end
    @test hyperbolic_excess_speed(_bp_sat(:hyp))^2 ≈ c3(_bp_sat(:hyp)) rtol = 1e-14
    @test_throws ArgumentError apoapsis_radius(_bp_sat(:hyp))
end

@testset "apoapsis_radius! and periapsis_radius! hold the other radius" begin
    sat = _bp_sat(:molniya)
    rp0, inc0 = periapsis_radius(sat), inclination(sat)
    apoapsis_radius!(sat; to = 40000.0)
    @test apoapsis_radius(sat) ≈ 40000.0 rtol = 1e-10
    @test periapsis_radius(sat) ≈ rp0 rtol = 1e-10
    @test inclination(sat) ≈ inc0 rtol = 1e-12

    ra0 = apoapsis_radius(sat)
    set_quantity!(sat, periapsis_radius; to = 7100.0)
    @test periapsis_radius(sat) ≈ 7100.0 rtol = 1e-10
    @test apoapsis_radius(sat) ≈ ra0 rtol = 1e-10

    @test is_settable(sat, apoapsis_radius) && !is_settable(sat, c3)
    before = position_vector(sat)
    @test_throws ArgumentError periapsis_radius!(sat; to = -5.0)
    @test position_vector(sat) ≈ before
end

@testset "geodetic readers invert the closed-form forward conversion" begin
    # Mars is oblate, so geodetic and geocentric latitude differ and the iteration is exercised.
    @test mars.flattening > 0
    for (lat, lon, alt) in ((0.0, 0.0, 400.0), (35.0, 120.0, 250.0), (-62.5, -45.0, 1200.0),
                            (89.9, 10.0, 300.0), (-10.0, 179.0, 30000.0))
        r = AstroModels.geodetic_to_body_fixed(mars, lat, lon, alt, Ellipsoid())
        c = Coordinate([collect(r); 0.0; 0.0; 0.0], _BP_MARS_FIXED, _BP_EPOCH)
        @test altitude(c, mars) ≈ alt atol = 1e-6
        @test rad2deg(geodetic_latitude(c, mars)) ≈ lat atol = 1e-9
        @test rad2deg(longitude(c, mars)) ≈ lon atol = 1e-9
    end
end

@testset "velocity azimuth is measured clockwise from local north" begin
    lat, lon = deg2rad(30.0), deg2rad(50.0)
    r = AstroModels.geodetic_to_body_fixed(mars, 30.0, 50.0, 500.0, Ellipsoid())
    north = [-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)]
    east  = [-sin(lon), cos(lon), 0.0]
    for (v, az) in ((north, 0.0), (east, π / 2), (-north, π), (-east, 3π / 2),
                    (normalize(north + east), π / 4))
        c = Coordinate([collect(r); 2.0 .* v], _BP_MARS_FIXED, _BP_EPOCH)
        Δ = rem(velocity_azimuth(c, mars) - az, 2π, RoundNearest)    # 0 and 2π are the same heading
        @test abs(Δ) < 1e-10
    end
end

@testset "apsis altitudes on a spherical body are radius less the radius of the body" begin
    sat = _bp_lunar()
    @test moon.flattening == 0
    @test apoapsis_altitude(sat) ≈ apoapsis_radius(sat) - moon.equatorial_radius rtol = 1e-10
    @test periapsis_altitude(sat) ≈ periapsis_radius(sat) - moon.equatorial_radius rtol = 1e-10
    @test periapsis_altitude(sat) <= altitude(sat) <= apoapsis_altitude(sat)
end

@testset "apsis altitudes on an oblate body sit below the radius less the equatorial radius" begin
    # Off the equator the ellipsoid is closer to the centre, so the height is greater.
    sat = _bp_sat(:molniya)
    @test periapsis_altitude(sat) > periapsis_radius(sat) - earth.equatorial_radius
    @test periapsis_altitude(sat) - (periapsis_radius(sat) - earth.equatorial_radius) < 21.4
    @test_throws ArgumentError apoapsis_altitude(_bp_sat(:hyp))
end

@testset "sidereal time and geographic longitude" begin
    sat = _bp_lunar()
    @test local_sidereal_time(sat) ≈ mod(meridian_hour_angle(sat) + longitude(sat), 2π) atol = 1e-12

    # A geostationary spacecraft above 75° E reads a mean geographic longitude of 75° E.
    λ = deg2rad(75.0)
    rgeo = 42164.17
    ω = 7.2921150e-5
    c = Coordinate([rgeo * cos(λ), rgeo * sin(λ), 0.0, 0.0, 0.0, 0.0], EarthFixed, _BP_EPOCH)
    @test rad2deg(mean_geographic_longitude(c)) ≈ 75.0 atol = 0.05
    @test rad2deg(longitude(c)) ≈ 75.0 atol = 1e-9
end

@testset "away from the Earth the hour angle is the IAU prime-meridian angle W" begin
    # Truth: the IAU 2015 rotation model, whose W is accumulated in radians.
    for body in (mars, jupiter, venus)
        c = Coordinate([4000.0, 300.0, 1200.0, -0.4, 3.4, 0.9], CoordinateSystem(body, ICRF()), _BP_EPOCH)
        W = mod(iau2015_orientation(body, _BP_JD).prime_meridian, 2π)
        @test abs(rem(meridian_hour_angle(c, body) - W, 2π, RoundNearest)) < 1e-9
    end
end

@testset "mean geographic longitude is measured in the body's own equator" begin
    # An orbit in Mars's equator, 25° off the ICRF equator. There, mean geographic longitude is
    # the longitude below the spacecraft advanced by M − ν, with no plane mixing in between.
    c0  = Coordinate([4000.0, 300.0, 1200.0, -0.4, 3.4, 0.9], CoordinateSystem(mars, ICRF()), _BP_EPOCH)
    x̂, ŷ, _ = AstroCallbacks._equator_axes(c0, mars)
    a, e, ν = 5000.0, 0.05, 0.7
    p   = a * (1 - e^2)
    μ   = gravitational_parameter(mars)
    r   = p / (1 + e * cos(ν))
    pos = r * (cos(ν) .* x̂ .+ sin(ν) .* ŷ)
    vel = sqrt(μ / p) * (-sin(ν) .* x̂ .+ (e + cos(ν)) .* ŷ)
    c   = Coordinate([pos; vel], CoordinateSystem(mars, ICRF()), _BP_EPOCH)
    E   = 2atan(sqrt((1 - e) / (1 + e)) * tan(ν / 2))
    M   = E - e * sin(E)
    @test abs(rem(mean_geographic_longitude(c, mars) - (longitude(c, mars) + M - ν), 2π,
                  RoundNearest)) < 1e-9
end

@testset "nodal precession: a sun-synchronous orbit precesses at the mean solar rate" begin
    # 700 km circular at the sun-synchronous inclination for that altitude, 98.19°.
    a, i = 7078.137, deg2rad(98.186)
    v = sqrt(398600.4418 / a)
    x = [a, 0.0, 0.0, 0.0, v * cos(i), v * sin(i)]
    sat = Spacecraft(state = CartesianState(x), time = _BP_EPOCH, coord_sys = EarthMJ2000Eq)
    @test nodal_precession_rate(sat) ≈ 2π / (365.2422 * 86400) rtol = 2e-3
    # A prograde orbit regresses.
    @test nodal_precession_rate(_bp_sat(:leo)) < 0
    @test_throws ArgumentError nodal_precession_rate(_bp_lunar())        # no default J2 for the Moon
    @test nodal_precession_rate(_bp_lunar(); j2 = 2.03e-4) < 0
end

@testset "Sun and Moon angles" begin
    sat  = _bp_sat(:leo)
    r    = position_vector(sat, EarthICRF)          # the ephemeris is in ICRF axes
    s    = translate(earth, sun,  _BP_JD)
    m    = translate(earth, moon, _BP_JD)
    ang(a, b) = acos(clamp(dot(a, b) / (norm(a) * norm(b)), -1, 1))

    @test sun_earth_spacecraft_angle(sat)  ≈ ang(s, r)           atol = 1e-12
    @test sun_spacecraft_earth_angle(sat)  ≈ ang(s - r, -r)       atol = 1e-12
    @test sun_spacecraft_moon_angle(sat)   ≈ ang(s - r, m - r)    atol = 1e-12
    @test moon_earth_sun_angle(sat)        ≈ ang(m, s)            atol = 1e-12
    @test moon_earth_spacecraft_angle(sat) ≈ ang(m, r)            atol = 1e-12
    @test moon_spacecraft_earth_angle(sat) ≈ ang(m - r, -r)       atol = 1e-12

    # The angles of the Sun–Earth–spacecraft triangle sum to π.
    @test sun_earth_spacecraft_angle(sat) + sun_spacecraft_earth_angle(sat) + ang(-s, r - s) ≈ π

    # From near the Earth the Sun lies along almost the same line as from its centre.
    h = cross(r, velocity_vector(sat, EarthICRF))
    @test lunar_beta_angle(sat) ≈ asin(dot(h, m) / (norm(h) * norm(m))) atol = 1e-12
    @test orbit_normal_sun_angle(sat) ≈ π / 2 - beta_angle(sat) atol = 2 * norm(r) / norm(s)
    @test_throws ArgumentError lunar_beta_angle(_bp_lunar())
end

@testset "shadow" begin
    ŝ = normalize(translate(earth, sun, _BP_JD))
    behind = Coordinate([-7000.0 .* ŝ; 0.0; 0.0; 0.0], EarthMJ2000Eq, _BP_EPOCH)
    sunlit = Coordinate([ 7000.0 .* ŝ; 0.0; 0.0; 0.0], EarthMJ2000Eq, _BP_EPOCH)
    @test earth_shadow_fraction(behind) == 1
    @test earth_shadow_fraction(sunlit) == 0

    # Crossing the shadow edge sideways, the hidden fraction rises continuously from 0 to 1.
    side = normalize(cross(ŝ, [0.0, 0.0, 1.0]))
    fracs = map(range(6300.0, 6460.0; length = 81)) do y
        earth_shadow_fraction(Coordinate([-7000.0 .* ŝ .+ y .* side; 0.0; 0.0; 0.0],
                                         EarthMJ2000Eq, _BP_EPOCH))
    end
    @test fracs[1] == 1 && fracs[end] == 0
    @test issorted(fracs; rev = true)
    @test any(f -> 0 < f < 1, fracs)

    # 2000 km behind the Moon, in line with the Sun: inside its umbra, which is 374,000 km long.
    m  = translate(earth, moon, _BP_JD)
    ŝm = normalize(translate(moon, sun, _BP_JD))
    lunar_shadow = Coordinate([m .- 2000.0 .* ŝm; 0.0; 0.0; 0.0], EarthMJ2000Eq, _BP_EPOCH)
    @test moon_shadow_fraction(lunar_shadow) == 1
    @test moon_shadow_fraction(_bp_sat(:leo)) == 0
end

@testset "local time" begin
    # Below the Sun it is noon, opposite it midnight.
    s = translate(moon, sun, _BP_JD)
    lunar_cs = CoordinateSystem(moon, ICRF())
    noon     = Coordinate([ 1900.0 .* normalize(s); 0.0; 1.6; 0.0], lunar_cs, _BP_EPOCH)
    midnight = Coordinate([-1900.0 .* normalize(s); 0.0; 1.6; 0.0], lunar_cs, _BP_EPOCH)
    @test true_local_time(noon) ≈ 12.0 atol = 1e-6
    @test min(true_local_time(midnight), 24 - true_local_time(midnight)) < 1e-6

    # On the Earth, mean and true local time differ by the equation of time, under 17 minutes,
    # and a point on the Greenwich meridian reads UTC.
    sat = _bp_sat(:leo)
    Δ = abs(rem(mean_local_time(sat) - true_local_time(sat) + 12, 24) - 12)
    @test Δ < 17 / 60
    greenwich = Coordinate([7000.0, 0.0, 0.0, 0.0, 0.0, 7.5], EarthFixed, _BP_EPOCH)
    @test mean_local_time(greenwich) ≈ 12.0 atol = 1e-9
end

# --- differentiability --------------------------------------------------------

const _BP_H = [1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6]

function _bp_central(f, x)
    map(eachindex(x)) do j
        δ = zeros(length(x)); δ[j] = _BP_H[j]
        (f(x .+ δ) - f(x .- δ)) / (2 * _BP_H[j])
    end
end

@testset "differentiable: ForwardDiff matches central differences" begin
    xe = _BP_STATES[:molniya]
    earth_sub(x) = Coordinate(x, EarthMJ2000Eq, _BP_EPOCH)
    for q in (apoapsis_radius, periapsis_radius, semilatus_rectum, c3, true_longitude, beta_angle,
              lunar_beta_angle, orbit_normal_sun_angle, sun_earth_spacecraft_angle,
              sun_spacecraft_moon_angle, moon_spacecraft_earth_angle)
        @testset "$(nameof(q))" begin
            f = x -> q(earth_sub(x))
            @test ForwardDiff.gradient(f, xe) ≈ _bp_central(f, xe) rtol = 1e-6
        end
    end
    xl = [1900.0, 100.0, 300.0, -0.1, 1.6, 0.3]
    lunar_sub(x) = Coordinate(x, CoordinateSystem(moon, ICRF()), _BP_EPOCH)
    for q in (altitude, geodetic_latitude, longitude, velocity_azimuth, periapsis_altitude,
              apoapsis_altitude, local_sidereal_time, mean_geographic_longitude, true_local_time)
        @testset "$(nameof(q)), lunar orbit" begin
            f = x -> q(lunar_sub(x))
            @test ForwardDiff.gradient(f, xl) ≈ _bp_central(f, xl) rtol = 1e-6
        end
    end
end

@testset "labels, units and wrapping" begin
    @test label(altitude) == "Altitude" && unit(altitude) == "km"
    @test label(beta_angle) == "Beta angle" && unit(beta_angle) == "rad"
    @test unit(true_local_time) == "h" && unit(nodal_precession_rate) == "rad/s"
    @test unit(earth_shadow_fraction) == ""
    @test is_cyclic(longitude) && is_cyclic(true_local_time)
    @test !is_cyclic(altitude) && !is_cyclic(beta_angle)
    @test beta_angle in quantities(Spacecraft) && earth_shadow_fraction in quantities(Coordinate)
end

nothing
