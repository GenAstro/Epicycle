# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Named shortcuts, and Sun and Moon geometry (#11, sections 2 and 4).
#
# Truth, strongest first:
# - External reference: GMAT R2026a, from truth/body_parameters_truth.script, for the apsis
#   radii, semi-latus rectum, C3 and beta angle. Three orbits at 20 Mar 2024 12:00 UTC, with the
#   Sun from the same DE440 kernel in both tools. The beta angle's tolerance is loose because
#   the two tools' frame models differ; it confirms the definition, and the self-consistency
#   test (test_correctness_sun_moon_frames.jl) is where its accuracy and precision are tested.
# - Geometry with a known answer: a spacecraft behind the Earth or the Moon, the angles of a
#   triangle summing to π.
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

# A lunar orbiter, which has no lunar beta angle.
_bp_lunar() = Spacecraft(state = CartesianState([1900.0, 100.0, 300.0, -0.1, 1.6, 0.3]),
                         time = _BP_EPOCH, coord_sys = CoordinateSystem(moon, ICRF()))

# GMAT output, angles in degrees.
const _BP_TRUTH = Dict(
    :leo     => (ra = 7144.525448883372, rp = 6985.629545395033, p = 7064.184091930428,
                 c3 = -56.4184104082937, beta = -9.327980902878894),
    :molniya => (ra = 46283.99999999989, rp = 6915.999999999998, p = 12033.84,
                 c3 = -14.98497901503762, beta = 35.06183703901353),
    :polar   => (ra = 8250.0, rp = 6750.0, p = 7424.999999999999, c3 = -53.14672557333333,
                 beta = 60.03377047858275))

@testset "agree with an external reference (GMAT): inertial" begin
    for name in (:leo, :molniya, :polar)
        @testset "$name" begin
            sat, t = _bp_sat(name), _BP_TRUTH[name]
            @test apoapsis_radius(sat)  ≈ t.ra rtol = 1e-12
            @test periapsis_radius(sat) ≈ t.rp rtol = 1e-12
            @test semilatus_rectum(sat) ≈ t.p  rtol = 1e-12
            @test c3(sat)               ≈ t.c3 rtol = 1e-12
            # The tolerance is loose on purpose. Both tools read the Sun from the same kernel, but
            # their frame models differ: how each ties the ephemeris frame to MJ2000Eq leaves a
            # residual of up to 7e-6° (0.025″), the size of the ICRF–MJ2000Eq frame bias. This
            # check confirms the definition matches GMAT's. Accuracy and precision are tested in
            # test_correctness_sun_moon_frames.jl, against Epicycle's own models, to 1e-12.
            @test rad2deg(beta_angle(sat)) ≈ t.beta atol = 2e-5
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
end

@testset "labels, units and wrapping" begin
    @test label(beta_angle) == "Beta angle" && unit(beta_angle) == "rad"
    @test label(earth_shadow_fraction) == "Earth shadow fraction" && unit(earth_shadow_fraction) == ""
    @test !is_cyclic(beta_angle)
    @test beta_angle in quantities(Spacecraft) && earth_shadow_fraction in quantities(Coordinate)
end

nothing
