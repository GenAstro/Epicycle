# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Each orbit parameter reads the frame, origin and gravitational parameter it is asked for.
#
# The math is checked elsewhere: against GMAT in test_correctness_orbit_parameters.jl, and in
# AstroRoutines and AstroStates. What is checked here is the wiring. A spacecraft is created
# from a state, its state is converted to a target frame with `Coordinate`, the underlying
# routine or state representation is applied to that converted state with the target origin's
# gravitational parameter, and the quantity read in the target frame must give the same value,
# or throw the same ArgumentError where the value is undefined.
#
# Frames: Earth MJ2000Eq, MJ2000Ec, ICRF and body-fixed (a rotating frame, where the velocity
# changes); Moon-centred MJ2000Eq and body-fixed; Mars-centred MJ2000Eq and body-fixed. The frame
# theory is IAU2006, Epicycle's default, set explicitly so the checks do not depend on what ran
# before. Each subject is both a Spacecraft and a Coordinate with its own epoch. The setters are
# checked the same way, in a lunar frame.
# =============================================================================

using Test
using AstroCallbacks
using AstroFrames
using AstroModels
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroRoutines

const _PF_EPOCH = Time("2024-01-01T12:00:00", UTC(), ISOT())
const _PF_ORIGINAL_THEORY = frame_theory()
set_frame_theory!(IAU2006())

# An orbit about each body, in that body's MJ2000Eq frame, and the frames to read it in.
const _PF_CASES = (
    (name = "Molniya about the Earth", x = [4637.031328726552, 178.53697947901998, -5679.055240387162,
                                            6.252424682730313, 6.928411997008257, 2.5730558589825416],
     from = CoordinateSystem(earth, MJ2000Eq()),
     frames = (EarthMJ2000Eq, EarthMJ2000Ec, EarthICRF, EarthFixed)),
    (name = "hyperbola about the Earth", x = [-3701.281373968369, 7485.70648182266, 4011.578013161604,
                                              -9.268115532432553, -1.7246408988602304, 4.13619682261008],
     from = CoordinateSystem(earth, MJ2000Eq()),
     frames = (EarthMJ2000Eq, EarthMJ2000Ec, EarthFixed)),
    (name = "lunar orbit", x = [1900.0, 300.0, -200.0, 0.2, 1.55, 0.4],
     from = CoordinateSystem(moon, MJ2000Eq()),
     frames = (CoordinateSystem(moon, MJ2000Eq()), MoonFixed)),
    (name = "Mars orbit", x = [4000.0, -1500.0, 900.0, 1.0, 2.9, 1.1],
     from = CoordinateSystem(mars, MJ2000Eq()),
     frames = (CoordinateSystem(mars, MJ2000Eq()), CoordinateSystem(mars, CelestialBodyFixed()))),
    (name = "Mars flyby", x = [-8000.0, 3000.0, 1500.0, 2.5, -4.0, 1.0],
     from = CoordinateSystem(mars, MJ2000Eq()),
     frames = (CoordinateSystem(mars, MJ2000Eq()), CoordinateSystem(mars, CelestialBodyFixed()))))

_pf_subjects(case) = (Spacecraft(state = CartesianState(case.x), time = _PF_EPOCH, coord_sys = case.from),
                      Coordinate(case.x, case.from, _PF_EPOCH))

# --- the underlying math, applied to a state already in the target frame ----------------------

_pf_wrap(x) = mod(x, 2π)

function _pf_b_plane(x, μ)
    asymptote = IncomingAsymptoteState(CartesianState(x), μ)
    e = 1 + asymptote.rp * asymptote.c3 / μ
    e > 1 || throw(ArgumentError("not a hyperbola"))
    magnitude = sqrt(asymptote.rp^2 + 2μ * asymptote.rp / asymptote.c3)
    return magnitude, magnitude * sin(asymptote.bpa), magnitude * cos(asymptote.bpa)
end

function _pf_mean_anomaly(x, μ)
    kep = KeplerianState(CartesianState(x), μ)
    kep.ecc < 1 && return _pf_wrap(true_to_mean_anomaly(kep.ta, kep.ecc))
    return hyperbolic_to_mean_anomaly(true_to_hyperbolic_anomaly(kep.ta, kep.ecc), kep.ecc)
end

const _PF_EXPECTED = (
    mean_anomaly                    => _pf_mean_anomaly,
    eccentric_anomaly               => (x, μ) -> (k = KeplerianState(CartesianState(x), μ);
                                                  _pf_wrap(true_to_eccentric_anomaly(k.ta, k.ecc))),
    hyperbolic_anomaly              => (x, μ) -> (k = KeplerianState(CartesianState(x), μ);
                                                  true_to_hyperbolic_anomaly(k.ta, k.ecc)),
    specific_energy                 => (x, μ) -> specific_energy(x[1:3], x[4:6], μ),
    hyperbolic_excess_speed         => (x, μ) -> hyperbolic_excess_speed(x[1:3], x[4:6], μ),
    apoapsis_speed                  => (x, μ) -> apoapsis_speed(x[1:3], x[4:6], μ),
    periapsis_speed                 => (x, μ) -> periapsis_speed(x[1:3], x[4:6], μ),
    angular_momentum                => (x, μ) -> angular_momentum(x[1:3], x[4:6]),
    angular_momentum_magnitude      => (x, μ) -> angular_momentum_magnitude(x[1:3], x[4:6]),
    orbit_normal                    => (x, μ) -> orbit_normal(x[1:3], x[4:6]),
    orbital_angular_velocity        => (x, μ) -> orbital_angular_velocity(x[1:3], x[4:6]),
    eccentricity_vector             => (x, μ) -> eccentricity_vector(x[1:3], x[4:6], μ),
    orbital_period                  => (x, μ) -> orbital_period(x[1:3], x[4:6], μ),
    mean_motion                     => (x, μ) -> mean_motion(x[1:3], x[4:6], μ),
    argument_of_latitude            => (x, μ) -> (k = KeplerianState(CartesianState(x), μ); _pf_wrap(k.aop + k.ta)),
    longitude_of_periapsis          => (x, μ) -> (k = KeplerianState(CartesianState(x), μ); _pf_wrap(k.raan + k.aop)),
    descending_node_right_ascension => (x, μ) -> (k = KeplerianState(CartesianState(x), μ); _pf_wrap(k.raan + π)),
    apoapsis_radius                 => (x, μ) -> (m = ModifiedKeplerianState(CartesianState(x), μ);
                                                  m.ra > 0 ? m.ra : throw(ArgumentError("open orbit"))),
    periapsis_radius                => (x, μ) -> ModifiedKeplerianState(CartesianState(x), μ).rp,
    semilatus_rectum                => (x, μ) -> ModifiedEquinoctialState(CartesianState(x), μ).p,
    c3                              => (x, μ) -> 2 * specific_energy(x[1:3], x[4:6], μ),
    true_longitude                  => (x, μ) -> _pf_wrap(ModifiedEquinoctialState(CartesianState(x), μ).L),
    b_dot_t                         => (x, μ) -> _pf_b_plane(x, μ)[2],
    b_dot_r                         => (x, μ) -> _pf_b_plane(x, μ)[3],
    b_magnitude                     => (x, μ) -> _pf_b_plane(x, μ)[1],
    b_plane_angle                   => (x, μ) -> (b = _pf_b_plane(x, μ); _pf_wrap(atan(b[3], b[2]))))

# The value, or :throws where the math says the quantity is undefined.
_pf_value(f) = try
    f()
catch err
    err isa ArgumentError || rethrow()
    :throws
end

_pf_same(a, b) = a === :throws || b === :throws ? a === b :
                 isapprox(a, b; rtol = 1e-12, atol = 1e-12 * max(1, maximum(abs, a)))

@testset "every orbit parameter reads the frame, origin and μ it is asked for" begin
    for case in _PF_CASES, frame in case.frames, subject in _pf_subjects(case)
        @testset "$(case.name) in $(nameof(typeof(frame.axes))) about $(frame.origin.name), $(nameof(typeof(subject)))" begin
            x = to_vector(CartesianState(Coordinate(subject, frame)))
            μ = get_gravparam(frame.origin)
            for (quantity, math) in _PF_EXPECTED
                expected = _pf_value(() -> math(x, μ))
                got      = _pf_value(() -> quantity(subject, frame))
                @test _pf_same(got, expected) || (@info "$(nameof(quantity))" got expected; false)
            end
        end
    end
end

@testset "the frame changes the answer where it should" begin
    # Guards the test above against a frame argument that is ignored: if every frame gave the
    # same number, the comparison would still pass.
    sat = first(_pf_subjects(_PF_CASES[1]))
    @test !(angular_momentum(sat, EarthMJ2000Ec) ≈ angular_momentum(sat, EarthMJ2000Eq))
    @test !(true_longitude(sat, EarthMJ2000Ec) ≈ true_longitude(sat, EarthMJ2000Eq))
    @test !(specific_energy(sat, EarthFixed) ≈ specific_energy(sat, EarthMJ2000Eq))
    lunar = first(_pf_subjects(_PF_CASES[3]))
    @test !(specific_energy(lunar, EarthMJ2000Eq) ≈ specific_energy(lunar, CoordinateSystem(moon, MJ2000Eq())))
end

@testset "setters write in the frame they are given, about the Moon" begin
    lunar_cs = CoordinateSystem(moon, MJ2000Eq())
    μ = get_gravparam(moon)
    radii(sat) = ModifiedKeplerianState(CartesianState(Coordinate(sat, lunar_cs)), μ)

    sat = first(_pf_subjects(_PF_CASES[3]))
    before = radii(sat)
    apoapsis_radius!(sat, lunar_cs; to = 2600.0)
    after = radii(sat)
    @test after.ra ≈ 2600.0 rtol = 1e-12
    @test [after.rp, after.inc, after.raan, after.aop, after.ta] ≈
          [before.rp, before.inc, before.raan, before.aop, before.ta] rtol = 1e-10

    sat = first(_pf_subjects(_PF_CASES[3]))
    periapsis_radius!(sat, lunar_cs; to = 1800.0)
    @test radii(sat).rp ≈ 1800.0 rtol = 1e-12

    sat = first(_pf_subjects(_PF_CASES[3]))
    kep_before = KeplerianState(CartesianState(Coordinate(sat, lunar_cs)), μ)
    mean_anomaly!(sat, lunar_cs; to = 2.0)
    kep_after = KeplerianState(CartesianState(Coordinate(sat, lunar_cs)), μ)
    @test mean_anomaly(sat, lunar_cs) ≈ 2.0 rtol = 1e-12
    @test to_vector(kep_after)[1:5] ≈ to_vector(kep_before)[1:5] rtol = 1e-10
end

set_frame_theory!(_PF_ORIGINAL_THEORY)
