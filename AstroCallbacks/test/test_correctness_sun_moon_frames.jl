# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# The Sun and Moon quantities, rebuilt from Epicycle's own models (#11, section 4).
#
# Each quantity is a combination of a few pieces: the subject's position and velocity about a
# body in ICRF axes, from the frame conversion (`Coordinate`); the Sun and Moon positions about
# that body, from the planetary ephemeris (`translate`) at the subject's TDB epoch; and the
# orbit normal, from AstroRoutines. Here the pieces are computed directly with those models and
# combined by each quantity's definition, and the quantity must agree. GMAT checks the beta
# angle externally in test_correctness_body_parameters.jl; this checks the wiring, for subjects
# defined in inertial, rotating and other-body frames, as a Spacecraft and as a Coordinate, and
# for the beta angles about a body given explicitly. The frame theory is IAU2006, set here.
# =============================================================================

using Test
using LinearAlgebra: dot, norm, cross
using AstroCallbacks
using AstroFrames
using AstroModels
using AstroStates
using AstroEpochs
using AstroUniverse
using AstroRoutines: orbit_normal

const _SM_ORIGINAL_THEORY = frame_theory()
set_frame_theory!(IAU2006())

# Subjects defined in a variety of frames, each about the body it orbits, at different epochs.
const _SM_CASES = (
    (name = "LEO in EarthMJ2000Eq", body = earth, cs = EarthMJ2000Eq,
     x = [7000.0, 300.0, 1200.0, -0.4, 7.4, 0.9], epoch = Time("2024-03-20T12:00:00", UTC(), ISOT())),
    (name = "LEO in EarthFixed", body = earth, cs = EarthFixed,
     x = [6500.0, -2200.0, 1500.0, 1.9, 6.6, 2.4], epoch = Time("2025-07-04T03:30:00", UTC(), ISOT())),
    (name = "lunar orbit in MoonFixed", body = moon, cs = MoonFixed,
     x = [1900.0, 100.0, 300.0, -0.1, 1.6, 0.3], epoch = Time("2024-11-15T18:00:00", UTC(), ISOT())),
    (name = "Mars orbit in Mars MJ2000Eq", body = mars, cs = CoordinateSystem(mars, MJ2000Eq()),
     x = [4000.0, -1500.0, 900.0, 1.0, 2.9, 1.1], epoch = Time("2026-01-10T00:00:00", UTC(), ISOT())))

_sm_subjects(c) = (Spacecraft(state = CartesianState(c.x), time = c.epoch, coord_sys = c.cs),
                   Coordinate(c.x, c.cs, c.epoch))

# --- the pieces, from Epicycle's models ------------------------------------------------------

"Position and velocity of `subject` about `body`, ICRF axes, from the frame conversion."
function _sm_rv(subject, body)
    x = to_vector(CartesianState(Coordinate(subject, CoordinateSystem(body, ICRF()))))
    return x[1:3], x[4:6]
end

_sm_jd(subject) = epoch_of(subject).tdb.jd
_sm_sun(body, subject)  = translate(body, sun,  _sm_jd(subject))
_sm_moon(body, subject) = translate(body, moon, _sm_jd(subject))
_sm_angle(a, b) = acos(clamp(dot(a, b) / (norm(a) * norm(b)), -1, 1))

function _sm_beta(subject, body, toward)
    r, v = _sm_rv(subject, body)
    return π / 2 - _sm_angle(orbit_normal(r, v), toward(body, subject))
end

@testset "Sun and Moon quantities match their pieces from Epicycle's models" begin
    for c in _SM_CASES, subject in _sm_subjects(c)
        @testset "$(c.name), $(nameof(typeof(subject)))" begin
            # About the central body of the subject's own frame, and about each body given explicitly
            @test beta_angle(subject) ≈ _sm_beta(subject, c.body, _sm_sun) atol = 1e-12
            for body in (earth, moon, mars)
                @test beta_angle(subject, body) ≈ _sm_beta(subject, body, _sm_sun) atol = 1e-12
                r, v = _sm_rv(subject, body)
                @test orbit_normal_sun_angle(subject, body) ≈
                      _sm_angle(orbit_normal(r, v), _sm_sun(body, subject) .- r) atol = 1e-12
            end
            for body in (earth, mars)
                @test lunar_beta_angle(subject, body) ≈ _sm_beta(subject, body, _sm_moon) atol = 1e-12
            end

            # The six angles among the Sun, the Earth, the Moon and the spacecraft, about the Earth
            r, _ = _sm_rv(subject, earth)
            s, m = _sm_sun(earth, subject), _sm_moon(earth, subject)
            @test sun_earth_spacecraft_angle(subject)  ≈ _sm_angle(s, r)          atol = 1e-12
            @test sun_spacecraft_earth_angle(subject)  ≈ _sm_angle(s .- r, -r)    atol = 1e-12
            @test sun_spacecraft_moon_angle(subject)   ≈ _sm_angle(s .- r, m .- r) atol = 1e-12
            @test moon_earth_sun_angle(subject)        ≈ _sm_angle(m, s)          atol = 1e-12
            @test moon_earth_spacecraft_angle(subject) ≈ _sm_angle(m, r)          atol = 1e-12
            @test moon_spacecraft_earth_angle(subject) ≈ _sm_angle(m .- r, -r)    atol = 1e-12
        end
    end
end

@testset "the body and the epoch change the answer" begin
    # Guards the test above against an ignored body or epoch argument.
    c = _SM_CASES[1]
    sat = first(_sm_subjects(c))
    @test !(beta_angle(sat, earth) ≈ beta_angle(sat, mars))
    later = Spacecraft(state = CartesianState(c.x), time = c.epoch + 90.0, coord_sys = c.cs)
    @test !(beta_angle(later) ≈ beta_angle(sat))
    @test !(moon_earth_sun_angle(later) ≈ moon_earth_sun_angle(sat))
end

set_frame_theory!(_SM_ORIGINAL_THEORY)

nothing
