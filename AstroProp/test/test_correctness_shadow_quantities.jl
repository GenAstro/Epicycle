# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# earth_shadow_fraction and moon_shadow_fraction (#11, section 4).
#
# Truth:
# - Geometry with a known answer: a spacecraft behind the Earth or the Moon in line with the Sun
#   is in umbra, one on the sunward side in full sun, and the fraction rises continuously from 0
#   to 1 across the shadow edge.
# - Self-consistency with the force: at points across the penumbra the quantity, given the same
#   eclipse model, equals one minus the ratio of the SRP acceleration SolarRadiationPressure
#   applies there to the same acceleration in full sun, for DualCone and SmoothedConical, so a stopping condition on the quantity and the force see the same
#   shadow. Subjects are Spacecraft held in a frame other than the occulter's ICRF, and
#   Coordinates, so the frame conversion is exercised.
# =============================================================================

using Test
using LinearAlgebra: norm, normalize, cross
using AstroProp
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using AstroCallbacks: quantities, subject_type, Calc
using EpicycleBase: label

const _SQ_EPOCH = Time("2024-03-20T12:00:00", UTC(), ISOT())
const _SQ_JD    = _SQ_EPOCH.tdb.jd

_sq_coord(r, cs) = Coordinate([r; 0.0; 0.0; 0.0], cs, _SQ_EPOCH)

@testset "shadow fraction: umbra, full sun and a continuous edge" begin
    ŝ = normalize(translate(earth, sun, _SQ_JD))
    @test earth_shadow_fraction(_sq_coord(-7000.0 .* ŝ, EarthICRF)) == 1
    @test earth_shadow_fraction(_sq_coord( 7000.0 .* ŝ, EarthICRF)) == 0

    # Crossing the shadow edge sideways, the hidden fraction rises continuously from 0 to 1.
    side = normalize(cross(ŝ, [0.0, 0.0, 1.0]))
    fracs = [earth_shadow_fraction(_sq_coord(-7000.0 .* ŝ .+ y .* side, EarthICRF))
             for y in range(6300.0, 6460.0; length = 81)]
    @test fracs[1] == 1 && fracs[end] == 0
    @test issorted(fracs; rev = true)
    @test any(f -> 0 < f < 1, fracs)

    # 2000 km behind the Moon, in line with the Sun: inside its umbra, which is 374,000 km long.
    m  = translate(earth, moon, _SQ_JD)
    ŝm = normalize(translate(moon, sun, _SQ_JD))
    @test moon_shadow_fraction(_sq_coord(m .- 2000.0 .* ŝm, EarthICRF)) == 1
    @test moon_shadow_fraction(_sq_coord([7000.0, 300.0, 1200.0], EarthICRF)) == 0
end

# One minus the SRP acceleration the force applies at `x` (about `body`, ICRF axes), over the
# same acceleration with no shadow.
function _sq_from_force(body, x, model)
    sc = Spacecraft(state = CartesianState(x), time = _SQ_EPOCH, mass = 1000.0,
                    coord_sys = CoordinateSystem(body, ICRF()),
                    srp = SphericalSRP(c_r = 1.5, srp_area = 6.0))
    srp = SolarRadiationPressure(body; shadow = model)
    dy = zeros(6)
    AstroProp.accel_eval!(srp, _SQ_EPOCH, x, dy, sc, [])
    d = norm(x[1:3] .- translate(body, sun, _SQ_JD))
    full = 1.5 * 6.0 / 1000.0 * srp.solar_flux / 299792458.0 * (srp.nominal_sun / d)^2 / 1.0e3
    return 1 - norm(dy[4:6]) / full
end

@testset "shadow fraction is the shadow SolarRadiationPressure applies, for each eclipse model" begin
    for (body, quantity, R) in ((earth, earth_shadow_fraction, 6378.137),
                                (moon,  moon_shadow_fraction,  1738.0)),
        model in (DualCone(), SmoothedConical())
        ŝ = normalize(translate(body, sun, _SQ_JD))
        side = normalize(cross(ŝ, [0.0, 0.0, 1.0]))
        held_in = CoordinateSystem(body, MJ2000Eq())          # not ICRF, so it is converted
        for y in range(0.9R, 1.1R; length = 21)               # across the umbra, penumbra and sun
            r = -(R + 1000.0) .* ŝ .+ y .* side
            x = [r; 0.0; 0.0; 0.0]
            expected = _sq_from_force(body, x, model)
            r_held = to_vector(CartesianState(Coordinate(x, CoordinateSystem(body, ICRF()), _SQ_EPOCH),
                                              held_in))
            sat = Spacecraft(state = CartesianState(r_held), time = _SQ_EPOCH, coord_sys = held_in)
            @test quantity(sat, model) ≈ expected atol = 1e-12
            @test quantity(Coordinate(x, CoordinateSystem(body, ICRF()), _SQ_EPOCH), model) ≈ expected atol = 1e-12
            model isa DualCone && @test quantity(sat) == quantity(sat, model)      # the default
        end
    end
end

@testset "shadow fractions are quantities" begin
    # The eclipse model is a positional argument, so a Calc, and with it StopAt and history, carry it
    sat = Spacecraft(state = CartesianState([-7000.0, 0.0, 0.0, 0.0, -7.5, 0.0]), time = _SQ_EPOCH)
    @test Calc(earth_shadow_fraction, sat, SmoothedConical())() == earth_shadow_fraction(sat, SmoothedConical())
    @test label(earth_shadow_fraction) == "Earth shadow fraction"
    @test label(moon_shadow_fraction) == "Moon shadow fraction"
    @test subject_type(earth_shadow_fraction) == Union{Spacecraft, Coordinate}
    # quantities looks in one module; these are defined in AstroProp
    @test earth_shadow_fraction in quantities(Spacecraft, AstroProp)
    @test moon_shadow_fraction in quantities(Coordinate, AstroProp)
end

nothing
