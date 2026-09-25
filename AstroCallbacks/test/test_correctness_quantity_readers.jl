# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# =============================================================================
# Quantity readers, and `Calc`.
#
# The load-bearing test is the first one: **the new readers agree with the
# shipped `OrbitCalc` stack**. Everything else in the refactor depends on being
# able to switch call sites one at a time, and that is only safe if both paths
# give the same number. When `AbstractCalc` is finally deleted this testset
# goes with it, having done its job.
#
# The rest covers what `Calc` is for — a zero-arg callable that keeps its parts,
# so a history walk can re-apply it to each recorded sample.
# =============================================================================

using Test
using AstroCallbacks
using AstroModels
using AstroStates
using AstroEpochs
using AstroFrames
using AstroUniverse
using LinearAlgebra
using EpicycleBase: label, is_cyclic, cycle, is_settable

const _EPOCH = Time(2458849.5, 0.0, :tdb, :jd)

_subject() = Spacecraft(CartesianState([7000.0, 1000.0, 2000.0, 1.5, 6.5, 1.0]), _EPOCH;
                        coord_sys = CoordinateSystem(earth, ICRF()))

@testset "the new readers agree with the OrbitCalc stack" begin
    # Both must give the same answer for the whole refactor to proceed one call
    # site at a time. Delete this testset with AbstractCalc.
    sat = _subject()
    for (reader, tag) in ((semi_major_axis, SMA()), (eccentricity, Ecc()),
                          (inclination, Inc()),     (raan, RAAN()),
                          (true_anomaly, TA()),     (position_magnitude, PosMag()),
                          (velocity_magnitude, VelMag()))
        @test reader(sat) ≈ get_calc(OrbitCalc(sat, tag)) rtol = 1e-12
    end

    @test position_vector(sat) ≈ get_calc(OrbitCalc(sat, PositionVector())) rtol = 1e-12
    @test position_x(sat)      ≈ get_calc(OrbitCalc(sat, PosX()))           rtol = 1e-12
    @test position_z(sat)      ≈ get_calc(OrbitCalc(sat, PosZ()))           rtol = 1e-12
end

@testset "omitting the frame uses the subject's own" begin
    # Not an inference: a Spacecraft declares its frame through `frame_of`, so
    # the one-argument form reads a stated fact.
    sat = _subject()
    own = AstroFrames.frame_of(sat)
    for reader in (semi_major_axis, eccentricity, inclination, raan,
                   argument_of_periapsis, true_anomaly, position_magnitude)
        @test reader(sat) == reader(sat, own)
    end
end

@testset "the frame changes the answer, and it is meant to" begin
    sat = _subject()
    eq  = CoordinateSystem(earth, MJ2000Eq())
    ec  = CoordinateSystem(earth, MJ2000Ec())

    # An angle measured in a differently-tilted plane is a different angle.
    @test !isapprox(inclination(sat, eq), inclination(sat, ec); atol = 1e-6)
    @test !isapprox(raan(sat, eq),        raan(sat, ec);        atol = 1e-6)

    # But a magnitude about the same origin is not.
    @test position_magnitude(sat, eq) ≈ position_magnitude(sat, ec) rtol = 1e-12
end

@testset "semi-major axis read in a rotating frame is meaningless" begin
    # The hazard FR-CALL-7 accepts, documented so it is not discovered in
    # anger. Earth-fixed velocity is short by ω × r — about 0.4 km/s at LEO —
    # so the energy, and with it the semi-major axis, moves by hundreds of km.
    sat  = _subject()
    itrf = CoordinateSystem(earth, ITRF())

    inertial = semi_major_axis(sat, CoordinateSystem(earth, MJ2000Eq()))
    fixed    = semi_major_axis(sat, itrf)

    @test abs(fixed - inertial) > 100.0          # not a rounding difference
    @test velocity_magnitude(sat, itrf) < velocity_magnitude(sat)
end

@testset "element reaches any representation" begin
    # No registry and no declared state tag: `element` takes the
    # representation as an argument, and the quantity names it in its body.
    sat = _subject()
    eq  = CoordinateSystem(earth, MJ2000Eq())

    @test element(KeplerianState, sat, eq, :sma) ≈ semi_major_axis(sat, eq) rtol = 1e-12
    @test isfinite(element(OutGoingAsymptoteState, sat, eq, :c3))
    @test isfinite(element(ModifiedEquinoctialState, sat, eq, :p))
end

@testset "Calc evaluates to the same number as the direct call" begin
    sat = _subject()
    ec  = CoordinateSystem(earth, MJ2000Ec())

    @test Calc(raan, sat, ec)()          == raan(sat, ec)
    @test Calc(semi_major_axis, sat)()   == semi_major_axis(sat)
    @test Calc(epoch, sat)()             == AstroFrames.epoch_of(sat)
end

@testset "Calc keeps its parts, so it can be re-applied" begin
    # This is what a history walk does with each recorded sample, and what a
    # closure cannot do — its subject is sealed in.
    ec   = CoordinateSystem(earth, MJ2000Ec())
    sat  = _subject()
    other = Spacecraft(CartesianState([8000.0, 0.0, 0.0, 0.0, 7.0, 1.0]), _EPOCH;
                       coord_sys = CoordinateSystem(earth, ICRF()))

    q = Calc(raan, sat, ec)
    @test reapply(q, other) == raan(other, ec)
    @test reapply(q, other) != q()               # it really moved

    # The dependency travels with it; only the subject is replaced.
    @test reapply(Calc(inclination, sat, ec), other) == inclination(other, ec)
end

@testset "traits pass through a Calc" begin
    sat = _subject()
    ec  = CoordinateSystem(earth, MJ2000Ec())

    @test label(Calc(raan, sat, ec))     == label(raan)
    @test is_cyclic(Calc(raan, sat, ec)) == true
    @test cycle(Calc(raan, sat, ec))     ≈ 2π
    @test is_cyclic(Calc(semi_major_axis, sat)) == false

    # It names itself, which is what a report needs.
    @test occursin("RAAN", sprint(show, Calc(raan, sat, ec)))
end

@testset "shipped quantities carry labels" begin
    for f in (semi_major_axis, eccentricity, inclination, raan,
              argument_of_periapsis, true_anomaly, position_vector,
              velocity_vector, position_x, position_magnitude, epoch)
        @test label(f) != "quantity"          # not the default
        @test !isempty(label(f))
    end
end

@testset "position_dot_velocity" begin
    sat = _subject()
    eq  = CoordinateSystem(earth, MJ2000Eq())

    # Agrees with the shipped tag, and with the dot product it claims to be.
    @test position_dot_velocity(sat) ≈ get_calc(OrbitCalc(sat, PosDotVel())) rtol = 1e-12
    @test position_dot_velocity(sat) ≈ dot(position_vector(sat), velocity_vector(sat)) rtol = 1e-12
    @test position_dot_velocity(sat, eq) ≈ dot(position_vector(sat, eq),
                                               velocity_vector(sat, eq)) rtol = 1e-12

    # It does not wrap, which is the reason it exists: it is the stopping
    # condition for periapsis and apoapsis now that angles are off the table.
    @test is_cyclic(position_dot_velocity) == false
    @test label(position_dot_velocity) != "quantity"

    # Zero exactly at an apsis, and it changes sign across one. A circular
    # orbit is r ⟂ v everywhere, so use an eccentric one.
    peri = Spacecraft(KeplerianState(8000.0, 0.2, 0.5, 0.3, 0.4, 0.0), _EPOCH;
                      coord_sys = CoordinateSystem(earth, MJ2000Eq()))
    apo  = Spacecraft(KeplerianState(8000.0, 0.2, 0.5, 0.3, 0.4, float(π)), _EPOCH;
                      coord_sys = CoordinateSystem(earth, MJ2000Eq()))
    @test abs(position_dot_velocity(peri)) < 1e-9
    @test abs(position_dot_velocity(apo))  < 1e-9

    # Rising through zero at periapsis, falling at apoapsis — which is what
    # `direction = 1` and `-1` pick out.
    after_peri = Spacecraft(KeplerianState(8000.0, 0.2, 0.5, 0.3, 0.4, 0.1), _EPOCH;
                            coord_sys = CoordinateSystem(earth, MJ2000Eq()))
    after_apo  = Spacecraft(KeplerianState(8000.0, 0.2, 0.5, 0.3, 0.4, π + 0.1), _EPOCH;
                            coord_sys = CoordinateSystem(earth, MJ2000Eq()))
    @test position_dot_velocity(after_peri) > 0
    @test position_dot_velocity(after_apo)  < 0
end
