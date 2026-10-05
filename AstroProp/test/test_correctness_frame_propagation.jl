# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Propagation respects the spacecraft's coordinate system. The forces integrate in their own frame,
# the central body with ICRF axes; a spacecraft held in any other coordinate system is converted in
# at the start, out wherever a stop condition reads it, and back at the end.
#
# The end-to-end case is a translunar trajectory flown by two propagators in turn, an Earth-centred
# one to the Moon's sphere of influence and a Moon-centred one to lunar periapsis, checked against
# GMAT R2026A (gmat_earth_moon_transfer.script) and flown twice, once with the spacecraft held in
# Earth-centred coordinates and once in Moon-centred ones.

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames, AstroCallbacks
using AstroUniverse
using LinearAlgebra: norm, I
using Test

# ── GMAT truth ───────────────────────────────────────────────────────────────────────────────────
# Point masses Earth, Moon and Sun, DE440 through SPICE, Epicycle's gravitational parameters,
# RungeKutta89 at 1e-13 (unchanged to 0.1 mm from 1e-11 to 1e-14, and with PrinceDormand78).
#
# The initial state is GMAT's EarthICRF report of the script's EarthMJ2000Eq input. GMAT's
# MJ2000Eq-to-ICRF conversion is not the IERS frame bias: at this state it is 0.90 m from SOFA's
# `rb` matrix, which Epicycle's MJ2000Eq matches. GMAT is consistent with itself, since it rotates
# the ephemeris into MJ2000Eq with the same conversion, so the comparison starts both tools from
# the same ICRF state. Epochs are TDB Julian dates minus 2430000; states are Moon-centred ICRF.
const _GMAT_T0_UTC   = "2026-01-01T00:00:00"
const _GMAT_X0_EARTH = [2491.242350649993, -5494.004399970941, -2864.491425983717,
                        10.13590350345997, 3.454065939484507, 2.19038539314474]
const _GMAT_SOI_TDB  = 31042.976737907
const _GMAT_SOI_MOON = [-32908.58164832491, -50059.24794239621, -27933.97482163025,
                        0.7341632940627261, 1.393964538525669, 0.7724730480025555]
const _GMAT_PERI_TDB  = 31043.40009532778
const _GMAT_PERI_MOON = [-3116.19122212095, 2507.018807988788, 1259.603427548602,
                         1.52751725273102, 1.476519159122127, 0.8402442614710658]

const _EARTH_ICRF = CoordinateSystem(earth, ICRF())
const _MOON_ICRF  = CoordinateSystem(moon, ICRF())

_tdb_mjd(sc) = (sc.time.tdb.jd1 - 2430000.0) + sc.time.tdb.jd2

# The same physical state as `x` in Earth ICRF, held in coordinate system `cs`.
function _sc_held_in(x, cs; utc = _GMAT_T0_UTC)
    t0 = Time(utc, UTC(), ISOT())
    sc = Spacecraft(state = CartesianState(x), time = t0, coord_sys = _EARTH_ICRF)
    return Spacecraft(state = CartesianState(to_vector(CartesianState(sc, cs))), time = t0, coord_sys = cs)
end

_integ() = IntegratorConfig(Vern9(); reltol = 1e-13, abstol = 1e-13)

@testset "Earth departure to lunar periapsis, two propagators, vs GMAT" begin
    # The script's gravitational parameters, which are Epicycle's defaults. runtests.jl sets the
    # suite to GMAT's defaults, and those move the events here by 2 ms.
    saved = (earth.mu, moon.mu, sun.mu)
    earth.mu, moon.mu, sun.mu = 398600.4418, 4902.8, 1.32712440018e11
    try
        prop_earth = OrbitPropagator(ForceModel(PointMassGravity(earth, (moon, sun))), _integ())
        prop_moon  = OrbitPropagator(ForceModel(PointMassGravity(moon, (earth, sun))), _integ())

        finals = Dict{String, Vector{Float64}}()
        for (label, cs) in (("held in Earth ICRF", _EARTH_ICRF), ("held in Moon ICRF", _MOON_ICRF))
            @testset "$label" begin
                sc = _sc_held_in(_GMAT_X0_EARTH, cs)

                # Earth-centred to the Moon's sphere of influence
                propagate!(prop_earth, sc, StopAt(position_magnitude, sc, _MOON_ICRF;
                                                  equals = 66100.0, direction = -1))
                x = to_vector(CartesianState(sc, _MOON_ICRF))
                @test sc.coord_sys === cs
                @test abs(norm(x[1:3]) - 66100.0) < 1e-5                         # km
                @test abs(_tdb_mjd(sc) - _GMAT_SOI_TDB) * 86400 < 1e-3           # s; 0.04 ms
                @test norm(x[1:3] - _GMAT_SOI_MOON[1:3]) < 2e-3                  # km; 0.50 m
                @test norm(x[4:6] - _GMAT_SOI_MOON[4:6]) < 1e-7                  # km/s; 0.003 mm/s

                # Moon-centred to periapsis
                propagate!(prop_moon, sc, StopAt(position_dot_velocity, sc, _MOON_ICRF;
                                                 equals = 0.0, direction = 1))
                x = to_vector(CartesianState(sc, _MOON_ICRF))
                @test sc.coord_sys === cs
                @test abs(_tdb_mjd(sc) - _GMAT_PERI_TDB) * 86400 < 1e-3          # s; 0.06 ms
                @test norm(x[1:3] - _GMAT_PERI_MOON[1:3]) < 2e-3                 # km; 0.45 m
                @test norm(x[4:6] - _GMAT_PERI_MOON[4:6]) < 1e-6                 # km/s; 0.07 mm/s
                @test abs(norm(x[1:3]) - norm(_GMAT_PERI_MOON[1:3])) < 1e-3      # km; 0.14 m
                finals[label] = x

                # One history segment per leg, each in the frame its forces integrated in
                segs = sc.history.segments
                @test segs[end-1].coordinate_system.origin === earth
                @test segs[end].coordinate_system.origin === moon
            end
        end
        # The coordinate system the spacecraft is held in does not change the trajectory. What
        # remains is rounding in the conversions and the stop's bisection, 1.2 mm.
        Δ = finals["held in Earth ICRF"] - finals["held in Moon ICRF"]
        @test norm(Δ[1:3]) < 1e-5                                                # km
        @test norm(Δ[4:6]) < 1e-8                                                # km/s
    finally
        earth.mu, moon.mu, sun.mu = saved
    end
end

@testset "a spacecraft held in another frame than the forces'" begin
    x0 = [7000.0, 0.0, 1000.0, 0.0, 7.5, 1.0]
    forces = ForceModel(PointMassGravity(earth, (moon, sun)))
    prop = OrbitPropagator(forces, _integ())

    @testset "the state is converted in and out" begin
        sc_e = _sc_held_in(x0, _EARTH_ICRF)
        sc_m = _sc_held_in(x0, _MOON_ICRF)
        propagate!(prop, sc_e, StopAt(sc_e, PropDurationSeconds(), 3600.0))
        propagate!(prop, sc_m, StopAt(sc_m, PropDurationSeconds(), 3600.0))
        @test sc_m.coord_sys === _MOON_ICRF
        # Before the fix these were 34,130 km apart: the Moon-centred numbers went into the
        # Earth-centred integration as they stood.
        @test norm(to_vector(CartesianState(sc_m, _EARTH_ICRF)) - to_posvel(sc_e)) < 1e-8
        @test sc_m.time.tt == sc_e.time.tt
    end

    @testset "a stop condition sees the epoch it tests" begin
        # The distance from the Moon at a candidate point needs the Moon where it is then. Before
        # the fix the check read the spacecraft's starting epoch, and this stop landed 3,300 km off.
        sc = _sc_held_in([6678.0, 0.0, 0.0, 0.0, 10.9, 0.5], _EARTH_ICRF)
        target = norm(to_vector(CartesianState(sc, _MOON_ICRF))[1:3]) - 5000.0
        propagate!(prop, sc, StopAt(position_magnitude, sc, _MOON_ICRF; equals = target))
        @test abs(norm(to_vector(CartesianState(sc, _MOON_ICRF))[1:3]) - target) < 1e-5
    end

    @testset "an epoch stop is timed from when the propagation starts" begin
        sc = _sc_held_in(x0, _EARTH_ICRF)
        stop = StopAt(sc, Time("2026-01-01T01:00:00", UTC(), ISOT()))
        propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 600.0))    # moves sc after the stop was built
        propagate!(prop, sc, stop)
        @test abs(sc.time.utc - Time("2026-01-01T01:00:00", UTC(), ISOT())) * 86400 < 1e-6
    end

    @testset "OrbitODEProblem" begin
        sc_e = _sc_held_in(x0, _EARTH_ICRF)
        sc_m = _sc_held_in(x0, _MOON_ICRF)
        stm  = STMConfig(Φ = true)
        r_e = solve(OrbitODEProblem(prop, sc_e; duration_s = 3600.0, stm = stm))
        r_m = solve(OrbitODEProblem(prop, sc_m; duration_s = 3600.0, stm = stm))

        # y_final comes back in the spacecraft's coordinate system, at the propagation's end epoch
        tf = Time(sc_e.time.tt.jd1, sc_e.time.tt.jd2 + 3600.0 / 86400.0, TT(), JD())
        back = to_vector(CartesianState(Coordinate(r_m.y_final, _MOON_ICRF, tf), _EARTH_ICRF))
        @test norm(back - r_e.y_final) < 1e-8
        # An origin change leaves the STM as it is
        @test r_m.Φ ≈ r_e.Φ rtol = 1e-10
        # and matches propagate!
        propagate!(prop, sc_e, StopAt(sc_e, PropDurationSeconds(), 3600.0))
        @test norm(to_posvel(sc_e) - r_e.y_final) < 1e-8

        # What cannot be handed back in the spacecraft's frame is refused
        @test_throws ArgumentError solve(OrbitODEProblem(prop, sc_m; duration_s = 60.0, dense = true))
        sc_j2k = _sc_held_in(x0, CoordinateSystem(earth, MJ2000Eq()))
        @test_throws ArgumentError solve(OrbitODEProblem(prop, sc_j2k; duration_s = 60.0, stm = stm))
        @test solve(OrbitODEProblem(prop, sc_j2k; duration_s = 60.0)).y_final isa Vector{Float64}
    end

    @testset "GCRF is ICRF's orientation, and is not converted" begin
        sc = _sc_held_in(x0, CoordinateSystem(earth, GCRF()))
        @test AstroProp._integration_frame(forces, sc) === sc.coord_sys
    end
end
