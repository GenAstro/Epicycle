# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Solar radiation pressure about the Sun. The Sun cannot occult its own light, so the spacecraft is
# in full sunlight everywhere and there are no shadow edges to stop at. Before this was handled, the
# Sun was taken as its own occulting body and the acceleration was NaN.

using AstroProp
using AstroModels, AstroStates, AstroEpochs, AstroFrames
using AstroUniverse
using LinearAlgebra: norm
using Test

@testset "SRP about the Sun" begin
    t  = Time("2025-01-01T00:00:00", UTC(), ISOT())
    x  = [1.5e8, 2.0e7, -1.0e7, 0.0, 30.0, 1.0]
    sc = Spacecraft(state = CartesianState(x), time = t, coord_sys = CoordinateSystem(sun, ICRF()),
                    mass = 1000.0, srp = SphericalSRP(c_r = 1.8, srp_area = 10.0))
    for shadow in (DualCone(), SmoothedConical())
        f = SolarRadiationPressure(sun; shadow = shadow)
        a = AstroProp.accel_eval!(f, t, x, zeros(6), sc, nothing)[4:6]
        r = x[1:3]
        # Full sunlight: Cr·A/m · Φ/c · (AU/r)², directed away from the Sun.
        expected = 1.8 * 10.0 / 1000.0 * 1367.0 / 2.99792458e8 *
                   (149597870.691 / norm(r))^2 / 1.0e3 .* r ./ norm(r)
        @test all(isfinite, a)
        @test a ≈ expected rtol = 1e-14
        @test AstroProp._n_kinks(f) == 0
    end

    # A propagation about the Sun with SRP runs, with no shadow stops.
    forces = ForceModel(PointMassGravity(sun, (earth, jupiter)), SolarRadiationPressure(sun))
    prop   = OrbitPropagator(forces, IntegratorConfig(Vern9(); reltol = 1e-10, abstol = 1e-10,
                                                      dt = 3600.0))
    propagate!(prop, sc, StopAt(sc, PropDurationDays(), 30.0))
    @test all(isfinite, sc.state.state)
end
