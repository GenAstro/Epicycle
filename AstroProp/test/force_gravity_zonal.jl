# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

# Zonal (J2–J5) gravity — open, native. Standalone validation script; kept dependency-light
# (Base @assert + println only) so it runs directly in the AstroProp environment, e.g.
#   julia --project=<...>/Epicycle/AstroProp  AstroProp/test/force_gravity_zonal.jl
#
# Two checks:
#   1. J2 self-check (no GMAT): native Zonal acceleration at degree 2 vs the closed-form J2.
#   2. GMAT J2–J5: 1-day LEO final state vs a GMAT run with EGM96, degree = 5, order = 0.
#      Fill `gmat_j2j5` with the truth vector and the assertions activate.

using AstroProp
using AstroProp: geopotential_data, geopotential_accel     # seam internals for the self-check
using Test
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth, frame_theory, set_frame_theory!, FK5
using OrdinaryDiffEqVerner: Vern9
using LinearAlgebra: norm


# ── 1. Native J2 self-check — closed form, no GMAT ────────────────────────────
function _j2_closed_form(r, μ, Re, J2)
    rn = norm(r); s = r[3] / rn
    k  = -1.5 * J2 * μ * (Re / rn)^2 / rn^3
    return [k * (1 - 5s^2) * r[1], k * (1 - 5s^2) * r[2], k * (3 - 5s^2) * r[3]]
end

let
    z    = Zonal()
    data = geopotential_data(z, earth, 2, 0)
    r    = [6878.137e3, 1200.0e3, 2200.0e3]                      # arbitrary Earth-fixed point [m]
    a    = geopotential_accel(z, data, r, 0.0, 2, 0)            # central + J2
    aref = (-data.mu_m .* r ./ norm(r)^3) .+
           _j2_closed_form(r, data.mu_m, data.Re_m, 1.0826266835531513e-3)
    relerr = norm(a .- aref) / norm(aref)
    relerr < 1e-13 || @warn "Zonal J2 acceleration does not match the closed form" relerr
    @test relerr < 1e-13

end

# ── 2. Full propagation vs GMAT (J2–J5) ───────────────────────────────────────
sc = Spacecraft(;
    state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
    time  = Time("2020-10-20T12:00:00", UTC(), ISOT()),
    mass  = 1000.0,
    name  = "LEO",
)

# In Epicycle's FK5 Earth axes. GMAT's leave out the IERS celestial-pole offsets, so its pole sits
# about 47 mas from the IERS pole in 2020; for a zonal field only the pole matters, and that moves
# the final state by the difference this test allows.
original_theory = frame_theory()
set_frame_theory!(FK5())
gravity = try
    HarmonicGravity(earth; degree = 5, order = 0, model = Zonal())
finally
    set_frame_theory!(original_theory)
end
forces  = ForceModel(gravity)

integ = IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0)
prop  = OrbitPropagator(forces, integ)
sol   = propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 86400.0))
yf    = sol.u[end]


# GMAT truth — EGM96 coefficients, zonal J2–J5 (Degree = 5, Order = 0), same epoch/state/duration.
gmat_j2j5 = [ 4634.0492943971,  2919.8542494901,  4172.5997431458,
                -5.6004330855833, 3.4747294213693, 3.8223866712728]

Δr = norm(gmat_j2j5[1:3] .- yf[1:3]) * 1e3      # m
Δv = norm(gmat_j2j5[4:6] .- yf[4:6]) * 1e6      # mm/s

# Measured 11.8 cm and 0.11 mm/s on 2026-10-04 in Epicycle's FK5 axes. In GMAT's own axes this test
# held 5 cm, so most of the difference is GMAT's pole.
@testset "Zonal gravity vs GMAT" begin
    @test Δr < 0.15
    @test Δv < 0.15
end

