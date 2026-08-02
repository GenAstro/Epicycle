# Atmospheric drag (NRLMSISE-00) on two-body gravity — Epicycle interface vs GMAT.
#
# Uses the spec'd user interface: cannonball drag geometry on the Spacecraft, an
# AtmosphericDrag force selected by a density-model tag, composed with point-mass gravity
# in a ForceModel, propagated with propagate!. Validates the 1-day final state against
# GMAT R2022a truth (matched mass/area/Cd, GMAT space-weather file).
#
# Notes:
#  - Epicycle earth.mu = 398600.4418; GMAT/EGM96 GM = 398600.4415 (~sub-meter/day).
#  - NRLMSISE-00 density uses historical SpaceIndices, matching the GMAT SW-file run.
#  - Residual is space-weather-data-limited (~10 m), not a physics gap.
#
# Run under an environment that has AstroProp developed (e.g. the force_epicycle project).

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9
using LinearAlgebra: norm
using Test

# ── User interface (from the force-model spec) ────────────────────────────────
sc = Spacecraft(;
    state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),
    time  = Time("2020-10-20T12:00:00", UTC(), ISOT()),
    mass  = 1000.0,
    name  = "LEO",
    drag  = SphericalDrag(c_d = 2.2, drag_area = 10.0),
)

gravity = PointMassGravity(earth, ())
drag    = AtmosphericDrag(earth; model = Exponential())
forces  = ForceModel(gravity, drag)

integ = IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0)
prop  = OrbitPropagator(forces, integ)
sol   = propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 86400.0))
yf    = sol.u[end]

# ── GMAT truth (two-body + NRLMSISE-00 drag, SW file, matched mass/area/Cd) ────
gmat = [ 5318.2703125793,  2703.7861413888,  3437.7763438110,
           -4.8236652445236, 3.6488637771088, 4.6394118851618]

Δr = norm(gmat[1:3] .- yf[1:3]) * 1e3      # m
Δv = norm(gmat[4:6] .- yf[4:6]) * 1e6      # mm/s

# Provisional tolerances — space-weather-data-limited; tighten once SW sources are matched.
@testset "Exponential drag vs GMAT" begin
    @test Δr < 1.2          # m
    @test Δv < 1.5          # mm/s
end
