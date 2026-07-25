# Spherical-harmonic gravity (EGM96 70x70) — Epicycle interface vs GMAT.
#
# Uses the spec'd user interface: a Spacecraft, a HarmonicGravity force selected by a
# model tag, composed in a ForceModel, propagated with propagate!. Validates the 1-day
# final state against GMAT R2022a truth.
#
# Run under an environment that has AstroProp developed (e.g. the force_epicycle project).

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9
using LinearAlgebra: norm
using Printf
using Test

# ── User interface (from the force-model spec) ────────────────────────────────
sc = Spacecraft(;
    state = CartesianState([6878.137, 0.0, 0.0, 0.0, 4.71754, 5.99820]),  # J2000, km/km·s⁻¹
    time  = Time("2020-10-20T12:00:00", UTC(), ISOT()),
    mass  = 1000.0,
    name  = "LEO",
)

gravity = HarmonicGravity(earth; degree = 70, order = 70, model = EGM96())
forces  = ForceModel(gravity)

integ = IntegratorConfig(Vern9(); reltol = 1e-12, abstol = 1e-12, dt = 60.0)
prop  = OrbitPropagator(forces, integ)
sol   = propagate!(prop, sc, StopAt(sc, PropDurationSeconds(), 86400.0))
yf    = sol.u[end]

# ── GMAT truth (EGM96 70x70, 1 day, UTC epoch) ────────────────────────────────
gmat = [ 4622.2958511411,  2927.2178362928,  4180.4514680258,
           -5.6123246712046, 3.4671360280002, 3.8116856436408]

Δr = norm(gmat[1:3] .- yf[1:3]) * 1e3      # m
Δv = norm(gmat[4:6] .- yf[4:6]) * 1e6      # mm/s

println("\n===== SH GRAVITY (EGM96 70x70) — Epicycle vs GMAT =====")
for (i, lab) in enumerate(("x","y","z","vx","vy","vz"))
    @printf("%-3s  epi = % .10f   gmat = % .10f   Δ = % .3e\n", lab, yf[i], gmat[i], gmat[i]-yf[i])
end
@printf("|Δr| = %.4f m    |Δv| = %.4f mm/s\n", Δr, Δv)

# Provisional tolerances — catch gross regressions; tighten once locked.
@testset "SH gravity vs GMAT" begin
    @test Δr < 0.05          # m
    @test Δv < 0.05          # mm/s
end
