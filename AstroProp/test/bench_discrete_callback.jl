# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# DIAGNOSTIC — not a production callback path.
#
# Purpose: measure how much of the callback bench's cost is DifferentialEquations'
# ContinuousCallback machinery (Vern9 dense-output evaluation + root-finding) versus
# our own calc cost. To isolate this, we monkey-patch `AstroProp._build_callback` to
# return a `DiscreteCallback` instead of a `ContinuousCallback`.
#
# What we lose: root-finding precision. DiscreteCallback fires at the accepted-step
# endpoint where MeanSMA has already crossed below the target — overshoot up to one
# step-length of orbit time (~5–20 s at LEO with tol=1e-10). Fine for a diagnostic.
#
# What we gain: no dense-output interpolant evaluations, no sub-step polling. Should
# be ~1 calc call per accepted step (vs 3 with the fixed ContinuousCallback).
#
# The monkey-patch persists for the REPL session — restart Julia between this and
# bench_prop_callback.jl so results aren't cross-contaminated.
#
# Run manually:
#   using Pkg; Pkg.activate("c:/Users/steve/Dev/epicycle-dev")
#   include(raw"c:\Users\steve\Dev\Epicycle\AstroProp\test\bench_discrete_callback.jl")

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroCallbacks: MeanSMA, OrbitCalc, get_calc
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9, DiscreteCallback, terminate!
using BenchmarkTools

# ── Monkey-patch: DiscreteCallback in place of ContinuousCallback ────────────
# Same signature as the original; body checks the direction and returns Bool.
function AstroProp._build_callback(cond::AstroProp.StopAt, dynsys)
    subject = cond.subject
    var     = cond.var
    target  = cond.target
    dir     = cond.direction
    calc    = AstroProp.make_calc(subject, var)

    function cond_fn(u, t, integ)
        AstroProp._subject_update_from_u!(subject, dynsys, u)
        val = get_calc(calc)
        return dir < 0 ? val <= target :
               dir > 0 ? val >= target :
                         false                     # dir==0 unsupported by this diagnostic
    end
    term!(integ) = terminate!(integ)

    return DiscreteCallback(cond_fn, term!; save_positions=(false, false))
end

const MEAN_SMA_TARGET = 6708.0    # km

function make_sat()
    kep = KeplerianState(6778.137, 0.0, deg2rad(51.6), 0.0, 0.0, 0.0)
    Spacecraft(;
        state = CartesianState(kep, earth.mu),
        time  = Time("2020-01-01T00:00:00", UTC(), ISOT()),
        name  = "SKSat",
        mass  = 1000.0,
        drag  = SphericalDrag(c_d = 2.2, drag_area = 10.0),
        save_history = false,
    )
end

function make_prop()
    gravity = PointMassGravity(earth, ())
    drag    = AtmosphericDrag(earth; model = Exponential())
    forces  = ForceModel(gravity, drag)
    integ   = IntegratorConfig(Vern9(); reltol = 1e-10, abstol = 1e-10, dt = 60.0)
    return OrbitPropagator(forces, integ)
end

function run_once(sat, prop)
    return propagate!(prop, sat,
        StopAt(sat, MeanSMA(), MEAN_SMA_TARGET; direction=-1))
end

# ── Warmup + reference ───────────────────────────────────────────────────────
println("Warming up (DiscreteCallback variant)...")
prop = make_prop(); sat = make_sat()
sol  = run_once(sat, prop)
nsteps_ref  = length(sol.t)
tf_ref_days = sol.t[end] / 86400.0
println("  ODE steps:  ", nsteps_ref)
println("  Final time: ", round(tf_ref_days; digits=3), " days")

# ── Benchmark ────────────────────────────────────────────────────────────────
println("\nBenchmarking propagate! (DiscreteCallback) to MeanSMA = $(MEAN_SMA_TARGET) km ...")
b = @benchmark run_once(sat, prop) setup=(sat = make_sat(); prop = make_prop()) evals=1 samples=5 seconds=600

show(stdout, MIME"text/plain"(), b)
println()

t_median_s = median(b.times) / 1e9
mem_median = median(b.memory)
allocs_med = median(b.allocs)
println("\nSUMMARY: median = $(round(t_median_s; digits=3)) s, " *
        "memory = $(round(mem_median / 2^20; digits=2)) MiB, " *
        "allocs = $(allocs_med), ODE steps = $(nsteps_ref)")

println("\nBENCH_DISCRETE_CALLBACK_DONE")
