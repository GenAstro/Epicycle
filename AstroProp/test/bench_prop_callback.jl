# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Benchmark: same LEO + Exponential-drag prop as bench_prop_exponential.jl, but propagated
# until the Brouwer mean SMA decays to 6708 km. The MeanSMA callback is evaluated at every
# integration step and fires exactly once (at the stop condition).
#
# Not 1:1 with GMAT: this runs to a physical condition instead of a fixed 30-day window.
# The intent is to measure the cost of the mean-SMA calc in a realistic station-keeping
# workload, where the trigger event determines the propagation length.
#
# Run manually (not wired into runtests.jl):
#   using Pkg; Pkg.activate("c:/Users/steve/Dev/epicycle-dev")
#   include(joinpath(pkgdir(AstroProp), "test", "bench_prop_callback.jl"))

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroCallbacks: MeanSMA
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9
using BenchmarkTools

const MEAN_SMA_TARGET = 6708.0    # km, decay trigger

# Build a fresh spacecraft for each sample — propagate! mutates sc.state/sc.time.
# NOTE: state is CartesianState (not KeplerianState) so `set_posvel!` inside the
# callback path hits the fast path — avoids a redundant Cart→Kep conversion per
# callback invocation. Orbit is unchanged.
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

# ── Warmup + a reference run for step count and end conditions ───────────────
println("Warming up...")
prop = make_prop()
sat  = make_sat()
sol  = run_once(sat, prop)
nsteps_ref  = length(sol.t)
tf_ref_days = sol.t[end] / 86400.0
yf_ref      = sol.u[end]
println("  ODE steps:     ", nsteps_ref)
println("  Final time:    ", round(tf_ref_days; digits=3), " days")
println("  Final state:   ", yf_ref)

# ── Benchmark ────────────────────────────────────────────────────────────────
println("\nBenchmarking propagate! to MeanSMA = $(MEAN_SMA_TARGET) km ...")
b = @benchmark run_once(sat, prop) setup=(sat = make_sat(); prop = make_prop()) evals=1 samples=5 seconds=600

show(stdout, MIME"text/plain"(), b)
println()

t_median_s = median(b.times) / 1e9
mem_median = median(b.memory)
allocs_med = median(b.allocs)
println("\nSUMMARY: median = $(round(t_median_s; digits=3)) s, " *
        "memory = $(round(mem_median / 2^20; digits=2)) MiB, " *
        "allocs = $(allocs_med), ODE steps = $(nsteps_ref)")

println("\nBENCH_PROP_CALLBACK_DONE")
