# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Benchmark: pure propagation with two-body gravity + Exponential atmospheric drag.
#
# Config mirrors examples/Ex_StationKeeping.jl (LEO, i = 51.6°, SMA = 6778.137 km,
# Cd = 2.2, A = 10 m², m = 1000 kg, Vern9, tol = 1e-10) but strips the solver/calcs
# infrastructure — just one call to `propagate!` for 30 days.
#
# Run manually (not wired into runtests.jl):
#   using Pkg; Pkg.activate("AstroProp"); using TestEnv
#   TestEnv.activate("AstroProp") do
#       include(joinpath(pkgdir(AstroProp), "test", "bench_prop_exponential.jl"))
#   end

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9
using BenchmarkTools

const PROP_DURATION_DAYS = 30.0
const PROP_DURATION_SEC  = PROP_DURATION_DAYS * 86400.0

# Build a fresh spacecraft for each sample — propagate! mutates sc.state/sc.time.
function make_sat()
    Spacecraft(;
        state = KeplerianState(6778.137, 0.0, deg2rad(51.6), 0.0, 0.0, 0.0),
        time  = Time("2020-01-01T00:00:00", UTC(), ISOT()),
        name  = "SKSat",
        mass  = 1000.0,
        drag  = SphericalDrag(c_d = 2.2, drag_area = 10.0),
        save_history = false,   # benchmark pure integration, not ephemeris save
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
    return propagate!(prop, sat, StopAt(sat, PropDurationSeconds(), PROP_DURATION_SEC))
end

# ── Warmup + a single reference run for step count / final state ─────────────
println("Warming up...")
prop = make_prop()
sat  = make_sat()
sol  = run_once(sat, prop)
nsteps_ref = length(sol.t)
yf_ref     = sol.u[end]
println("  ODE steps:     ", nsteps_ref)
println("  Final state:   ", yf_ref)

# ── Benchmark ────────────────────────────────────────────────────────────────
println("\nBenchmarking propagate! for $(PROP_DURATION_DAYS) days ...")
# `setup` rebuilds sat and prop on every sample so drag/EOP construction cost
# doesn't leak into the timing, and each sample starts from the same t0.
b = @benchmark run_once(sat, prop) setup=(sat = make_sat(); prop = make_prop()) evals=1 samples=5 seconds=120

show(stdout, MIME"text/plain"(), b)
println()

# ── Summary line for quick eyeballing across runs ────────────────────────────
t_median_s   = median(b.times) / 1e9
mem_median   = median(b.memory)
allocs_med   = median(b.allocs)
println("\nSUMMARY: median = $(round(t_median_s; digits=3)) s, " *
        "memory = $(round(mem_median / 2^20; digits=2)) MiB, " *
        "allocs = $(allocs_med), ODE steps = $(nsteps_ref)")

println("\nBENCH_PROP_EXPONENTIAL_DONE")
