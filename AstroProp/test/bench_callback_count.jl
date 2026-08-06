# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# Instrumentation: count how many times MeanSMA is evaluated per accepted ODE step
# during a propagate! with a MeanSMA StopAt callback.
#
# Method: monkey-patch AstroCallbacks._evaluate(::MeanSMA, ::BrouwerMeanLongState) to
# increment a Ref. This is the same one-liner AstroCallbacks defines in its own source;
# we just interpose a counter. Runs one full propagation (no benchmarking wrapper) and
# prints:
#   - wall time
#   - final time (days)

#   - accepted ODE steps
#   - total MeanSMA evaluations
#   - calls per accepted step   <-- the number we care about
#
# Hypothesis: this ratio is ~28, based on the earlier ratio of alloc counts.
# If confirmed, the fix is in the callback wiring (ContinuousCallback config or the
# root-finding polling schedule), not in the MeanSMA calc itself.
#
# Run manually:
#   using Pkg; Pkg.activate("c:/Users/steve/Dev/epicycle-dev")
#   include(raw"c:\Users\steve\Dev\Epicycle\AstroProp\test\bench_callback_count.jl")

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroCallbacks
using AstroCallbacks: MeanSMA
using AstroStates: BrouwerMeanLongState
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9

const MEAN_SMA_TARGET = 6768.0    # km
const _MEANSMA_CALLS  = Ref(0)

# Monkey-patch: same body as the original method, plus a counter bump.
# Original (AstroCallbacks/src/orbitcalc_meansma.jl):
#     _evaluate(::MeanSMA, bml::BrouwerMeanLongState) = bml.sma
function AstroCallbacks._evaluate(::MeanSMA, bml::BrouwerMeanLongState)
    _MEANSMA_CALLS[] += 1
    return bml.sma
end


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

# ── Warmup (compile everything) ──────────────────────────────────────────────
println("Warmup pass (short prop, discards counter)...")
sat = make_sat(); prop = make_prop()
propagate!(prop, sat,
    StopAt(sat, MeanSMA(), MEAN_SMA_TARGET; direction=-1),
    StopAt(sat, PropDurationSeconds(), 3600.0))     # 1 hour, just to compile

# ── Instrumented run ─────────────────────────────────────────────────────────
println("\nInstrumented propagation to MeanSMA = $(MEAN_SMA_TARGET) km ...")
sat = make_sat(); prop = make_prop()
_MEANSMA_CALLS[] = 0
t0  = time()
sol = propagate!(prop, sat, StopAt(sat, MeanSMA(), MEAN_SMA_TARGET; direction=-1))
elapsed_s = time() - t0

nsteps          = length(sol.t)
ncalls          = _MEANSMA_CALLS[]
calls_per_step  = ncalls / nsteps
tf_days         = sol.t[end] / 86400.0
us_per_call     = 1e6 * elapsed_s / ncalls
us_per_step     = 1e6 * elapsed_s / nsteps

println()
println("═══════════════════════════════════════════════════════════════")
println("  wall time                = ", round(elapsed_s;      digits=2), " s")
println("  final time               = ", round(tf_days;        digits=3), " days")
println("  ODE steps (accepted)     = ", nsteps)
println("  MeanSMA evaluations      = ", ncalls)
println("  ── calls per accepted step = ", round(calls_per_step; digits=2))
println("  wall time per step       = ", round(us_per_step;   digits=1), " μs")
println("  wall time per calc call  = ", round(us_per_call;   digits=1), " μs")
println("═══════════════════════════════════════════════════════════════")

println("\nBENCH_CALLBACK_COUNT_DONE")
