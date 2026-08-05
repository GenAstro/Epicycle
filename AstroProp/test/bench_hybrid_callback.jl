# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LGPL-3.0-only OR LicenseRef-GenAstro-Commercial OR LicenseRef-GenAstro-Evaluation
#
# DIAGNOSTIC — hybrid DiscreteCallback + one-shot interpolant root-find.
#
# Design:
#   - Per-step polling uses DiscreteCallback → one cheap g(integrator.u) eval per accepted
#     step, no interpolant construction.
#   - When a sign change is detected between consecutive step endpoints, `affect!` runs a
#     bisection on the Vern9 dense-output interpolant to nail down t_root and u_root, then
#     terminates. This is the ONE step where we pay for the interpolant, not every step.
#
# Expected result: same wall-time as bench_discrete_callback.jl (~11 s at full scale) but
# with the exact-root precision of ContinuousCallback (overshoot < bisection tol).
#
# What's exposed for verification:
#   _ROOT_T[]      — interpolated root time (elapsed sec from integration t0)
#   _ROOT_U[]      — interpolated 6-elt state at root
#   _ROOT_ITERS[]  — bisection iterations to converge (should be ~30–40)
#   sol.t[end]     — accepted-step endpoint where affect! actually fired (overshoot side)
#
# The monkey-patch persists for the REPL session — restart Julia between this and the
# other bench_*.jl files.
#
# Run manually:
#   using Pkg; Pkg.activate("c:/Users/steve/Dev/epicycle-dev")
#   include(raw"c:\Users\steve\Dev\Epicycle\AstroProp\test\bench_hybrid_callback.jl")

using AstroProp
using AstroModels, AstroStates, AstroEpochs
using AstroCallbacks: MeanSMA, OrbitCalc, get_calc
using AstroUniverse: earth
using OrdinaryDiffEq: Vern9, DiscreteCallback, terminate!
using BenchmarkTools

# Populated by affect! on each run; reset in run_once.
const _ROOT_T     = Ref(NaN)
const _ROOT_U     = Ref{Vector{Float64}}(zeros(6))
const _ROOT_ITERS = Ref(0)

# ── Monkey-patch: DiscreteCallback + on-crossing root-find ───────────────────
function AstroProp._build_callback(cond::AstroProp.StopAt, dynsys)
    subject = cond.subject
    var     = cond.var
    target  = cond.target
    dir     = cond.direction
    calc    = AstroProp.make_calc(subject, var)

    g_prev = Ref(NaN)     # value of (calc - target) at previous accepted step
    t_prev = Ref(NaN)     # elapsed integration time at that step

    # g(u,t) = calc_value - target. Mutates subject as a side effect (unavoidable
    # given the get_calc(calc) contract; only touches subject.state).
    function g_at(u)
        AstroProp._subject_update_from_u!(subject, dynsys, u)
        return get_calc(calc) - target
    end

    # Per-step: 1 direct u eval, no interpolant.
    function cond_fn(u, t, _integ)
        g_now = g_at(u)
        if isnan(g_prev[])
            g_prev[] = g_now; t_prev[] = t
            return false
        end
        crossed = dir < 0 ? (g_prev[] > 0 && g_now ≤ 0) :
                  dir > 0 ? (g_prev[] < 0 && g_now ≥ 0) :
                            (sign(g_prev[]) != sign(g_now))
        if !crossed
            g_prev[] = g_now; t_prev[] = t
        end
        return crossed
    end

    # One-time bisection on the Vern9 interpolant between t_prev and current t.
    function affect!(integ)
        tl, th = t_prev[], integ.t
        gl     = g_prev[]                 # sign known: opposite of gh
        iters  = 0
        t_mid  = 0.5 * (tl + th)
        u_mid  = integ(t_mid)             # Vern9 dense output at t_mid
        while iters < 60
            iters += 1
            gm = g_at(u_mid)
            if abs(gm) < 1e-9 || (th - tl) < 1e-6
                break
            end
            if (gl > 0) == (gm > 0)
                tl = t_mid; gl = gm
            else
                th = t_mid
            end
            t_mid = 0.5 * (tl + th)
            u_mid = integ(t_mid)
        end
        _ROOT_T[]     = t_mid
        _ROOT_U[]     = collect(u_mid)
        _ROOT_ITERS[] = iters
        terminate!(integ)
    end

    return DiscreteCallback(cond_fn, affect!; save_positions=(false, false))
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
    _ROOT_T[]     = NaN
    _ROOT_U[]     = zeros(6)
    _ROOT_ITERS[] = 0
    return propagate!(prop, sat, StopAt(sat, MeanSMA(), MEAN_SMA_TARGET; direction=-1))
end

# ── Warmup + reference run ───────────────────────────────────────────────────
println("Warming up (hybrid variant)...")
prop = make_prop(); sat = make_sat()
sol  = run_once(sat, prop)
overshoot_s = sol.t[end] - _ROOT_T[]
println("  ODE steps:            ", length(sol.t))
println("  Final t (accepted):   ", round(sol.t[end] / 86400.0; digits=6), " days")
println("  Root t (interpolated):", round(_ROOT_T[]   / 86400.0; digits=6), " days")
println("  Overshoot:            ", round(overshoot_s; digits=6), " s of orbit time")
println("  Bisection iterations: ", _ROOT_ITERS[])

# Verify root state gives MeanSMA ≈ target
sat_v = make_sat()
AstroModels.set_posvel!(sat_v, _ROOT_U[])
sma_at_root = get_calc(OrbitCalc(sat_v, MeanSMA()))
println("  MeanSMA at root_u:    ", round(sma_at_root; digits=6), " km  (target = $MEAN_SMA_TARGET)")
println("  Residual:             ", round(sma_at_root - MEAN_SMA_TARGET; digits=8), " km")

# ── Benchmark ────────────────────────────────────────────────────────────────
println("\nBenchmarking propagate! (hybrid) to MeanSMA = $(MEAN_SMA_TARGET) km ...")
b = @benchmark run_once(sat, prop) setup=(sat = make_sat(); prop = make_prop()) evals=1 samples=5 seconds=600

show(stdout, MIME"text/plain"(), b)
println()

t_median_s = median(b.times) / 1e9
mem_median = median(b.memory)
allocs_med = median(b.allocs)
println("\nSUMMARY: median = $(round(t_median_s; digits=3)) s, " *
        "memory = $(round(mem_median / 2^20; digits=2)) MiB, " *
        "allocs = $(allocs_med), ODE steps = $(length(sol.t))")

println("\nBENCH_HYBRID_CALLBACK_DONE")
