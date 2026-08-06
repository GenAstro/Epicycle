# Station-keeping via periodic Hohmann re-boosts.
#
# Nominal circular LEO at mean SMA = 6778 km. Drag pulls the mean SMA down;
# when it drops to 6768 km (10 km), a 2-burn Hohmann raises it back to 6778 km.
# Analytic Hohmann ΔVs are used as the initial guess to the solver, and the
# solver refines them against a mean-SMA constraint at MOI.
#
# Loop terminates after 1 year of simulated time (365.25 days).
#
# NOTE: uses the exponential atmosphere for portability. For GMAT comparison,
# switch to MSISE (see EpicycleEnterprise/test/drag_msise_gmat.jl) and match
# the space-weather file.

using Epicycle
using OrdinaryDiffEq: Vern9

# ============================================================================
# Constants
# ============================================================================

const MEAN_SMA_NOMINAL = 6778.0          # km, mission mean SMA
const MEAN_SMA_TRIGGER = 6768.0          # km, decay stop threshold
const SIM_BUDGET_SEC   = 365.25 * 86400  # 1 year

# ============================================================================
# Spacecraft, forces, propagator
# ============================================================================

sat = Spacecraft(
    state = KeplerianState(6778.137, 0.0, deg2rad(51.6), 0.0, 0.0, 0.0),
    time  = Time("2020-01-01T00:00:00", UTC(), ISOT()),
    name  = "SKSat",
    mass  = 1000.0,
    drag  = SphericalDrag(c_d = 2.2, drag_area = 10.0),
)

const T_START_TT_JD = sat.time.tt.jd    # frozen at t0 for elapsed-time bookkeeping

gravity = PointMassGravity(earth, ())
drag    = AtmosphericDrag(earth; model = Exponential())
forces  = ForceModel(gravity, drag)
integ   = IntegratorConfig(Vern9(); reltol = 1e-10, abstol = 1e-10, dt = 60.0)
prop    = OrbitPropagator(forces, integ)

# ============================================================================
# Helpers
# ============================================================================

# Elapsed sim time since t0, in seconds (dynamical time — matches propagator).
elapsed_sec(sc) = (sc.time.tt.jd - T_START_TT_JD) * 86400.0

# Analytic Hohmann ΔVs (tangential) from circular r1 to circular r2.
function hohmann_dvs(r1, r2, μ)
    a_trans = 0.5 * (r1 + r2)
    v1_c    = sqrt(μ / r1)
    v2_c    = sqrt(μ / r2)
    vp_t    = sqrt(μ * (2/r1 - 1/a_trans))
    va_t    = sqrt(μ * (2/r2 - 1/a_trans))
    return (vp_t - v1_c, v2_c - va_t)
end

# Build and solve one Hohmann re-boost sequence at the current spacecraft state.
# Uses analytic ΔVs as the initial guess for both maneuvers.
function reboost!(sat, prop; target_mean_sma)

    # ------------------------------------------------------------------------
    # Analytic Hohmann guess — seeds the solver for both burns
    # ------------------------------------------------------------------------
    r1 = get_calc(OrbitCalc(sat, MeanSMA()))
    r2 = target_mean_sma
    μ  = earth.mu
    dv1_guess, dv2_guess = hohmann_dvs(r1, r2, μ)

    # ------------------------------------------------------------------------
    # TOI event — Transfer Orbit Insertion
    # ------------------------------------------------------------------------
    toi = ImpulsiveManeuver(axes = VNB(),
                            element1 = dv1_guess, element2 = 0.0, element3 = 0.0)

    toi_var = SolverVariable(
        calc        = ManeuverCalc(toi, sat, DeltaVVector()),
        name        = "toi",
        lower_bound = [0.0,  0.0, 0.0],
        upper_bound = [0.05, 0.0, 0.0],   # 50 m/s cap — comfortably above analytic
    )

    toi_fun()   = maneuver!(sat, toi)
    toi_event   = Event(name = "TOI", event = toi_fun, vars = [toi_var])

    # ------------------------------------------------------------------------
    # Coast event — propagate to apoapsis of the transfer orbit
    # ------------------------------------------------------------------------
    coast_fun() = propagate!(prop, sat, StopAt(sat, PosDotVel(), 0.0; direction=-1))
    coast_event = Event(name = "Coast to Apo", event = coast_fun)

    # ------------------------------------------------------------------------
    # MOI event — Mission Orbit Insertion, constrained to hit target mean SMA
    # ------------------------------------------------------------------------
    moi = ImpulsiveManeuver(axes = VNB(),
                            element1 = dv2_guess, element2 = 0.0, element3 = 0.0)

    moi_var = SolverVariable(
        calc        = ManeuverCalc(moi, sat, DeltaVVector()),
        name        = "moi",
        lower_bound = [0.0,  0.0, 0.0],
        upper_bound = [0.05, 0.0, 0.0],
    )

    sma_con = Constraint(
        calc         = OrbitCalc(sat, MeanSMA()),
        lower_bounds = [target_mean_sma],
        upper_bounds = [target_mean_sma],
        scale        = [1.0],
    )

    moi_fun()   = maneuver!(sat, moi)
    moi_event   = Event(name = "MOI", event = moi_fun,
                        vars = [moi_var], funcs = [sma_con])

    # ------------------------------------------------------------------------
    # Assemble and solve the sequence
    # ------------------------------------------------------------------------
    seq = Sequence()
    add_sequence!(seq, toi_event, coast_event, moi_event)
    solve_trajectory!(seq)

    return (dv1_guess, dv2_guess, toi.element1, moi.element1)
end

# ============================================================================
# Station-keeping loop
# ============================================================================

cycle = 0
history = NamedTuple[]   # per-cycle record for inspection
while true
    t_used     = elapsed_sec(sat)
    t_remain   = SIM_BUDGET_SEC - t_used
    t_remain > 0.0 || break

    # Coast until mean SMA drops to the trigger, OR the sim budget runs out —
    # whichever comes first.
    propagate!(prop, sat,
        StopAt(sat, MeanSMA(), MEAN_SMA_TRIGGER; direction=-1),
        StopAt(sat, PropDurationSeconds(), t_remain))

    # If the budget bounded us, we're done.
    elapsed_sec(sat) >= SIM_BUDGET_SEC - 1.0 && break

    sma_before = get_calc(OrbitCalc(sat, MeanSMA()))
    dv1g, dv2g, dv1, dv2 = reboost!(sat, prop; target_mean_sma = MEAN_SMA_NOMINAL)
    sma_after = get_calc(OrbitCalc(sat, MeanSMA()))

    global cycle += 1
    push!(history, (cycle=cycle, t_day=elapsed_sec(sat)/86400.0,
                    sma_before=sma_before, sma_after=sma_after,
                    dv1_guess_mps=1e3*dv1g, dv1_solved_mps=1e3*dv1,
                    dv2_guess_mps=1e3*dv2g, dv2_solved_mps=1e3*dv2))
    println("cycle $cycle: t=$(round(elapsed_sec(sat)/86400.0; digits=3)) d  " *
            "SMA $(round(sma_before; digits=3)) → $(round(sma_after; digits=3)) km  " *
            "Δv₁ $(round(1e3*dv1; digits=3)) (guess $(round(1e3*dv1g; digits=3)))  " *
            "Δv₂ $(round(1e3*dv2; digits=3)) (guess $(round(1e3*dv2g; digits=3))) m/s")
end

println("done: $(cycle) cycles over $(round(elapsed_sec(sat)/86400.0; digits=2)) days")
