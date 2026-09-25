# Copyright (C) 2026 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Station Keeping
#'
#' Hold a low Earth orbit against drag for one year by re-boosting whenever its mean semi-major axis
#' reaches a lower threshold. Each cycle coasts to the trigger, then solves a two-burn Hohmann
#' transfer back to the nominal orbit.
#'
#' The example uses the exponential atmosphere so that it runs anywhere. The MSISE model, which
#' EpicycleEnterprise provides, is what a comparison against another tool would use.

using Epicycle

#' ## Configuration
#'
#' Define the mission thresholds, the spacecraft with its drag geometry, and the propagator.

# Set the nominal orbit, re-boost trigger, and mission duration
const MEAN_SMA_NOMINAL = 6778.0             # km, the mission mean semi-major axis
const MEAN_SMA_TRIGGER = 6768.0             # km, where a re-boost is due
const SIM_BUDGET_SEC = 365.25 * 86400       # one year

# Configure the spacecraft mass and drag geometry
sat = Spacecraft(state = KeplerianState(6778.137, 0.0, deg2rad(51.6), 0.0, 0.0, 0.0),
                 time = Time("2020-01-01T00:00:00", UTC(), ISOT()),
                 name = "SKSat",
                 mass = 1000.0,
                 drag = SphericalDrag(c_d = 2.2, drag_area = 10.0))

# Configure point-mass gravity and atmospheric drag
gravity = PointMassGravity(earth, ())
drag = AtmosphericDrag(earth; model = Exponential())
forces = ForceModel(gravity, drag)
integ = IntegratorConfig(Vern9();
                         reltol = 1e-10,
                         abstol = 1e-10,
                         dt = 60.0)
prop = OrbitPropagator(forces, integ)

# Record the epoch used to measure elapsed mission time
const T_START_TT_JD = sat.time.tt.jd

elapsed_sec(sc) = (sc.time.tt.jd - T_START_TT_JD) * 86400.0

#' ## The re-boost problem
#'
#' A re-boost is the Hohmann transfer of the targeting examples, solved where the spacecraft
#' happens to be. The analytic two-burn solution is the guess, and the mean semi-major axis after
#' the second burn is the constraint.

# Compute the analytic Hohmann burns between two circular radii
function hohmann_dvs(r1, r2, μ)
    a_trans = 0.5 * (r1 + r2)
    v1_c = sqrt(μ / r1)
    v2_c = sqrt(μ / r2)
    vp_t = sqrt(μ * (2 / r1 - 1 / a_trans))
    va_t = sqrt(μ * (2 / r2 - 1 / a_trans))
    return (vp_t - v1_c, v2_c - va_t)
end

# Solve one re-boost at the spacecraft's current state
function reboost!(sat, prop; target_mean_sma)

    # Seed both burns with the analytic Hohmann solution
    dv1_guess, dv2_guess = hohmann_dvs(mean_long_sma(sat), target_mean_sma, earth.mu)

    # Seed the two VNB maneuvers
    toi = ImpulsiveManeuver(axes = VNB(),
                            element1 = dv1_guess,
                            element2 = 0.0,
                            element3 = 0.0)

    moi = ImpulsiveManeuver(axes = VNB(),
                            element1 = dv2_guess,
                            element2 = 0.0,
                            element3 = 0.0)

    # Assemble the re-boost sequence
    seq = Sequence()
    add_sequence!(seq,

        # Apply the orbit-raising burn, capped at 50 m/s
        Event(name = "TOI",
              event = () -> maneuver!(sat, toi),
              vars = [Vary(delta_v, toi;
                           lower_bound = [0.0, 0.0, 0.0],
                           upper_bound = [0.05, 0.0, 0.0],
                           name = "toi")]),

        # Coast to transfer-orbit apoapsis
        Event(name = "Coast to apoapsis",
              event = () -> propagate!(prop, sat,
                                       StopAt(position_dot_velocity, sat;
                                              equals = 0.0,
                                              direction = -1))),

        # Circularize and constrain the restored mean semi-major axis
        Event(name = "MOI",
              event = () -> maneuver!(sat, moi),
              vars = [Vary(delta_v, moi;
                           lower_bound = [0.0, 0.0, 0.0],
                           upper_bound = [0.05, 0.0, 0.0],
                           name = "moi")],
              funcs = [Constraint(mean_long_sma, sat; equals = target_mean_sma)]))

    # Solve the re-boost
    solve!(seq; method = Optimize(derivatives = :fd, print_level = 5))

    return (dv1_guess, dv2_guess, toi.element1, moi.element1)
end

#' ## The station keeping loop
#'
#' Each cycle coasts until the mean semi-major axis reaches the trigger, or until the simulation
#' budget runs out, whichever comes first. Two stopping conditions in one call are what express
#' that: the propagation ends on whichever is met first.

# Run the mission, one re-boost cycle at a time
cycle = 0
log = NamedTuple[]

while true
    t_remain = SIM_BUDGET_SEC - elapsed_sec(sat)
    t_remain > 0.0 || break

    # Coast to the trigger, or to the end of the simulation
    propagate!(prop, sat,
               StopAt(mean_long_sma, sat; equals = MEAN_SMA_TRIGGER, direction = -1),
               StopAt(sat, PropDurationSeconds(), t_remain))

    # Stop when the budget, rather than the trigger, ended the coast
    elapsed_sec(sat) >= SIM_BUDGET_SEC - 1.0 && break

    # Solve the re-boost and record what it took
    sma_before = mean_long_sma(sat)
    dv1_guess, dv2_guess, dv1, dv2 = reboost!(sat, prop; target_mean_sma = MEAN_SMA_NOMINAL)

    global cycle += 1
    push!(log, (cycle = cycle,
                day = elapsed_sec(sat) / 86400.0,
                sma_before = sma_before,
                sma_after = mean_long_sma(sat),
                dv1_mps = 1e3 * dv1,
                dv2_mps = 1e3 * dv2))
end

#' ## Report the mission
#'
#' What the year cost: the re-boosts flown, when each was due, and the total delta-v spent.

# Report each cycle and the total delta-v
for c in log
    println("cycle ", c.cycle,
            ": day ", round(c.day, digits = 2),
            "   sma ", round(c.sma_before, digits = 3), " -> ", round(c.sma_after, digits = 3), " km",
            "   delta-v ", round(c.dv1_mps, digits = 3), " + ", round(c.dv2_mps, digits = 3), " m/s")
end

total_dv = sum(c.dv1_mps + c.dv2_mps for c in log; init = 0.0)
println("re-boosts : ", cycle, " over ", round(elapsed_sec(sat) / 86400.0, digits = 2), " days")
println("total dv  : ", round(total_dv, digits = 3), " m/s")
