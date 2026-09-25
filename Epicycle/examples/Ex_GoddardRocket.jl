# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # The Goddard Rocket Problem
#'
#' Maximize the altitude reached by a vertical sounding rocket with drag and
#' limited propellant. The solution uses three Hermite-Simpson phases connected
#' by state and time continuity constraints.
#'
#' The optimal control has three arcs: full thrust, a singular-thrust arc, and a
#' coast to apogee. Separate phases allow each arc to use its own thrust bounds
#' and boundary conditions.

using Epicycle

#' ## Problem Formulation
#'
#' The maximum-altitude problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{h(t),\,v(t),\,m(t),\,T(t),\,t_f}{\operatorname{maximize}}
#' \quad & h(t_f) \\
#' \text{subject to} \quad
#' & \dot{h} = v, \\
#' & \dot{v} = \frac{T-D(h,v)}{m} - g, \\
#' & \dot{m} = -\frac{T}{c}, \\
#' & D(h,v) = D_0 v^2 e^{-Bh}, \\
#' & h(0)=0,\quad v(0)=0,\quad m(0)=3, \\
#' & v(t_f)=0,\quad m(t_f)=1, \\
#' & 0 \le T(t) \le T_{\max}.
#' \end{aligned}
#' ```
#'
#' where ``h`` is altitude, ``v`` is vertical speed, ``m`` is mass, ``T`` is
#' thrust, ``D`` is atmospheric drag, ``g`` is gravitational acceleration,
#' and ``c`` is effective exhaust velocity.
#'
#' The optimal thrust history contains full-thrust, singular-thrust, and coast
#' arcs. The three phases below impose that structure and solve for the two
#' switching times.

#' ## Configuration
#'
#' Define the gravity, atmosphere, propulsion model, rocket state, and thrust
#' control. The state contains altitude, vertical speed, and mass.

# Define the rocket and atmosphere model
struct GRKModel
    g::Float64                # ft/s^2
    D0::Float64               # drag coefficient
    B::Float64                # inverse density scale height, 1/ft
    c::Float64                # effective exhaust velocity, ft/s
    T_max::Float64            # lbf
end

const grk = GRKModel(32.174,
                     5.49153484923381e-5,
                     1.0 / 23800.0,
                     sqrt(3.264 * 32.174 * 23800.0),
                     193.044)

# Define the state and control variables
struct GRKState{T} <: AbstractState
    h::T
    v::T
    m::T
end

struct GRKControl{T} <: AbstractControl
    T_::T
end

#' ## Equations of Motion
#'
#' Use the same vertical dynamics in every phase. Atmospheric drag decreases
#' exponentially with altitude, and propellant flow is proportional to thrust.

# Evaluate altitude, speed, and mass rates
function grk_dynamics!(dy, y::GRKState, u::GRKControl, p, t, md)
    D = md.D0 * y.v^2 * exp(-md.B * y.h)
    dy[1] = y.v
    dy[2] = (u.T_ - D) / y.m - md.g
    dy[3] = -u.T_ / md.c
end

# Supply the dynamics partial with respect to state
@partial(grk_dynamics!, state) do dF, y, u, p, t, md
    D = md.D0 * y.v^2 * exp(-md.B * y.h)
    dF[1, 2] = 1.0
    dF[2, 1] = md.B * D / y.m
    dF[2, 2] = -2.0 * md.D0 * y.v * exp(-md.B * y.h) / y.m
    dF[2, 3] = -(u.T_ - D) / y.m^2
end

# Supply the dynamics partial with respect to control
@partial(grk_dynamics!, control) do dF, y, u, p, t, md
    dF[2, 1] = 1.0 / y.m
    dF[3, 1] = -1.0 / md.c
end

#' ## Build the Initial Guess
#'
#' Generate a dynamically consistent guess by integrating the expected control
#' structure with RK4: full thrust to the first switch, singular thrust to
#' propellant depletion, and zero thrust to apogee. The propagated switching
#' times and states initialize the three phases.

# Define drag, singular thrust, and the switching condition
grk_drag(h, v) = grk.D0 * v^2 * exp(-grk.B * h)

function grk_singular_thrust(h, v, m)
    D = grk_drag(h, v)
    vc = v / grk.c
    mg = m * grk.g
    t1 = grk.c^2 * (1.0 + vc) * grk.B / grk.g - 1.0 - 2.0 / vc
    t2 = mg / (1.0 + 4.0 / vc + 2.0 / vc^2)
    return clamp(D + mg + t1 * t2, 0.0, grk.T_max)
end

grk_switch(h, v, m) = m * grk.g - (1.0 + v / grk.c) * grk_drag(h, v)

# Propagate the three arcs and return their boundary states and times
function grk_flown(; dt = 0.01)
    f(y, T) = [y[2], (T - grk_drag(y[1], y[2])) / y[3] - grk.g, -T / grk.c]
    step(y, T) = (k1 = f(y, T); k2 = f(y .+ dt / 2 .* k1, T); k3 = f(y .+ dt / 2 .* k2, T);
                  k4 = f(y .+ dt .* k3, T); y .+ dt / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4))
    y, t = [0.0, 0.0, 3.0], 0.0

    # Apply full thrust to the first switch
    while grk_switch(y...) > 0 && y[3] > 1.0
        y, t = step(y, grk.T_max), t + dt
    end
    y1, t1 = copy(y), t

    # Follow the singular-thrust law to propellant depletion
    while y[3] > 1.0
        y, t = step(y, grk_singular_thrust(y...)), t + dt
    end
    y[3] = 1.0
    y2, t2 = copy(y), t

    # Coast to apogee
    while y[2] > 0
        y, t = step(y, 0.0), t + dt
    end
    return (y1, t1), (y2, t2), (copy(y), t)
end

(y1, t1), (y2, t2), (y3, t3) = grk_flown()
println("flown guess: peak altitude ", round(y3[1], digits = 1), " ft at t = ", round(t3, digits = 2), " s")

#' ## Build the Three Phases
#'
#' Build each phase with the same dynamics and different thrust bounds. The
#' propagated arc endpoints provide state and time guesses.

# Configure one ascent phase
function grk_phase(name; n_steps, h_wp, v_wp, m_wp, T_lb, T_ub, T_guess = (T_lb + T_ub) / 2, t0, tf)
    ph = CollocationPhase(name = name,
                          transcription = HermiteSimpson(n_steps = n_steps),
                          dynamics = grk_dynamics!,
                          model = grk,
                          state = GRKState,
                          control = GRKControl,
                          tspan = (Float64(t0), Float64(tf)))

    Vary(state, ph;
         guess = [h_wp[1] h_wp[2]; v_wp[1] v_wp[2]; m_wp[1] m_wp[2]],
         lower_bound = [0.0, 0.0, 1.0],
         upper_bound = [1e10, 1e10, 3.0])

    Vary(control, ph;
         guess = fill(T_guess, 1, 2),
         lower_bound = [T_lb],
         upper_bound = [T_ub])

    return ph
end

# Create the boost, singular-thrust, and coast phases
boost = grk_phase(:boost;
                  n_steps = 12, t0 = 0.0, tf = t1,
                  h_wp = [0.0, y1[1]], v_wp = [0.0, y1[2]], m_wp = [3.0, y1[3]],
                  T_lb = grk.T_max, T_ub = grk.T_max)

singular = grk_phase(:singular;
                     n_steps = 12, t0 = t1, tf = t2,
                     h_wp = [y1[1], y2[1]], v_wp = [y1[2], y2[2]], m_wp = [y1[3], y2[3]],
                     T_lb = 0.0, T_ub = grk.T_max, T_guess = grk_singular_thrust(y1...))

coast = grk_phase(:coast;
                  n_steps = 12, t0 = t2, tf = t3,
                  h_wp = [y2[1], y3[1]], v_wp = [y2[2], 0.0], m_wp = [1.0, 1.0],
                  T_lb = 0.0, T_ub = 0.0)

# Vary the phase boundaries so the solver determines the switching times
Vary(final_time, boost; guess = t1, lower_bound = 0.0, upper_bound = 100.0)
Vary(initial_time, singular; guess = t1, lower_bound = 0.0, upper_bound = 100.0)
Vary(final_time, singular; guess = t2, lower_bound = 0.0, upper_bound = 100.0)
Vary(initial_time, coast; guess = t2, lower_bound = 0.0, upper_bound = 100.0)
Vary(final_time, coast; guess = t3, lower_bound = 0.1, upper_bound = 100.0)

#' ## Apply the Phase Constraints
#'
#' Fix the liftoff state, enforce the singular-thrust law throughout the middle
#' phase, and end the final coast at apogee with the propellant depleted.

# Fix the liftoff state
liftoff(c) = [state(c).h, state(c).v, state(c).m]
Constraint(liftoff, boost; equals = [0.0, 0.0, 3.0], at = Initial())

# Enforce Betts's singular-arc condition
function singular_arc(c)
    y = state(c)
    D = grk.D0 * y.v^2 * exp(-grk.B * y.h)
    vc = y.v / grk.c
    mg = y.m * grk.g
    t1 = grk.c^2 * (1.0 + vc) * grk.B / grk.g - 1.0 - 2.0 / vc
    t2 = mg / (1.0 + 4.0 / vc + 2.0 / vc^2)
    [control(c).T_ - D - mg - t1 * t2]
end
Constraint(singular_arc, singular; equals = [0.0], at = Path())

# Locate the transition from singular thrust to coast
function switching(c)
    y = state(c)
    D = grk.D0 * y.v^2 * exp(-grk.B * y.h)
    [y.m * grk.g - (1.0 + y.v / grk.c) * D]
end
Constraint(switching, singular; equals = [0.0], at = Final())

# End the coast at apogee with dry mass remaining
burnout(c) = [state(c).v, state(c).m]
Constraint(burnout, coast; equals = [0.0, 1.0], at = Final())

#' ## Join the Phases
#'
#' Apply `continuity` to each `Link` so state and time agree across both phase
#' boundaries.

# Assemble and connect the phase sequence
seq = Sequence()
add_sequence!(seq, boost)
add_sequence!(seq, singular)
add_sequence!(seq, coast)

Constraint(continuity, Link(boost, singular; name = :boost_to_singular))
Constraint(continuity, Link(singular, coast; name = :singular_to_coast))

# Maximize altitude at the end of the coast
peak_altitude(c) = state(c).h
Objective(peak_altitude, coast; sense = Max())

#' ## Solve the Optimal Control Problem
#'
#' Solve the three-phase nonlinear program and compare its peak altitude with
#' the trajectory used to initialize the phases.

# Solve and report the peak altitude
result = solve!(seq; method = Optimize(max_iter = 3000, tol = 1e-5, print_level = 5))
println("status        : ", result.info)
println("peak altitude : ", round(get_final_state(coast).h, digits = 2),
        " ft   (the same control law flown: ", round(y3[1], digits = 1), " ft)")
