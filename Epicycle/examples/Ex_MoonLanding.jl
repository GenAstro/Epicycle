# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Soft Lunar Landing
#'
#' Land on the Moon with as much propellant remaining as possible. A `HermiteSimpson` collocation
#' phase solves the altitude, velocity, mass, thrust history, and landing time.
#'
#' The problem is Meditch's, in canonical units where the lunar gravity is one. The answer is a
#' bang-bang control: coast, then full thrust to touchdown. Dymos reports a final mass of 0.3953
#' after 1.397 units of time.

using Epicycle

#' ## Problem Formulation
#'
#' The maximum-final-mass problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{h(t),\,v(t),\,m(t),\,T(t),\,t_f}{\operatorname{maximize}}
#' \quad & m(t_f) \\
#' \text{subject to} \quad
#' & \dot{h}=v, \\
#' & \dot{v}=-1+\frac{T}{m}, \\
#' & \dot{m}=-\frac{T}{c}, \\
#' & h(0)=1,\quad v(0)=-0.783,\quad m(0)=1, \\
#' & h(t_f)=0,\quad v(t_f)=0, \\
#' & 0 \le T(t) \le 1.227.
#' \end{aligned}
#' ```
#'
#' where ``h`` is altitude, ``v`` is vertical velocity, ``m`` is mass, ``T``
#' is thrust, and ``c=2.349`` is the mass-flow coefficient. Time and the
#' remaining quantities are expressed in the canonical units of the problem.

#' ## Configuration
#'
#' Use altitude, velocity, and mass as the state, with thrust as the control.

# Set the mass-flow coefficient and thrust limit
struct MoonModel
    mdot_coeff::Float64
    T_max::Float64
end

model = MoonModel(2.349, 1.227)

# Define the state and control representations
struct MoonState{T} <: AbstractState
    h::T
    v::T
    m::T
end

struct MoonControl{T} <: AbstractControl
    T_::T
end

#' ## Equations of motion
#'
#' Altitude changes with velocity, velocity with the thrust acceleration against gravity, and mass
#' with the thrust.

# Evaluate the state rates
function moon_dynamics!(dy, y::MoonState, u::MoonControl, p, t, model)
    dy[1] = y.v
    dy[2] = -1.0 + u.T_ / y.m
    dy[3] = -u.T_ / model.mdot_coeff
end

# Supply the state partials
@partial(moon_dynamics!, state) do dF, y, u, p, t, model
    dF[1, 2] = 1.0
    dF[2, 3] = -u.T_ / y.m^2
end

# Supply the control partials
@partial(moon_dynamics!, control) do dF, y, u, p, t, model
    dF[2, 1] = 1.0 / y.m
    dF[3, 1] = -1.0 / model.mdot_coeff
end

#' ## Configure the optimal control problem
#'
#' The lander starts at altitude one, descending, with full mass, and must reach the surface at
#' rest. The final time is free, so the descent finds its own duration.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :moon_landing,
                         transcription = HermiteSimpson(n_steps = 30),
                         dynamics = moon_dynamics!,
                         model = model,
                         state = MoonState,
                         control = MoonControl,
                         tspan = (0.0, 1.4))

# Vary the state history from the endpoint guesses
Vary(state, phase;
     guess = [1.0 0.0
              -0.783 0.0
              1.0 0.4],
     lower_bound = [0.0, -5.0, 0.001],
     upper_bound = [5.0, 5.0, 2.0])

# Vary the thrust history and landing time
Vary(control, phase;
     guess = reshape([model.T_max / 2 model.T_max / 2], 1, 2),
     lower_bound = [0.0],
     upper_bound = [model.T_max])

Vary(final_time, phase;
     guess = 1.4,
     lower_bound = 0.5,
     upper_bound = 5.0)

# Define the boundary quantities and objective
launch_state(c) = [state(c).h, state(c).v, state(c).m]

touchdown(c) = [state(c).h, state(c).v]

final_mass(c) = state(c).m

# Supply their analytic partials
@partial(launch_state, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0
     0.0 0.0 1.0]
end

@partial(touchdown, state) do c
    [1.0 0.0 0.0
     0.0 1.0 0.0]
end

@partial(final_mass, state) do c
    [0.0 0.0 1.0]
end

# Fix the initial state and require touchdown at rest
Constraint(launch_state, phase; equals = [1.0, -0.783, 1.0], at = Initial())
Constraint(touchdown, phase; equals = [0.0, 0.0], at = Final())

# Maximize the final mass
Objective(final_mass, phase; sense = Max())

#' ## Solve the optimal control problem

# Check the declared partials against finite differences
check_partials(phase)

# Solve and compare the landing with the published values
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
yf = get_final_state(phase)
println("status     : ", result.info)
println("altitude   : ", round(yf.h, digits = 6), "   (touchdown 0)")
println("velocity   : ", round(yf.v, digits = 6), "   (at rest 0)")
println("final mass : ", round(yf.m, digits = 6), "   (Dymos 0.3953)")
println("tf         : ", round(get_final_time(phase), digits = 6), "   (Dymos 1.397)")
