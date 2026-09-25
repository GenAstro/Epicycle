# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Low-Thrust Orbit Raising with Collocation
#'
#' Maximize the radius reached by a continuously thrusting spacecraft in a fixed time. A
#' `HermiteSimpson` collocation phase solves the state history and thrust direction.
#'
#' This is Bryson's orbit-raising problem in canonical units, with an initial radius and
#' gravitational parameter of one. CSALT reports a final radius of 1.5230.

using Epicycle
using LinearAlgebra

#' ## Problem Formulation
#'
#' The fixed-time orbit-raising problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{r(t),\,\theta(t),\,v_r(t),\,v_\theta(t),\,m(t),
#' \,u_r(t),\,u_\theta(t)}{\operatorname{maximize}} \quad & r(t_f) \\
#' \text{subject to} \quad
#' & \dot{r}=v_r, \\
#' & \dot{\theta}=\frac{v_\theta}{r}, \\
#' & \dot{v}_r=\frac{v_\theta^2}{r}-\frac{\mu}{r^2}
#'   +\frac{T}{m}u_r, \\
#' & \dot{v}_\theta=-\frac{v_rv_\theta}{r}+\frac{T}{m}u_\theta, \\
#' & \dot{m}=-\dot{m}_p, \\
#' & r(0)=1,\quad \theta(0)=0,\quad v_r(0)=0,
#'   \quad v_\theta(0)=1,\quad m(0)=1, \\
#' & v_r(t_f)=0,\quad v_\theta(t_f)=\sqrt{\frac{\mu}{r(t_f)}}, \\
#' & u_r^2+u_\theta^2=1,\quad t_f=3.32.
#' \end{aligned}
#' ```
#'
#' where ``r`` and ``\theta`` are polar position, ``v_r`` and ``v_\theta``
#' are radial and transverse velocity, ``m`` is mass, and ``u_r`` and
#' ``u_\theta`` are the radial and transverse components of thrust direction.
#' The canonical model uses ``\mu=1``, ``T=0.1405``, and
#' ``\dot{m}_p=0.0749``.

#' ## Configuration
#'
#' Use constant thrust and mass flow, leaving thrust direction as the control.

# Set the gravitational parameter, thrust, and mass flow
struct RaisingModel
    mu    ::Float64
    thrust::Float64
    mdot  ::Float64
end

model = RaisingModel(1.0, 0.1405, 0.0749)

# Set the initial state and transfer duration
y0 = [1.0, 0.0, 0.0, 1.0, 1.0]        # radius, angle, radial and transverse speed, mass
tf = 3.32

#' ## Equations of motion
#'
#' The state is radius, polar angle, radial speed, transverse speed and mass. The control is the
#' thrust direction, given by its radial and transverse components.

# Define the state and control representations
struct RaisingState{T} <: AbstractState
    r ::T
    θ ::T
    vr::T
    vθ::T
    m ::T
end

struct RaisingControl{T} <: AbstractControl
    u_r::T
    u_θ::T
end

# Evaluate the state rates
function raising!(dy, y::RaisingState, u::RaisingControl, p, t, model)
    a = model.thrust / y.m
    dy[1] = y.vr
    dy[2] = y.vθ / y.r
    dy[3] = y.vθ^2 / y.r - model.mu / y.r^2 + a * u.u_r
    dy[4] = -y.vr * y.vθ / y.r + a * u.u_θ
    dy[5] = -model.mdot
end

# Supply the state partials
@partial(raising!, state) do dF, y, u, p, t, model
    a = model.thrust / y.m
    dF[1, 3] = 1.0
    dF[2, 1] = -y.vθ / y.r^2
    dF[2, 4] = 1.0 / y.r
    dF[3, 1] = -y.vθ^2 / y.r^2 + 2.0 * model.mu / y.r^3
    dF[3, 4] = 2.0 * y.vθ / y.r
    dF[3, 5] = -a * u.u_r / y.m
    dF[4, 1] = y.vr * y.vθ / y.r^2
    dF[4, 3] = -y.vθ / y.r
    dF[4, 4] = -y.vr / y.r
    dF[4, 5] = -a * u.u_θ / y.m
end

# Supply the control partials
@partial(raising!, control) do dF, y, u, p, t, model
    a = model.thrust / y.m
    dF[3, 1] = a
    dF[4, 2] = a
end

#' ## The initial guess
#'
#' Generate a dynamically consistent guess by integrating transverse thrust with RK4 and sampling
#' the result at the mesh nodes.

# Fly the initial guess
function flown_guess(times; steps = 2000)
    f(y) = [y[3],
            y[4] / y[1],
            y[4]^2 / y[1] - model.mu / y[1]^2,
            -y[3] * y[4] / y[1] + model.thrust / y[5],
            -model.mdot]
    h = tf / steps
    y = copy(y0)
    t = 0.0
    out = zeros(5, length(times))
    next = 1
    for _ in 0:steps
        while next <= length(times) && t >= times[next] - 1e-9
            out[:, next] = y
            next += 1
        end
        k1 = f(y)
        k2 = f(y .+ h / 2 .* k1)
        k3 = f(y .+ h / 2 .* k2)
        k4 = f(y .+ h .* k3)
        y = y .+ h / 6 .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        t += h
    end
    return out
end

guess = flown_guess([0.0, tf / 2, tf])

#' ## Configure the optimal control problem
#'
#' Fix the departure state, keep the thrust direction at unit magnitude, and require a circular
#' arrival orbit. The objective maximizes the final radius.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :orbit_raising,
                         transcription = HermiteSimpson(n_steps = 40),
                         dynamics = raising!,
                         model = model,
                         state = RaisingState,
                         control = RaisingControl,
                         tspan = (0.0, tf))

# Vary the state history using the flown guess
Vary(state, phase;
     guess = guess,
     lower_bound = [0.5, 0.0, -10.0, -10.0, 0.1],
     upper_bound = [5.0, 4π, 10.0, 10.0, 3.0])

# Vary the thrust direction
Vary(control, phase;
     guess = fill(1.0 / sqrt(2), 2, 3),
     lower_bound = [-10.0, -10.0],
     upper_bound = [10.0, 10.0])

# Define the boundary, path, and objective quantities
departure(c) = (y = state(c); [y.r, y.θ, y.vr, y.vθ, y.m])

unit_thrust(c) = [control(c).u_r^2 + control(c).u_θ^2]

circular_orbit(c) = (y = state(c); [y.vr, sqrt(c.model.mu / y.r) - y.vθ])

final_radius(c) = state(c).r

# Supply their analytic partials
@partial(departure, state) do c
    Matrix(1.0I, 5, 5)
end

@partial(unit_thrust, control) do c
    [2.0 * control(c).u_r  2.0 * control(c).u_θ]
end

@partial(circular_orbit, state) do c
    y = state(c)
    [0.0                                     0.0  1.0   0.0  0.0
     -0.5 * sqrt(c.model.mu / y.r) / y.r     0.0  0.0  -1.0  0.0]
end

@partial(final_radius, state) do c
    [1.0  0.0  0.0  0.0  0.0]
end

# Apply the departure, thrust-direction, and arrival constraints
Constraint(departure, phase; equals = y0, at = Initial())
Constraint(unit_thrust, phase; equals = 1.0, at = Path())
Constraint(circular_orbit, phase; equals = [0.0, 0.0], at = Final())

# Maximize the final radius
Objective(final_radius, phase; sense = Max())

#' ## Solve the optimal control problem
#'
#' `check_partials` compares every declared partial against a finite difference before the solve,
#' which is where a wrong derivative is cheap to find.

# Check the declared partials against finite differences
check_partials(phase)

# Solve the problem
result = solve!(Sequence(phase); method = Optimize(max_iter = 1000, tol = 1e-6, print_level = 5))

# Compare the solution with the CSALT reference
y_final = get_final_state(phase)
println("status       : ", result.info)
println("final radius : ", round(y_final.r, digits = 6), "   (CSALT reference 1.5230)")
println("final mass   : ", round(y_final.m, digits = 6))
println("radial speed : ", round(y_final.vr, sigdigits = 3), "   (circular orbit: 0)")
