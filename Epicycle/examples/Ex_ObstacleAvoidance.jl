# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Obstacle Avoidance
#'
#' Steer a vehicle between two points while avoiding two circular keep-out zones. A
#' `HermiteSimpson` collocation phase solves the route and heading history.
#'
#' The keep-out zones are the point of the example: a clearance is a path constraint, checked at
#' every mesh point rather than at an endpoint. Speed is constant here, so the cost is a fixed
#' value and the constraint does all the work.

using Epicycle
using LinearAlgebra

#' ## Problem Formulation
#'
#' The path-planning problem is
#'
#' ```math
#' \begin{aligned}
#' \underset{x(t),\,y(t),\,\theta(t)}{\operatorname{minimize}} \quad
#' & \int_0^1 V^2\,dt \\
#' \text{subject to} \quad
#' & \dot{x}=V\cos\theta,\qquad \dot{y}=V\sin\theta, \\
#' & x(0)=0,\quad y(0)=0, \\
#' & x(1)=1.2,\quad y(1)=1.6, \\
#' & (x-x_{c,i})^2+(y-y_{c,i})^2 \ge r_i^2,\qquad i=1,2.
#' \end{aligned}
#' ```
#'
#' where ``x`` and ``y`` are position, ``\theta`` is heading, and
#' ``V=2.138`` is the constant speed. The obstacle centers are ``(0.4,0.5)``
#' and ``(0.8,1.5)``, and both radii are ``\pi/10``.
#'
#' Because the duration and speed are fixed, the objective is constant. The
#' endpoint and clearance constraints determine the feasible route.

#' ## Configuration
#'
#' Define the cruise speed, the obstacles, and the two ends of the route.

# Set the vehicle and course geometry
const V = 2.138                       # cruise speed, held constant
const R2 = (π / 10)^2                 # squared radius of each keep-out zone
const XC1, YC1 = 0.4, 0.5             # centre of the first obstacle
const XC2, YC2 = 0.8, 1.5             # centre of the second obstacle
const X0, Y0 = 0.0, 0.0               # start
const XF, YF = 1.2, 1.6               # finish

# Define the state and control representations
struct OAState{T} <: AbstractState
    x::T
    y::T
end

struct OAControl{T} <: AbstractControl
    θ::T
end

#' ## Equations of motion
#'
#' The vehicle cruises at constant speed, so the heading is the only control.

# Evaluate the position rates
function oa_dynamics!(dy, y::OAState, u::OAControl, p, t, model)
    dy[1] = V * cos(u.θ)
    dy[2] = V * sin(u.θ)
end

# Supply the control partials
@partial(oa_dynamics!, control) do dF, y, u, p, t, model
    dF[1, 1] = -V * sin(u.θ)
    dF[2, 1] = V * cos(u.θ)
end

#' ## Configure the optimal control problem
#'
#' The guess is a set of waypoints around the obstacles, with the heading between consecutive
#' waypoints as the control guess, so the solver starts pointing the way it has to go.

# Build the Hermite-Simpson phase
phase = CollocationPhase(name = :obstacle_avoidance,
                         transcription = HermiteSimpson(n_steps = 35),
                         dynamics = oa_dynamics!,
                         state = OAState,
                         control = OAControl,
                         tspan = (0.0, 1.0))

# Route the initial guess around the obstacles
const WX = [X0, 0.3, 0.6, 0.9, XF]
const WY = [Y0, 0.4, 0.8, 1.2, YF]

const THETA = let θ = [atan(WY[i+1] - WY[i], WX[i+1] - WX[i]) for i in 1:4]
    reshape(push!(θ, θ[end]), 1, 5)
end

# Vary the position and heading histories
Vary(state, phase;
     guess = [WX'; WY'],
     lower_bound = [0.0, 0.0],
     upper_bound = [XF, YF],
     scale = [XF, YF])

Vary(control, phase;
     guess = THETA,
     lower_bound = [-10.0],
     upper_bound = [10.0])

# Define the boundary, clearance, and objective quantities
position(c) = [state(c).x, state(c).y]

clearance(c) = [(state(c).x - XC1)^2 + (state(c).y - YC1)^2,
                (state(c).x - XC2)^2 + (state(c).y - YC2)^2]

speed_squared(c) = (V * cos(control(c).θ))^2 + (V * sin(control(c).θ))^2

# Supply their analytic partials
@partial(position, state) do c
    Matrix{Float64}(I, 2, 2)
end

@partial(clearance, state) do c
    [2 * (state(c).x - XC1)  2 * (state(c).y - YC1)
     2 * (state(c).x - XC2)  2 * (state(c).y - YC2)]
end

@partial(speed_squared, control) do c
    θ = control(c).θ
    [-2 * V^2 * cos(θ) * sin(θ) + 2 * V^2 * sin(θ) * cos(θ)]
end

# Fix both endpoints and enforce obstacle clearance along the path
Constraint(position, phase; equals = [X0, Y0], at = Initial())
Constraint(position, phase; equals = [XF, YF], at = Final())
Constraint(clearance, phase; lower_bound = [R2, R2], at = Path())

# Minimize the integrated squared speed
Objective(speed_squared, phase; sense = Min(), at = Path())

#' ## Solve the optimal control problem

# Check the declared partials against finite differences
check_partials(phase)

# Solve and report the arrival and cost
result = solve!(Sequence(phase); method = Optimize(print_level = 5))
yf = get_final_state(phase)
println("status    : ", result.info)
println("x_f       : ", round(yf.x, digits = 6), "   (target ", XF, ")")
println("y_f       : ", round(yf.y, digits = 6), "   (target ", YF, ")")
println("objective : ", round(result.objective, digits = 6),
        "   (constant speed, so V^2 tf = ", round(V^2, digits = 6), ")")
