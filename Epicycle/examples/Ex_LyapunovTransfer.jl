# Copyright (C) 2025 Gen Astro LLC
# SPDX-License-Identifier: LicenseRef-GenAstro-SourceAvailable-1.0

#' # Earth-Moon Lyapunov Orbit Transfer
#'
#' Transfer from an L1 Lyapunov orbit to the L2 orbit at the same energy, using as little control
#' as the connection allows. In this example we define the state and control structs, write the
#' circular restricted three-body equations of motion and their partials, and create a
#' `CollocationPhase` using the `HermiteSimpson` transcription. Then we declare the optimization
#' variables with `Vary`, pin both ends with `Constraint`s, and minimize the integrated squared
#' control with `Objective`.
#'
#' The mass ratio, the libration points and both Lyapunov states are JPL's, and AstroRoutines
#' reproduces the published libration points and Jacobi constants from them. Two orbits of equal
#' energy are joined by their invariant manifolds, so a ballistic connection exists and the cost of
#' the transfer approaches zero.
#'
#' References:
#'
#' - JPL Solar System Dynamics, periodic orbit catalogue,
#'   <https://ssd-api.jpl.nasa.gov/doc/periodic_orbits.html>
#' - Koon, Lo, Marsden and Ross, *Chaos* 10(2), 2000.

using Epicycle
using LinearAlgebra

#' ## Configuration
#'
#' Define the Earth-Moon mass ratio and the two orbit states, in normalised units where the Earth
#' and Moon are one unit apart. Each state sits on the x-axis, moving in the +y direction.

# Define the published values, from JPL's periodic orbit catalogue
const MU = 1.215058560962404e-02

const L1_X, L1_VY = 8.1596252146384562e-01, 2.0722124749217649e-01
const L2_X, L2_VY = 1.1182825695532028e+00, 1.8601928389638619e-01
const L1_PERIOD, L2_PERIOD = 2.8447547942812985e+00, 3.4205694071950448e+00

const S_L1 = [L1_X, 0.0, 0.0, 0.0, L1_VY, 0.0]
const S_L2 = [L2_X, 0.0, 0.0, 0.0, L2_VY, 0.0]

# Reproduce the published benchmark values from the mass ratio alone
println("published benchmark, reproduced by AstroRoutines")
println("  L1 point          ", round(libration_point(MU, :L1)[1], digits = 12),
        "   JPL 0.836915125772357")
println("  L2 point          ", round(libration_point(MU, :L2)[1], digits = 12),
        "   JPL 1.15568216544488")
println("  C of the L1 orbit ", round(jacobi_constant(S_L1, MU), digits = 12),
        "   JPL 3.15001683280912")
println("  C of the L2 orbit ", round(jacobi_constant(S_L2, MU), digits = 12),
        "   JPL 3.15000013081292")

# Create the model, whose only constant is the mass ratio
struct CRModel
    mu::Float64
end

model = CRModel(MU)

# Define the state and the control
struct CRState{T} <: AbstractState
    x::T
    y::T
    z::T
    vx::T
    vy::T
    vz::T
end

struct CRControl{T} <: AbstractControl
    ux::T
    uy::T
    uz::T
end

#' ## Equations of motion
#'
#' The ballistic acceleration is AstroRoutines', and the control adds to it. The same package
#' supplies the Jacobian, so the state partial is one call.

# Define the equations of motion
function cr_dynamics!(dy, y::CRState, u::CRControl, p, t, model)
    a = cr3bp_accel([y.x, y.y, y.z, y.vx, y.vy, y.vz], model.mu)
    dy[1] = y.vx
    dy[2] = y.vy
    dy[3] = y.vz
    dy[4] = a[1] + u.ux
    dy[5] = a[2] + u.uy
    dy[6] = a[3] + u.uz
end

# Declare the partial of the equations of motion with respect to the state
@partial(cr_dynamics!, state) do dF, y, u, p, t, model
    dF .= cr3bp_jacobian([y.x, y.y, y.z, y.vx, y.vy, y.vz], model.mu)
end

# Declare the partial of the equations of motion with respect to the control
@partial(cr_dynamics!, control) do dF, y, u, p, t, model
    dF[4, 1] = 1.0
    dF[5, 2] = 1.0
    dF[6, 3] = 1.0
end

#' ## Configure the optimal control problem
#'
#' The transfer leaves a point on the L1 orbit and arrives at a point on the L2 orbit. The final
#' time is free, with both endpoints pinned, so shrinking the duration cannot make the cost vanish.

# Create the collocation phase using the HermiteSimpson transcription
const TF = (L1_PERIOD + L2_PERIOD) / 2

phase = CollocationPhase(name = :lyapunov_transfer,
                         transcription = HermiteSimpson(n_steps = 50),
                         dynamics = cr_dynamics!,
                         model = model,
                         state = CRState,
                         control = CRControl,
                         tspan = (0.0, TF))

# Declare the state, the control and the final time as optimization variables
Vary(state, phase;
     guess = hcat(S_L1, S_L2),
     lower_bound = [0.70, -0.45, -0.05, -3.0, -3.0, -0.5],
     upper_bound = [1.30, 0.45, 0.05, 3.0, 3.0, 0.5])

Vary(control, phase;
     guess = zeros(3, 2),
     lower_bound = [-0.5, -0.5, -0.1],
     upper_bound = [0.5, 0.5, 0.1])

Vary(final_time, phase;
     guess = TF,
     lower_bound = 0.5 * TF,
     upper_bound = 4.0 * TF)

# Define the state at each end and the integrated control effort
posvel(c) = (y = state(c); [y.x, y.y, y.z, y.vx, y.vy, y.vz])

effort(c) = control(c).ux^2 + control(c).uy^2 + control(c).uz^2

# Declare the partials of the constraint and objective functions
@partial(posvel, state) do c
    Matrix(1.0I, 6, 6)
end

@partial(effort, control) do c
    [2 * control(c).ux  2 * control(c).uy  2 * control(c).uz]
end

# Constrain departure from the L1 orbit and arrival on the L2 orbit
Constraint(posvel, phase; equals = S_L1, at = Initial())
Constraint(posvel, phase; equals = S_L2, at = Final())

# Spend as little control as the connection allows
Objective(effort, phase; sense = Min(), at = Path())

#' ## Solve the optimal control problem
#'
#' The energy the control supplied is the Jacobi constant's change across the transfer, which a
#' ballistic arc would leave at zero.

# Check the declared partials against finite differences
check_partials(phase)

# Solve the transfer problem
result = solve!(Sequence(phase); method = Optimize(print_level = 5))

# Report the transfer against the published orbits
y0 = get_initial_state(phase)
yf = get_final_state(phase)
s0 = [y0.x, y0.y, y0.z, y0.vx, y0.vy, y0.vz]
sf = [yf.x, yf.y, yf.z, yf.vx, yf.vy, yf.vz]

println()
println("transfer")
println("  status             ", result.info)
println("  transfer time      ", round(get_final_time(phase), digits = 6), " TU  (",
        round(get_final_time(phase) * 4.342, digits = 3), " days)")
println("  departure error    ", round(maximum(abs.(s0 .- S_L1)), sigdigits = 3))
println("  arrival error      ", round(maximum(abs.(sf .- S_L2)), sigdigits = 3))
println("  control effort     ", round(result.objective, sigdigits = 8))
println("  C at departure     ", round(jacobi_constant(s0, MU), digits = 9))
println("  C at arrival       ", round(jacobi_constant(sf, MU), digits = 9))
println("  energy supplied    ", round(jacobi_constant(sf, MU) - jacobi_constant(s0, MU),
                                       sigdigits = 3))
